#!/usr/bin/env python3
"""Join AlphaFold DB pLDDT summaries (+ low-pLDDT segments) into SQLite.

Uses the public AlphaFold prediction API for mean pLDDT / fraction bins, and
per-residue confidence JSON for contiguous low-pLDDT segments. Results are
cached under data/alphafold/ so runs are resumable.

v0.7 default for segments is the **scored proteome** (``--segments-tier ALL``):
every protein with an AFDB model gets confidence segments when reachable.
File downloads need a **browser-like User-Agent** (AFDB returns 403 otherwise).

Example::

    python3 scripts/enrich_alphafold_plddt.py --workers 16
    python3 scripts/enrich_alphafold_plddt.py --segments-only --segments-tier ALL --workers 24
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from _http_util import urlopen_bytes, urlopen_text  # noqa: E402

DATA = ROOT / "data"
DEFAULT_DB = DATA / "grape_condensatome.sqlite"
CACHE_DIR = DATA / "alphafold"
API = "https://alphafold.ebi.ac.uk/api/prediction/{pid}"

# AFDB file CDN rejects non-browser UAs with HTTP 403.
BROWSER_UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)
CONF_VERSIONS = (6, 5, 4, 3)


def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    for col, decl in (
        ("af_mean_plddt", "REAL"),
        ("af_frac_very_low", "REAL"),
        ("af_frac_low", "REAL"),
        ("af_frac_confident", "REAL"),
        ("af_frac_very_high", "REAL"),
        ("af_model_id", "TEXT"),
        ("af_disorder_segments_json", "TEXT"),
        ("af_n_disorder_residues", "INTEGER"),
    ):
        if col not in cols:
            conn.execute(f"ALTER TABLE proteins ADD COLUMN {col} {decl}")


def segments_from_confidence(payload: dict, thr: float = 50.0) -> tuple[list[dict], int]:
    nums = payload.get("residueNumber") or []
    scores = payload.get("confidenceScore") or payload.get("pLDDT") or payload.get("plddt")
    if not scores:
        return [], 0
    low_idx = [int(nums[i]) if i < len(nums) else i + 1 for i, v in enumerate(scores) if float(v) < thr]
    if not low_idx:
        return [], 0
    segs: list[dict] = []
    s = e = low_idx[0]
    for x in low_idx[1:]:
        if x == e + 1:
            e = x
        else:
            segs.append({"start": s, "end": e, "length": e - s + 1})
            s = e = x
    segs.append({"start": s, "end": e, "length": e - s + 1})
    return segs, len(low_idx)


def _unwrap_prediction(data):
    if isinstance(data, list):
        return data[0] if data else None
    return data if isinstance(data, dict) else None


def fetch_prediction(pid: str, cache: Path) -> dict | None:
    path = cache / f"{pid}.prediction.json"
    if path.exists() and path.stat().st_size >= 2:
        try:
            cached = _unwrap_prediction(json.loads(path.read_text()))
            if path.stat().st_size <= 5 and cached is None:
                return None
            if cached is not None:
                return cached
            if path.stat().st_size > 20:
                return None
        except json.JSONDecodeError:
            pass
    try:
        raw = urlopen_text(
            API.format(pid=pid),
            timeout=45,
            accept="application/json",
            ua=BROWSER_UA,
        )
    except Exception as e:
        code = getattr(e, "code", None)
        if code == 404:
            path.write_text("[]")
            return None
        raise
    path.write_text(raw)
    return _unwrap_prediction(json.loads(raw))


def _confidence_urls(pred: dict | None, pid: str, model_id: str | None) -> list[str]:
    urls: list[str] = []
    if pred:
        url = pred.get("plddtDocUrl") or ""
        if url:
            urls.append(url)
        mid = pred.get("modelEntityId") or model_id or f"AF-{pid}-F1"
        ver = pred.get("latestVersion")
        if ver:
            urls.append(f"https://alphafold.ebi.ac.uk/files/{mid}-confidence_v{ver}.json")
    mid = model_id or (pred or {}).get("modelEntityId") or f"AF-{pid}-F1"
    for ver in CONF_VERSIONS:
        urls.append(f"https://alphafold.ebi.ac.uk/files/{mid}-confidence_v{ver}.json")
    # de-dupe preserving order
    seen: set[str] = set()
    out: list[str] = []
    for u in urls:
        if u and u not in seen:
            seen.add(u)
            out.append(u)
    return out


def fetch_segments(
    pred: dict | None,
    pid: str,
    cache: Path,
    thr: float,
    model_id: str | None = None,
) -> tuple[list[dict], int, bool]:
    """Return (segments, n_low, ok). ok=False means download failed entirely."""
    cpath = cache / f"{pid}.confidence.json"
    if cpath.exists() and cpath.stat().st_size > 20:
        try:
            payload = json.loads(cpath.read_text())
            segs, n = segments_from_confidence(payload, thr)
            return segs, n, True
        except json.JSONDecodeError:
            pass
    last_err: Exception | None = None
    for url in _confidence_urls(pred, pid, model_id):
        try:
            raw = urlopen_bytes(url, timeout=60, accept="application/json", ua=BROWSER_UA)
            cpath.write_bytes(raw)
            payload = json.loads(raw.decode("utf-8"))
            segs, n = segments_from_confidence(payload, thr)
            return segs, n, True
        except Exception as e:
            code = getattr(e, "code", None)
            last_err = e
            if code == 404:
                continue
            # 403 / transient — try next version then give up for this ID
            continue
    if last_err is not None:
        return [], 0, False
    return [], 0, False


def summarize_row(pred: dict) -> dict:
    return {
        "af_mean_plddt": float(pred.get("globalMetricValue") or 0.0) or None,
        "af_frac_very_low": float(pred.get("fractionPlddtVeryLow") or 0.0),
        "af_frac_low": float(pred.get("fractionPlddtLow") or 0.0),
        "af_frac_confident": float(pred.get("fractionPlddtConfident") or 0.0),
        "af_frac_very_high": float(pred.get("fractionPlddtVeryHigh") or 0.0),
        "af_model_id": pred.get("modelEntityId") or pred.get("entryId"),
    }


def coverage_counts(conn: sqlite3.Connection) -> dict:
    row = conn.execute(
        """
        SELECT
          COUNT(*) AS n_proteins,
          SUM(af_mean_plddt IS NOT NULL) AS n_with_mean,
          SUM(af_disorder_segments_json IS NOT NULL) AS n_with_segments_field,
          SUM(
            af_disorder_segments_json IS NOT NULL
            AND af_disorder_segments_json != ''
            AND af_disorder_segments_json != '[]'
          ) AS n_with_nonempty_segments,
          SUM(IFNULL(af_n_disorder_residues, 0) > 0) AS n_with_low_residues
        FROM proteins
        """
    ).fetchone()
    keys = (
        "n_proteins",
        "n_with_mean",
        "n_with_segments_field",
        "n_with_nonempty_segments",
        "n_with_low_residues",
    )
    return dict(zip(keys, row))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--cache-dir", default=str(CACHE_DIR))
    ap.add_argument("--limit", type=int, default=0, help="Max proteins (0 = all in DB)")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--with-segments", action="store_true", help="Fetch confidence JSON segments")
    ap.add_argument(
        "--segments-only",
        action="store_true",
        help="Only fetch/join low-pLDDT segments (skip API summary refresh when model_id known)",
    )
    ap.add_argument(
        "--segments-tier",
        default="ALL",
        help="When fetching segments: tier letter, or ALL for every ID with an AF model (default ALL)",
    )
    ap.add_argument(
        "--force-segments",
        action="store_true",
        help="Re-fetch segments even when af_disorder_segments_json is already set",
    )
    ap.add_argument("--segment-threshold", type=float, default=50.0)
    ap.add_argument("--sleep", type=float, default=0.0, help="Extra delay between jobs")
    args = ap.parse_args()

    want_segments = bool(args.with_segments or args.segments_only)
    if args.segments_only:
        args.with_segments = True
        want_segments = True

    cache = Path(args.cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(args.db), timeout=120.0)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=60000")
    ensure_columns(conn)

    gold = DATA / "benchmarks" / "gold_standard_reconstructed.csv"
    extra = {"F6HET4"}
    if gold.exists():
        import csv

        for row in csv.DictReader(gold.open()):
            extra.add(row["protein_id"])

    if args.segments_only:
        # Prefer IDs that already have AF summaries but lack a segments field.
        q = (
            "SELECT protein_id, af_model_id FROM proteins "
            "WHERE af_model_id IS NOT NULL"
        )
        if not args.force_segments:
            q += " AND af_disorder_segments_json IS NULL"
        if args.segments_tier.upper() != "ALL":
            q += f" AND tier = '{args.segments_tier.upper()}'"
        q += " ORDER BY score DESC"
        rows = list(conn.execute(q))
        ids = [r[0] for r in rows]
        model_by_id = {r[0]: r[1] for r in rows}
        for pid in sorted(extra):
            if pid not in model_by_id:
                r = conn.execute(
                    "SELECT af_model_id FROM proteins WHERE protein_id = ?", (pid,)
                ).fetchone()
                if r and r[0]:
                    ids.append(pid)
                    model_by_id[pid] = r[0]
        segment_ids = set(ids)
    else:
        ids = [r[0] for r in conn.execute("SELECT protein_id FROM proteins ORDER BY score DESC")]
        if args.limit:
            ids = ids[: args.limit]
        for pid in sorted(extra):
            if pid not in ids:
                ids.append(pid)
        model_by_id = {
            r[0]: r[1]
            for r in conn.execute("SELECT protein_id, af_model_id FROM proteins")
        }
        segment_ids = set(extra)
        if want_segments:
            if args.segments_tier.upper() == "ALL":
                segment_ids.update(ids)
            else:
                tier = args.segments_tier.upper()
                for (pid,) in conn.execute(
                    "SELECT protein_id FROM proteins WHERE tier = ?", (tier,)
                ):
                    segment_ids.add(pid)
            if not args.force_segments:
                have = {
                    r[0]
                    for r in conn.execute(
                        "SELECT protein_id FROM proteins WHERE af_disorder_segments_json IS NOT NULL"
                    )
                }
                segment_ids -= have

    if args.limit and args.segments_only:
        ids = ids[: args.limit]
        segment_ids = set(ids)

    print(
        f"AlphaFold: n_ids={len(ids)} segments={want_segments} "
        f"segment_ids={len(segment_ids)} workers={args.workers} "
        f"segments_only={args.segments_only}",
        flush=True,
    )
    ok = fail = updated = seg_ok = seg_fail = 0

    def upsert(pid: str, out: dict, segments_only: bool) -> bool:
        exists = conn.execute(
            "SELECT 1 FROM proteins WHERE protein_id = ?", (pid,)
        ).fetchone()
        if not exists:
            return False
        if segments_only:
            conn.execute(
                """
                UPDATE proteins SET
                    af_disorder_segments_json = ?,
                    af_n_disorder_residues = ?
                WHERE protein_id = ?
                """,
                (
                    out.get("af_disorder_segments_json"),
                    out.get("af_n_disorder_residues"),
                    pid,
                ),
            )
        else:
            conn.execute(
                """
                UPDATE proteins SET
                    af_mean_plddt = COALESCE(?, af_mean_plddt),
                    af_frac_very_low = COALESCE(?, af_frac_very_low),
                    af_frac_low = COALESCE(?, af_frac_low),
                    af_frac_confident = COALESCE(?, af_frac_confident),
                    af_frac_very_high = COALESCE(?, af_frac_very_high),
                    af_model_id = COALESCE(?, af_model_id),
                    af_disorder_segments_json = COALESCE(?, af_disorder_segments_json),
                    af_n_disorder_residues = COALESCE(?, af_n_disorder_residues)
                WHERE protein_id = ?
                """,
                (
                    out.get("af_mean_plddt"),
                    out.get("af_frac_very_low"),
                    out.get("af_frac_low"),
                    out.get("af_frac_confident"),
                    out.get("af_frac_very_high"),
                    out.get("af_model_id"),
                    out.get("af_disorder_segments_json"),
                    out.get("af_n_disorder_residues"),
                    pid,
                ),
            )
        return True

    def job(pid: str) -> tuple[str, dict | None, str | None]:
        if args.sleep:
            time.sleep(args.sleep)
        try:
            mid = model_by_id.get(pid)
            pred = None
            out: dict = {}
            if not args.segments_only or not mid:
                pred = fetch_prediction(pid, cache)
                if not pred and not mid:
                    return pid, None, "empty"
                if pred:
                    out.update(summarize_row(pred))
                    mid = out.get("af_model_id") or mid
            if want_segments and pid in segment_ids:
                segs, n_low, sok = fetch_segments(
                    pred, pid, cache, args.segment_threshold, model_id=mid
                )
                if sok:
                    out["af_disorder_segments_json"] = json.dumps(segs)
                    out["af_n_disorder_residues"] = n_low
                else:
                    return pid, out if out else None, "confidence_fail"
            if not out:
                return pid, None, "empty"
            return pid, out, None
        except Exception as e:
            return pid, None, str(e)

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = [ex.submit(job, pid) for pid in ids]
        for i, fut in enumerate(as_completed(futs), 1):
            pid, out, err = fut.result()
            if out and "af_disorder_segments_json" in out:
                seg_ok += 1
            elif err == "confidence_fail":
                seg_fail += 1
            if out:
                ok += 1
                if upsert(pid, out, segments_only=args.segments_only):
                    updated += 1
            else:
                fail += 1
            if i % 200 == 0 or i == len(ids):
                conn.commit()
                print(
                    f"  progress {i}/{len(ids)} ok={ok} fail={fail} "
                    f"seg_ok={seg_ok} seg_fail={seg_fail} db={updated}",
                    flush=True,
                )

    cov = coverage_counts(conn)
    summary = {
        "n_requested": len(ids),
        "n_ok": ok,
        "n_fail": fail,
        "n_updated_db": updated,
        "with_segments": want_segments,
        "segments_only": bool(args.segments_only),
        "segments_tier": args.segments_tier,
        "segment_ids_targeted": len(segment_ids) if want_segments else 0,
        "n_segments_ok": seg_ok,
        "n_segments_fail": seg_fail,
        "cache_dir": str(cache),
        "user_agent": "browser-like",
        **cov,
    }
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("alphafold_plddt", json.dumps(summary)),
    )
    cov_path = cache / "coverage.json"
    cov_path.write_text(json.dumps(summary, indent=2) + "\n")
    conn.commit()
    print(json.dumps(summary, indent=2))
    for row in conn.execute(
        "SELECT protein_id, af_mean_plddt, af_frac_very_low, af_n_disorder_residues "
        "FROM proteins WHERE protein_id IN ('F6HET4') OR af_mean_plddt IS NOT NULL "
        "ORDER BY CASE protein_id WHEN 'F6HET4' THEN 0 ELSE 1 END LIMIT 5"
    ):
        print(" ", row)
    conn.close()


if __name__ == "__main__":
    main()
