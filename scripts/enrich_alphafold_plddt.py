#!/usr/bin/env python3
"""Join AlphaFold DB pLDDT summaries (+ optional disordered segments) into SQLite.

Uses the public AlphaFold prediction API for mean pLDDT / fraction bins, and
(when reachable) the per-residue confidence JSON for contiguous low-pLDDT
segments. Results are cached under data/alphafold/ so runs are resumable.

Full proteome (~26k) is feasible but network-bound; default scope is the
scored proteome with --limit for smoke tests. Segments default to gold +
anchors + Tier A (pass --segments-all to attempt every fetched ID).

Example::

    python3 scripts/enrich_alphafold_plddt.py
    python3 scripts/enrich_alphafold_plddt.py --limit 200 --with-segments
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
    # confidence JSON uses confidenceScore; some dumps use pLDDT
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
            # size-2 "[]" means previously confirmed missing
            if path.stat().st_size <= 5 and cached is None:
                return None
            if cached is not None:
                return cached
            if path.stat().st_size > 20:
                return None
        except json.JSONDecodeError:
            pass
    try:
        raw = urlopen_text(API.format(pid=pid), timeout=45, accept="application/json")
    except Exception as e:
        code = getattr(e, "code", None)
        if code == 404:
            path.write_text("[]")
            return None
        raise
    path.write_text(raw)
    return _unwrap_prediction(json.loads(raw))


def fetch_segments(pred: dict, pid: str, cache: Path, thr: float) -> tuple[list[dict], int]:
    url = pred.get("plddtDocUrl") or ""
    if not url:
        # fall back to conventional path
        mid = pred.get("modelEntityId") or f"AF-{pid}-F1"
        ver = pred.get("latestVersion") or 4
        url = f"https://alphafold.ebi.ac.uk/files/{mid}-confidence_v{ver}.json"
    cpath = cache / f"{pid}.confidence.json"
    if cpath.exists() and cpath.stat().st_size > 20:
        try:
            payload = json.loads(cpath.read_text())
            return segments_from_confidence(payload, thr)
        except json.JSONDecodeError:
            pass
    try:
        raw = urlopen_bytes(url, timeout=60, accept="application/json")
        cpath.write_bytes(raw)
        payload = json.loads(raw.decode("utf-8"))
        return segments_from_confidence(payload, thr)
    except Exception:
        return [], 0


def summarize_row(pred: dict) -> dict:
    return {
        "af_mean_plddt": float(pred.get("globalMetricValue") or 0.0) or None,
        "af_frac_very_low": float(pred.get("fractionPlddtVeryLow") or 0.0),
        "af_frac_low": float(pred.get("fractionPlddtLow") or 0.0),
        "af_frac_confident": float(pred.get("fractionPlddtConfident") or 0.0),
        "af_frac_very_high": float(pred.get("fractionPlddtVeryHigh") or 0.0),
        "af_model_id": pred.get("modelEntityId") or pred.get("entryId"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--cache-dir", default=str(CACHE_DIR))
    ap.add_argument("--limit", type=int, default=0, help="Max proteins (0 = all in DB)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--with-segments", action="store_true", help="Fetch confidence JSON segments")
    ap.add_argument(
        "--segments-tier",
        default="A",
        help="When --with-segments: also include this tier (default A); use ALL for every ID",
    )
    ap.add_argument("--segment-threshold", type=float, default=50.0)
    ap.add_argument("--sleep", type=float, default=0.0, help="Extra delay between jobs")
    args = ap.parse_args()

    cache = Path(args.cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(args.db))
    ensure_columns(conn)

    # Prefer proteome IDs; always include gold + F6HET4
    ids = [r[0] for r in conn.execute("SELECT protein_id FROM proteins ORDER BY score DESC")]
    if args.limit:
        ids = ids[: args.limit]
    gold = DATA / "benchmarks" / "gold_standard_reconstructed.csv"
    extra = {"F6HET4"}
    if gold.exists():
        import csv

        for row in csv.DictReader(gold.open()):
            extra.add(row["protein_id"])
    for pid in sorted(extra):
        if pid not in ids:
            ids.append(pid)

    segment_ids: set[str] = set(extra)
    if args.with_segments:
        if args.segments_tier.upper() == "ALL":
            segment_ids.update(ids)
        else:
            tier = args.segments_tier.upper()
            for (pid,) in conn.execute(
                "SELECT protein_id FROM proteins WHERE tier = ?", (tier,)
            ):
                segment_ids.add(pid)

    print(f"Fetching AlphaFold summaries for {len(ids)} IDs (workers={args.workers}) …", flush=True)
    ok = fail = updated = 0

    def upsert(pid: str, out: dict) -> bool:
        exists = conn.execute(
            "SELECT 1 FROM proteins WHERE protein_id = ?", (pid,)
        ).fetchone()
        if not exists:
            return False
        conn.execute(
            """
            UPDATE proteins SET
                af_mean_plddt = ?,
                af_frac_very_low = ?,
                af_frac_low = ?,
                af_frac_confident = ?,
                af_frac_very_high = ?,
                af_model_id = ?,
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
            pred = fetch_prediction(pid, cache)
            if not pred:
                return pid, None, "empty"
            out = summarize_row(pred)
            if args.with_segments and pid in segment_ids:
                segs, n_low = fetch_segments(pred, pid, cache, args.segment_threshold)
                out["af_disorder_segments_json"] = json.dumps(segs)
                out["af_n_disorder_residues"] = n_low
            return pid, out, None
        except Exception as e:
            return pid, None, str(e)

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = [ex.submit(job, pid) for pid in ids]
        for i, fut in enumerate(as_completed(futs), 1):
            pid, out, err = fut.result()
            if out:
                ok += 1
                if upsert(pid, out):
                    updated += 1
            else:
                fail += 1
            if i % 200 == 0 or i == len(ids):
                conn.commit()
                print(f"  progress {i}/{len(ids)} ok={ok} fail={fail} db={updated}", flush=True)

    summary = {
        "n_requested": len(ids),
        "n_ok": ok,
        "n_fail": fail,
        "n_updated_db": updated,
        "with_segments": bool(args.with_segments),
        "segment_ids": len(segment_ids) if args.with_segments else 0,
        "cache_dir": str(cache),
    }
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("alphafold_plddt", json.dumps(summary)),
    )
    conn.commit()
    n_af = conn.execute(
        "SELECT COUNT(*) FROM proteins WHERE af_mean_plddt IS NOT NULL"
    ).fetchone()[0]
    print(json.dumps({**summary, "n_with_af_mean_plddt": n_af}, indent=2))
    for row in conn.execute(
        "SELECT protein_id, af_mean_plddt, af_frac_very_low, af_n_disorder_residues "
        "FROM proteins WHERE protein_id IN ('F6HET4') OR af_mean_plddt IS NOT NULL "
        "ORDER BY CASE protein_id WHEN 'F6HET4' THEN 0 ELSE 1 END LIMIT 5"
    ):
        print(" ", row)
    conn.close()


if __name__ == "__main__":
    main()
