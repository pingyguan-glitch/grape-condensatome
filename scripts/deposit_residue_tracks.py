#!/usr/bin/env python3
"""Deposit per-residue metapredict disorder + AlphaFold pLDDT tracks.

Writes a sidecar SQLite (default ``data/tracks/grape_residue_tracks.sqlite``)
with uint8-quantized trajectories so the main atlas DB stays under GitHub limits.
Results are offline-first for NAR reproducibility; AFDB/metapredict are only
needed at deposit time.

Examples::

    # Grape: disorder (all scored IDs) + pLDDT (AFDB models)
    python3 scripts/deposit_residue_tracks.py --species grape --disorder --plddt \\
        --workers 16

    # Priority subset first (gold ∪ Tier A∪B), then resume full
    python3 scripts/deposit_residue_tracks.py --species grape --disorder --tier-ab-first
    python3 scripts/deposit_residue_tracks.py --species grape --plddt --workers 24

    # Fruit-crop path (documented subset; same schema)
    python3 scripts/deposit_residue_tracks.py --species soybean --tier A --limit 500 \\
        --plddt --disorder --tracks-db data/tracks/soybean_residue_tracks.sqlite
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sqlite3
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from tracks_lib import (  # noqa: E402
    coverage,
    open_tracks_db,
    quantize_plddt,
    quantize_unit_interval,
    set_meta,
    upsert_track,
)
from _http_util import urlopen_bytes, urlopen_text  # noqa: E402

DATA = ROOT / "data"
DEFAULT_GRAPE_DB = DATA / "grape_condensatome.sqlite"
DEFAULT_GRAPE_FASTA = DATA / "UP000009183.fasta.gz"
DEFAULT_TRACKS = DATA / "tracks" / "grape_residue_tracks.sqlite"
AF_CACHE = DATA / "alphafold"
API = "https://alphafold.ebi.ac.uk/api/prediction/{pid}"
BROWSER_UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)
CONF_VERSIONS = (6, 5, 4, 3)

SPECIES_DB = {
    "grape": DEFAULT_GRAPE_DB,
    "soybean": DATA / "second_crop" / "soybean_condensatome.sqlite",
    "tomato": DATA / "second_crop" / "tomato_condensatome.sqlite",
}


def iter_fasta(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        pid = None
        buf: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if pid is not None:
                    yield pid, "".join(buf)
                header = line[1:].split()[0]
                parts = header.split("|")
                pid = parts[1] if len(parts) > 1 else header
                buf = []
            else:
                buf.append(line.strip())
        if pid is not None:
            yield pid, "".join(buf)


def _meta_clean(sequence: str) -> str:
    return "".join(ch for ch in sequence.upper() if ch in "ACDEFGHIKLMNPQRSTVWY")


def predict_disorder_trajectory(seq: str) -> list[float]:
    import metapredict as meta

    clean = _meta_clean(seq)
    if not clean:
        return []
    scores = meta.predict_disorder(clean)
    return [float(x) for x in scores]


def predict_disorder_batch(
    records: list[tuple[str, str]],
    *,
    python_exe: str | None = None,
) -> dict[str, list[float]]:
    """In-process metapredict when available; else subprocess with METAPREDICT_PYTHON."""
    if not python_exe:
        try:
            import metapredict  # noqa: F401

            out: dict[str, list[float]] = {}
            for pid, seq in records:
                out[pid] = predict_disorder_trajectory(seq)
            return out
        except Exception:
            pass
    # Subprocess path returns trajectories as lists
    import subprocess
    import tempfile

    exe = python_exe or os.environ.get("METAPREDICT_PYTHON") or sys.executable
    with tempfile.TemporaryDirectory(prefix="meta_tracks_") as td:
        td_path = Path(td)
        fa = td_path / "in.fasta"
        out_json = td_path / "out.json"
        with fa.open("w") as fh:
            for pid, seq in records:
                clean = _meta_clean(seq)
                fh.write(f">{pid}\n")
                for i in range(0, len(clean), 80):
                    fh.write(clean[i : i + 80] + "\n")
        helper = r"""
import json, sys
from pathlib import Path
import metapredict as meta
fa, out = Path(sys.argv[1]), Path(sys.argv[2])
recs = []
pid = None
buf = []
for line in fa.read_text().splitlines():
    if line.startswith(">"):
        if pid is not None:
            recs.append((pid, "".join(buf)))
        pid = line[1:].split()[0]
        buf = []
    else:
        buf.append(line.strip())
if pid is not None:
    recs.append((pid, "".join(buf)))
result = {}
for p, s in recs:
    if not s:
        result[p] = []
        continue
    result[p] = [float(x) for x in meta.predict_disorder(s)]
out.write_text(json.dumps(result))
"""
        script = td_path / "run_meta_traj.py"
        script.write_text(helper)
        proc = subprocess.run(
            [exe, str(script), str(fa), str(out_json)],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0 or not out_json.exists():
            raise RuntimeError(
                "metapredict trajectory subprocess failed "
                f"(exe={exe}, code={proc.returncode}):\n{proc.stderr[-2000:]}"
            )
        return {k: [float(x) for x in v] for k, v in json.loads(out_json.read_text()).items()}


def _unwrap_prediction(data):
    if isinstance(data, list):
        return data[0] if data else None
    return data if isinstance(data, dict) else None


def fetch_prediction(pid: str, cache: Path) -> dict | None:
    path = cache / f"{pid}.prediction.json"
    if path.exists() and path.stat().st_size >= 2:
        try:
            cached = _unwrap_prediction(json.loads(path.read_text()))
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
    seen: set[str] = set()
    out: list[str] = []
    for u in urls:
        if u and u not in seen:
            seen.add(u)
            out.append(u)
    return out


def fetch_plddt_array(
    pid: str,
    *,
    pred: dict | None,
    model_id: str | None,
    cache: Path,
) -> tuple[list[float] | None, str | None]:
    """Return (plddt_list, source_url_or_cache) or (None, None)."""
    cpath = cache / f"{pid}.confidence.json"
    if cpath.exists() and cpath.stat().st_size > 20:
        try:
            payload = json.loads(cpath.read_text())
            scores = payload.get("confidenceScore") or payload.get("pLDDT") or payload.get("plddt")
            if scores:
                return [float(x) for x in scores], f"cache:{cpath.name}"
        except json.JSONDecodeError:
            pass
    for url in _confidence_urls(pred, pid, model_id):
        try:
            raw = urlopen_bytes(url, timeout=60, accept="application/json", ua=BROWSER_UA)
            cpath.write_bytes(raw)
            payload = json.loads(raw.decode("utf-8"))
            scores = payload.get("confidenceScore") or payload.get("pLDDT") or payload.get("plddt")
            if scores:
                return [float(x) for x in scores], url
        except Exception:
            continue
    return None, None


def select_ids(
    conn: sqlite3.Connection,
    *,
    tier: str | None,
    tier_ab_first: bool,
    limit: int,
    gold_extra: set[str],
) -> list[str]:
    if tier_ab_first:
        rows = list(
            conn.execute(
                "SELECT protein_id FROM proteins WHERE tier IN ('A','B') ORDER BY score DESC"
            )
        )
        ids = [r[0] for r in rows]
        for pid in sorted(gold_extra):
            if pid not in ids:
                ids.insert(0, pid)
        if limit:
            ids = ids[:limit]
        return ids
    q = "SELECT protein_id FROM proteins"
    params: list = []
    if tier:
        q += " WHERE tier = ?"
        params.append(tier.upper())
    q += " ORDER BY score DESC"
    ids = [r[0] for r in conn.execute(q, params)]
    for pid in sorted(gold_extra):
        if pid not in ids:
            ids.append(pid)
    if limit:
        ids = ids[:limit]
    return ids


def load_gold_ids() -> set[str]:
    extra = {"F6HET4"}
    gold = DATA / "benchmarks" / "gold_standard_reconstructed.csv"
    if gold.exists():
        import csv

        for row in csv.DictReader(gold.open()):
            if row.get("protein_id"):
                extra.add(row["protein_id"])
    curated = DATA / "curated" / "llps_experimental_gold.tsv"
    if curated.exists():
        import csv

        for row in csv.DictReader(curated.open(), delimiter="\t"):
            if row.get("protein_id"):
                extra.add(row["protein_id"])
    return extra


def resolve_species_paths(species: str, fruit_slug: str | None) -> tuple[Path, Path | None, Path]:
    species = species.lower().strip()
    if fruit_slug:
        slug = fruit_slug.lower().strip()
        db = DATA / "fruit_crops" / f"{slug}_condensatome.sqlite"
        tracks = DATA / "tracks" / f"{slug}_residue_tracks.sqlite"
        return db, None, tracks
    if species == "grape":
        return DEFAULT_GRAPE_DB, DEFAULT_GRAPE_FASTA, DEFAULT_TRACKS
    if species in SPECIES_DB:
        db = SPECIES_DB[species]
        tracks = DATA / "tracks" / f"{species}_residue_tracks.sqlite"
        return db, None, tracks
    raise SystemExit(f"Unknown species={species!r}; use grape/soybean/tomato or --fruit-slug")


def deposit_disorder(
    *,
    atlas: sqlite3.Connection,
    tracks: sqlite3.Connection,
    ids: list[str],
    fasta: Path | None,
    batch_size: int,
    python_exe: str | None,
    force: bool,
) -> dict:
    have = set()
    if not force:
        have = {
            r[0]
            for r in tracks.execute(
                "SELECT protein_id FROM residue_tracks WHERE disorder_u8 IS NOT NULL"
            )
        }
    pending_ids = [pid for pid in ids if pid not in have]
    print(f"Disorder pending: {len(pending_ids)} (skip {len(ids) - len(pending_ids)})", flush=True)

    seq_map: dict[str, str] = {}
    if fasta and fasta.exists():
        want = set(pending_ids)
        for pid, seq in iter_fasta(fasta):
            if pid in want:
                seq_map[pid] = seq
                if len(seq_map) >= len(want):
                    break
    else:
        # Fruit / second-crop DBs store no sequence — pull from UniProt FASTA stream is heavy;
        # require --fasta for disorder when not grape.
        raise SystemExit(
            "Disorder deposit needs --fasta (grape default UP000009183.fasta.gz). "
            "For fruit crops, download the proteome FASTA first."
        )

    t0 = time.time()
    n_ok = 0
    n_fail = 0
    for i in range(0, len(pending_ids), batch_size):
        batch_ids = pending_ids[i : i + batch_size]
        records = [(pid, seq_map[pid]) for pid in batch_ids if pid in seq_map]
        missing = [pid for pid in batch_ids if pid not in seq_map]
        if missing:
            n_fail += len(missing)
            print(f"  missing FASTA for {len(missing)} ids (e.g. {missing[:3]})", flush=True)
        if not records:
            continue
        print(
            f"  disorder batch {i // batch_size + 1}: {len(records)} seqs …",
            flush=True,
        )
        try:
            trajs = predict_disorder_batch(records, python_exe=python_exe)
        except Exception as e:
            print(f"  batch failed: {e}", flush=True)
            n_fail += len(records)
            continue
        for pid, seq in records:
            arr = trajs.get(pid) or []
            if not arr:
                n_fail += 1
                continue
            mean = sum(arr) / len(arr)
            upsert_track(
                tracks,
                protein_id=pid,
                length=len(arr),
                disorder_u8=quantize_unit_interval(arr),
                disorder_source="metapredict",
                disorder_mean=round(mean, 6),
            )
            n_ok += 1
        tracks.commit()
        elapsed = time.time() - t0
        done = min(i + batch_size, len(pending_ids))
        print(
            f"    ok={n_ok} fail={n_fail} ({done}/{len(pending_ids)}; "
            f"{done / max(elapsed, 1e-6):.1f} seq/s)",
            flush=True,
        )
    summary = {
        "n_requested": len(ids),
        "n_pending": len(pending_ids),
        "n_ok": n_ok,
        "n_fail": n_fail,
        "source": "metapredict",
        "elapsed_sec": round(time.time() - t0, 1),
    }
    set_meta(tracks, "disorder_deposit", summary)
    # Mirror mean into main atlas when column exists (grape).
    cols = {r[1] for r in atlas.execute("PRAGMA table_info(proteins)")}
    if "frac_disorder_metapredict" in cols and n_ok:
        updated = 0
        for (pid, mean) in tracks.execute(
            "SELECT protein_id, disorder_mean FROM residue_tracks WHERE disorder_mean IS NOT NULL"
        ):
            atlas.execute(
                "UPDATE proteins SET frac_disorder_metapredict = ? WHERE protein_id = ?",
                (mean, pid),
            )
            updated += 1
        atlas.commit()
        summary["atlas_frac_updated"] = updated
    return summary


def deposit_plddt(
    *,
    atlas: sqlite3.Connection,
    tracks: sqlite3.Connection,
    ids: list[str],
    workers: int,
    force: bool,
    cache_dir: Path,
) -> dict:
    cache_dir.mkdir(parents=True, exist_ok=True)
    have = set()
    if not force:
        have = {
            r[0]
            for r in tracks.execute(
                "SELECT protein_id FROM residue_tracks WHERE plddt_u8 IS NOT NULL"
            )
        }
    pending = [pid for pid in ids if pid not in have]
    model_by_id = {
        r[0]: r[1]
        for r in atlas.execute("SELECT protein_id, af_model_id FROM proteins")
        if True
    }
    # Prefer IDs that already have AF summaries when column present
    cols = {r[1] for r in atlas.execute("PRAGMA table_info(proteins)")}
    if "af_mean_plddt" in cols:
        with_af = {
            r[0]
            for r in atlas.execute(
                "SELECT protein_id FROM proteins WHERE af_mean_plddt IS NOT NULL"
            )
        }
        # Still attempt gold / selected IDs without a mean (may resolve via prediction API)
        pending = [pid for pid in pending if pid in with_af or force]
        pending_set = set(pending)
        for pid in ids:
            if pid in have or pid in pending_set:
                continue
            if pid in with_af:
                pending.append(pid)
                pending_set.add(pid)

    print(f"pLDDT pending: {len(pending)} (skip {len(ids) - len(pending)})", flush=True)
    t0 = time.time()
    n_ok = 0
    n_fail = 0

    def job(pid: str):
        pred = None
        mid = model_by_id.get(pid)
        try:
            pred = fetch_prediction(pid, cache_dir)
        except Exception:
            pred = None
        if pred and not mid:
            mid = pred.get("modelEntityId") or pred.get("entryId")
        arr, src = fetch_plddt_array(pid, pred=pred, model_id=mid, cache=cache_dir)
        return pid, arr, src, mid, pred

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futs = {pool.submit(job, pid): pid for pid in pending}
        done = 0
        for fut in as_completed(futs):
            done += 1
            pid = futs[fut]
            try:
                pid, arr, src, mid, pred = fut.result()
            except Exception as e:
                n_fail += 1
                if done % 200 == 0 or done == len(pending):
                    print(f"  pLDDT {done}/{len(pending)} ok={n_ok} fail={n_fail} ({e})", flush=True)
                continue
            if not arr:
                n_fail += 1
            else:
                mean = sum(arr) / len(arr)
                upsert_track(
                    tracks,
                    protein_id=pid,
                    length=len(arr),
                    plddt_u8=quantize_plddt(arr),
                    plddt_source=src or "alphafold_db",
                    plddt_mean=round(mean, 4),
                )
                n_ok += 1
                # Refresh AF summary columns when present and empty
                if pred and "af_mean_plddt" in cols:
                    row = atlas.execute(
                        "SELECT af_mean_plddt FROM proteins WHERE protein_id = ?", (pid,)
                    ).fetchone()
                    if row and row[0] is None:
                        atlas.execute(
                            """
                            UPDATE proteins SET
                              af_mean_plddt = ?,
                              af_frac_very_low = ?,
                              af_frac_low = ?,
                              af_frac_confident = ?,
                              af_frac_very_high = ?,
                              af_model_id = ?
                            WHERE protein_id = ?
                            """,
                            (
                                float(pred.get("globalMetricValue") or mean),
                                float(pred.get("fractionPlddtVeryLow") or 0.0),
                                float(pred.get("fractionPlddtLow") or 0.0),
                                float(pred.get("fractionPlddtConfident") or 0.0),
                                float(pred.get("fractionPlddtVeryHigh") or 0.0),
                                mid,
                                pid,
                            ),
                        )
            if done % 200 == 0 or done == len(pending):
                tracks.commit()
                atlas.commit()
                elapsed = time.time() - t0
                print(
                    f"  pLDDT {done}/{len(pending)} ok={n_ok} fail={n_fail} "
                    f"({done / max(elapsed, 1e-6):.1f}/s)",
                    flush=True,
                )
    tracks.commit()
    atlas.commit()
    summary = {
        "n_requested": len(ids),
        "n_pending": len(pending),
        "n_ok": n_ok,
        "n_fail": n_fail,
        "source": "alphafold.ebi.ac.uk confidence JSON",
        "cache_dir": str(cache_dir),
        "elapsed_sec": round(time.time() - t0, 1),
    }
    set_meta(tracks, "plddt_deposit", summary)
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species", default="grape", help="grape | soybean | tomato")
    ap.add_argument("--fruit-slug", default="", help="Fruit crop slug under data/fruit_crops/")
    ap.add_argument("--db", default="", help="Override atlas SQLite path")
    ap.add_argument("--tracks-db", default="", help="Override tracks SQLite path")
    ap.add_argument("--fasta", default="", help="FASTA for disorder (required for non-grape)")
    ap.add_argument("--disorder", action="store_true", help="Deposit metapredict trajectories")
    ap.add_argument("--plddt", action="store_true", help="Deposit AFDB pLDDT trajectories")
    ap.add_argument("--tier", default="", help="Restrict to tier letter")
    ap.add_argument(
        "--tier-ab-first",
        action="store_true",
        help="Only gold ∪ Tier A∪B (priority subset)",
    )
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--force", action="store_true", help="Recompute even if blob present")
    ap.add_argument("--metapredict-python", default=os.environ.get("METAPREDICT_PYTHON"))
    ap.add_argument("--cache-dir", default=str(AF_CACHE))
    args = ap.parse_args()

    if not args.disorder and not args.plddt:
        ap.error("Specify --disorder and/or --plddt")

    db_path, default_fasta, tracks_path = resolve_species_paths(
        args.species, args.fruit_slug or None
    )
    if args.db:
        db_path = Path(args.db)
    if args.tracks_db:
        tracks_path = Path(args.tracks_db)
    fasta = Path(args.fasta) if args.fasta else default_fasta

    if not db_path.exists():
        raise SystemExit(f"Atlas DB missing: {db_path}")

    atlas = sqlite3.connect(str(db_path), timeout=120.0)
    atlas.execute("PRAGMA busy_timeout=60000")
    tracks = open_tracks_db(tracks_path, writable=True)

    gold = load_gold_ids() if args.species == "grape" and not args.fruit_slug else set()
    ids = select_ids(
        atlas,
        tier=args.tier or None,
        tier_ab_first=args.tier_ab_first,
        limit=args.limit,
        gold_extra=gold,
    )
    print(f"Selected {len(ids)} proteins → {tracks_path}", flush=True)

    report: dict = {
        "species": args.fruit_slug or args.species,
        "atlas_db": str(db_path),
        "tracks_db": str(tracks_path),
        "n_selected": len(ids),
    }
    if args.disorder:
        if not fasta or not Path(fasta).exists():
            raise SystemExit(f"--disorder requires FASTA (missing: {fasta})")
        report["disorder"] = deposit_disorder(
            atlas=atlas,
            tracks=tracks,
            ids=ids,
            fasta=Path(fasta),
            batch_size=args.batch_size,
            python_exe=args.metapredict_python,
            force=args.force,
        )
    if args.plddt:
        report["plddt"] = deposit_plddt(
            atlas=atlas,
            tracks=tracks,
            ids=ids,
            workers=args.workers,
            force=args.force,
            cache_dir=Path(args.cache_dir),
        )

    cov = coverage(tracks)
    report["coverage"] = cov
    set_meta(tracks, "coverage", cov)
    set_meta(tracks, "last_deposit", report)
    tracks.commit()

    # Mirror coverage into atlas meta when table exists
    try:
        atlas.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            ("residue_tracks", json.dumps(report, sort_keys=True)),
        )
        atlas.commit()
    except sqlite3.Error:
        pass

    cov_path = tracks_path.with_suffix(".coverage.json")
    cov_path.write_text(json.dumps(report, indent=2, sort_keys=True))
    print(json.dumps(report, indent=2), flush=True)
    tracks.close()
    atlas.close()


if __name__ == "__main__":
    main()
