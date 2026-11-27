#!/usr/bin/env python3
"""Merge llphyscore_proteome.shard*.checkpoint.json into the main TSV + SQLite.

Used when scoring was fan-out across parallel shard workers::

    python3 scripts/merge_llphyscore_shards.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sqlite3
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
BENCH = DATA / "benchmarks"
DEFAULT_DB = DATA / "grape_condensatome.sqlite"
DEFAULT_FASTA = DATA / "UP000009183.fasta.gz"
DEFAULT_OUT = BENCH / "llphyscore_proteome.tsv"
DEFAULT_CKPT = BENCH / "llphyscore_proteome.checkpoint.json"


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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--checkpoint", default=str(DEFAULT_CKPT))
    ap.add_argument("--fasta", default=str(DEFAULT_FASTA))
    ap.add_argument("--glob", default="llphyscore_proteome.shard*.checkpoint.json")
    ap.add_argument("--skip-db-join", action="store_true")
    args = ap.parse_args()

    scores: dict[str, float] = {}
    # seed from main if present
    main_ckpt = Path(args.checkpoint)
    if main_ckpt.exists():
        scores.update(
            {k: float(v) for k, v in json.loads(main_ckpt.read_text()).get("scores", {}).items()}
        )
    shard_files = sorted(BENCH.glob(args.glob))
    print(f"Merging {len(shard_files)} shard checkpoints (+ main n={len(scores)})", flush=True)
    for path in shard_files:
        part = json.loads(path.read_text()).get("scores") or {}
        before = len(scores)
        scores.update({k: float(v) for k, v in part.items()})
        print(f"  {path.name}: +{len(scores) - before} → n={len(scores)}", flush=True)

    meta = {
        "n": len(scores),
        "merged_from": [p.name for p in shard_files],
        "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "native": True,
        "model": "human+PDB",
    }
    main_ckpt.write_text(json.dumps({"scores": scores, **meta}, indent=2))

    length_map: dict[str, int] = {}
    if scores and Path(args.fasta).exists():
        for pid, seq in iter_fasta(Path(args.fasta)):
            if pid in scores:
                length_map[pid] = len(seq)

    out = Path(args.out)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["protein_id", "length", "score_llphyscore_percentile"],
            delimiter="\t",
        )
        w.writeheader()
        for pid in sorted(scores):
            w.writerow(
                {
                    "protein_id": pid,
                    "length": length_map.get(pid, ""),
                    "score_llphyscore_percentile": round(scores[pid], 8),
                }
            )

    summary = {
        **meta,
        "out": str(out),
        "mean_percentile": round(sum(scores.values()) / max(1, len(scores)), 4) if scores else None,
    }
    Path(str(out) + ".summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    if not args.skip_db_join and Path(args.db).exists() and scores:
        conn = sqlite3.connect(str(args.db))
        cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
        if "score_llphyscore_percentile" not in cols:
            conn.execute("ALTER TABLE proteins ADD COLUMN score_llphyscore_percentile REAL")
        if "comparator_set" not in cols:
            conn.execute("ALTER TABLE proteins ADD COLUMN comparator_set TEXT")
        updated = 0
        for pid, sc in scores.items():
            cur = conn.execute(
                """
                UPDATE proteins SET
                    score_llphyscore_percentile = ?,
                    comparator_set = COALESCE(comparator_set, 'proteome_native')
                WHERE protein_id = ?
                """,
                (float(sc), pid),
            )
            updated += cur.rowcount
        conn.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            ("llphyscore_proteome", json.dumps(summary)),
        )
        conn.commit()
        print(f"Joined {updated} rows into {args.db}", flush=True)
        conn.close()


if __name__ == "__main__":
    main()
