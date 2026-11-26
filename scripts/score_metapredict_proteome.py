#!/usr/bin/env python3
"""Score the grapevine proteome with metapredict (batched + resumable).

Writes data/benchmarks/metapredict_proteome.tsv and joins fractions into
SQLite proteins.features_json + a dedicated column when present.

Default: full scored proteome FASTA. Use --batch-size / --limit for CPU-bound
machines; progress is checkpointed so re-runs skip completed IDs.

Example::

    python3 scripts/score_metapredict_proteome.py
    python3 scripts/score_metapredict_proteome.py --limit 500 --batch-size 50
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sqlite3
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plantllps.metapredict_bridge import resolve_metapredict_scores  # noqa: E402

DATA = ROOT / "data"
DEFAULT_DB = DATA / "grape_condensatome.sqlite"
DEFAULT_FASTA = DATA / "UP000009183.fasta.gz"
DEFAULT_OUT = DATA / "benchmarks" / "metapredict_proteome.tsv"
DEFAULT_CKPT = DATA / "benchmarks" / "metapredict_proteome.checkpoint.json"


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


def load_checkpoint(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text())
        return {k: float(v) for k, v in (data.get("scores") or {}).items()}
    except (json.JSONDecodeError, TypeError, ValueError):
        return {}


def save_checkpoint(path: Path, scores: dict[str, float], meta: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"scores": scores, **meta}, indent=2))


def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    if "frac_disorder_metapredict" not in cols:
        conn.execute("ALTER TABLE proteins ADD COLUMN frac_disorder_metapredict REAL")


def join_db(conn: sqlite3.Connection, scores: dict[str, float]) -> int:
    ensure_columns(conn)
    updated = 0
    for pid, frac in scores.items():
        row = conn.execute(
            "SELECT features_json FROM proteins WHERE protein_id = ?", (pid,)
        ).fetchone()
        if not row:
            continue
        try:
            feats = json.loads(row[0] or "{}")
        except json.JSONDecodeError:
            feats = {}
        feats["frac_disorder_metapredict"] = float(frac)
        conn.execute(
            "UPDATE proteins SET features_json = ?, frac_disorder_metapredict = ? WHERE protein_id = ?",
            (json.dumps(feats, sort_keys=True), float(frac), pid),
        )
        updated += 1
    return updated


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", default=str(DEFAULT_FASTA))
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--checkpoint", default=str(DEFAULT_CKPT))
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--limit", type=int, default=0, help="Max new sequences to score this run")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument(
        "--metapredict-python",
        default=os.environ.get("METAPREDICT_PYTHON"),
    )
    ap.add_argument("--db-ids-only", action="store_true", help="Only score IDs present in SQLite")
    args = ap.parse_args()

    ckpt_path = Path(args.checkpoint)
    scores = load_checkpoint(ckpt_path)
    print(f"Checkpoint has {len(scores)} scored IDs", flush=True)

    db_ids = None
    if args.db_ids_only and Path(args.db).exists():
        conn = sqlite3.connect(str(args.db))
        db_ids = {r[0] for r in conn.execute("SELECT protein_id FROM proteins")}
        conn.close()
        print(f"Restricting to {len(db_ids)} DB IDs", flush=True)

    pending: list[tuple[str, str]] = []
    for pid, seq in iter_fasta(Path(args.fasta)):
        if pid in scores:
            continue
        if db_ids is not None and pid not in db_ids:
            continue
        pending.append((pid, seq))
        if args.limit and len(pending) >= args.limit:
            break

    print(f"Pending this run: {len(pending)} (batch={args.batch_size})", flush=True)
    t0 = time.time()
    for i in range(0, len(pending), args.batch_size):
        batch = pending[i : i + args.batch_size]
        print(f"  batch {i // args.batch_size + 1}: {len(batch)} seqs …", flush=True)
        part = resolve_metapredict_scores(
            batch,
            threshold=args.threshold,
            python_exe=args.metapredict_python,
        )
        scores.update(part)
        save_checkpoint(
            ckpt_path,
            scores,
            {
                "threshold": args.threshold,
                "n": len(scores),
                "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
        )
        elapsed = time.time() - t0
        done = min(i + args.batch_size, len(pending))
        rate = done / max(elapsed, 1e-6)
        print(
            f"    checkpoint n={len(scores)} ({done}/{len(pending)} this run; {rate:.1f} seq/s)",
            flush=True,
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    # lengths from fasta for written rows
    lengths = {pid: len(seq) for pid, seq in pending}
    # also keep prior
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["protein_id", "length", "frac_disorder_metapredict"],
            delimiter="\t",
        )
        w.writeheader()
        # rewrite full checkpoint table; length unknown for prior → blank unless we rescan
        length_map: dict[str, int] = {}
        if scores:
            for pid, seq in iter_fasta(Path(args.fasta)):
                if pid in scores:
                    length_map[pid] = len(seq)
        for pid in sorted(scores):
            w.writerow(
                {
                    "protein_id": pid,
                    "length": length_map.get(pid, lengths.get(pid, "")),
                    "frac_disorder_metapredict": round(scores[pid], 6),
                }
            )

    summary = {
        "n": len(scores),
        "threshold": args.threshold,
        "out": str(out),
        "mean_frac": round(sum(scores.values()) / max(1, len(scores)), 4) if scores else None,
        "elapsed_sec_this_run": round(time.time() - t0, 1),
        "pending_this_run": len(pending),
    }
    Path(str(out) + ".summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    if Path(args.db).exists() and scores:
        conn = sqlite3.connect(str(args.db))
        n = join_db(conn, scores)
        conn.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            ("metapredict_proteome", json.dumps(summary)),
        )
        conn.commit()
        print(f"Joined {n} rows into {args.db}", flush=True)
        conn.close()


if __name__ == "__main__":
    main()
