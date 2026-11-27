#!/usr/bin/env python3
"""Score UniProt UP000009183 with native LLPhyScore (batched + checkpointed).

Requires an unpacked standalone tree and a Python env with its deps::

    export LLPHYSCORE_HOME=/path/to/standalone_package/LLPhyScore
    export LLPHYSCORE_PYTHON=/path/to/llphy39/bin/python   # optional
    python3 scripts/score_llphyscore_proteome.py --db-ids-only --batch-size 100

Writes:
  * data/benchmarks/llphyscore_proteome.tsv  (8-feature-sum human+PDB %ile)
  * data/benchmarks/llphyscore_proteome.checkpoint.json (gitignored)
  * joins score_llphyscore_percentile into SQLite when --db is set
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sqlite3
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DEFAULT_DB = DATA / "grape_condensatome.sqlite"
DEFAULT_FASTA = DATA / "UP000009183.fasta.gz"
DEFAULT_OUT = DATA / "benchmarks" / "llphyscore_proteome.tsv"
DEFAULT_CKPT = DATA / "benchmarks" / "llphyscore_proteome.checkpoint.json"
DEFAULT_HOME = os.environ.get(
    "LLPHYSCORE_HOME",
    "/tmp/LLPhyScore/standalone_package/LLPhyScore",
)


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


def run_llphy_batch(
    batch: list[tuple[str, str]],
    llphy_home: Path,
    python_exe: str,
    model: str,
    work: Path,
) -> dict[str, float]:
    script = llphy_home / "LLPhyScore_standalone.py"
    if not script.exists():
        raise FileNotFoundError(script)
    fa = work / "batch.fasta"
    out_csv = work / "scores.csv"
    with fa.open("w") as fh:
        for pid, seq in batch:
            clean = "".join(c for c in seq if c in "ACDEFGHIKLMNPQRSTVWY")
            fh.write(f">{pid}\n{clean}\n")
    cmd = [
        python_exe,
        str(script),
        "-i",
        str(fa),
        "-m",
        model,
        "-o",
        str(out_csv),
    ]
    proc = subprocess.run(
        cmd,
        cwd=str(llphy_home),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not out_csv.exists():
        raise RuntimeError(
            f"LLPhyScore failed ({proc.returncode}):\n{proc.stderr[-2000:]}\n{proc.stdout[-2000:]}"
        )
    scores: dict[str, float] = {}
    with out_csv.open() as fh:
        for row in csv.DictReader(fh):
            tag = row.get("tag") or row.get("protein_id")
            if not tag:
                continue
            # Prefer the 8-feature sum percentile (matches gold compare script)
            val = row.get("8-feature sum")
            if val is None or val == "":
                continue
            scores[tag] = float(val)
    if len(scores) < max(1, len(batch) // 2):
        raise RuntimeError(f"parsed only {len(scores)}/{len(batch)} LLPhyScore rows")
    return scores


def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    if "score_llphyscore_percentile" not in cols:
        conn.execute("ALTER TABLE proteins ADD COLUMN score_llphyscore_percentile REAL")
    if "comparator_set" not in cols:
        conn.execute("ALTER TABLE proteins ADD COLUMN comparator_set TEXT")


def join_db(conn: sqlite3.Connection, scores: dict[str, float]) -> int:
    ensure_columns(conn)
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
    return updated


def resolve_python(explicit: str | None) -> str:
    if explicit and Path(explicit).exists():
        return explicit
    for cand in (
        os.environ.get("LLPHYSCORE_PYTHON"),
        "/tmp/mamba/envs/llphy39/bin/python",
        sys.executable,
    ):
        if cand and Path(cand).exists():
            return cand
    return sys.executable


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", default=str(DEFAULT_FASTA))
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--checkpoint", default=str(DEFAULT_CKPT))
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--llphyscore-home", default=DEFAULT_HOME)
    ap.add_argument(
        "--llphyscore-python",
        default=os.environ.get("LLPHYSCORE_PYTHON"),
        help="Interpreter with numpy/scipy/pandas/sklearn/tqdm for standalone",
    )
    ap.add_argument("--model", default="human+PDB", choices=("human+PDB", "human", "PDB"))
    ap.add_argument("--db-ids-only", action="store_true")
    ap.add_argument("--skip-db-join", action="store_true")
    args = ap.parse_args()

    llphy_home = Path(args.llphyscore_home)
    python_exe = resolve_python(args.llphyscore_python)
    if not (llphy_home / "LLPhyScore_standalone.py").exists():
        raise SystemExit(
            f"LLPhyScore standalone not found at {llphy_home}. "
            "Clone julie-forman-kay-lab/LLPhyScore and unpack standalone_package."
        )

    ckpt_path = Path(args.checkpoint)
    scores = load_checkpoint(ckpt_path)
    print(
        f"Checkpoint has {len(scores)} scored IDs; home={llphy_home}; python={python_exe}",
        flush=True,
    )

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
        with tempfile.TemporaryDirectory(prefix="llphy_proteome_") as td:
            part = run_llphy_batch(batch, llphy_home, python_exe, args.model, Path(td))
        scores.update(part)
        save_checkpoint(
            ckpt_path,
            scores,
            {
                "n": len(scores),
                "model": args.model,
                "llphyscore_home": str(llphy_home),
                "llphyscore_python": python_exe,
                "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
        )
        elapsed = time.time() - t0
        done = min(i + args.batch_size, len(pending))
        rate = done / max(elapsed, 1e-6)
        eta = (len(pending) - done) / max(rate, 1e-6)
        print(
            f"    checkpoint n={len(scores)} ({done}/{len(pending)}; "
            f"{rate:.2f} seq/s; ETA {eta/60:.1f} min)",
            flush=True,
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    length_map: dict[str, int] = {}
    if scores:
        for pid, seq in iter_fasta(Path(args.fasta)):
            if pid in scores:
                length_map[pid] = len(seq)
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
        "n": len(scores),
        "out": str(out),
        "model": args.model,
        "mean_percentile": round(sum(scores.values()) / max(1, len(scores)), 4) if scores else None,
        "elapsed_sec_this_run": round(time.time() - t0, 1),
        "pending_this_run": len(pending),
        "llphyscore_home": str(llphy_home),
        "llphyscore_python": python_exe,
        "native": True,
    }
    Path(str(out) + ".summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    if not args.skip_db_join and Path(args.db).exists() and scores:
        conn = sqlite3.connect(str(args.db))
        n = join_db(conn, scores)
        conn.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            ("llphyscore_proteome", json.dumps(summary)),
        )
        conn.commit()
        print(f"Joined {n} rows into {args.db}", flush=True)
        conn.close()


if __name__ == "__main__":
    main()
