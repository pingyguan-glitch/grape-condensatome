#!/usr/bin/env python3
"""Score UniProt UP000009183 with native PSAP (batched + checkpointed).

Requires a Python ≤3.9 env with pinned deps (see docs/EXTERNAL_TOOLS.md)::

    export PSAP_PYTHON=/path/to/psap39/bin/python
    python3 scripts/score_psap_proteome.py --db-ids-only --batch-size 500

Writes:
  * data/benchmarks/psap_proteome.tsv
  * data/benchmarks/psap_proteome.checkpoint.json (gitignored)
  * joins score_psap into SQLite when --db is set
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
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
DEFAULT_OUT = DATA / "benchmarks" / "psap_proteome.tsv"
DEFAULT_CKPT = DATA / "benchmarks" / "psap_proteome.checkpoint.json"


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


def write_uniprot_fasta(batch: list[tuple[str, str]], path: Path) -> None:
    with path.open("w") as fh:
        for pid, seq in batch:
            clean = "".join(c for c in seq if c in "ACDEFGHIKLMNPQRSTVWY")
            # NAME == accession so prediction_out protein_name maps 1:1
            fh.write(f">sp|{pid}|{pid}\n")
            for i in range(0, len(clean), 60):
                fh.write(clean[i : i + 60] + "\n")


def run_psap_batch(batch: list[tuple[str, str]], psap_python: str, work: Path) -> dict[str, float]:
    fasta = work / "batch.fasta"
    out_dir = work / "out"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    write_uniprot_fasta(batch, fasta)

    psap_bin = Path(psap_python).resolve().parent / "psap"
    if psap_bin.exists():
        cmd = [str(psap_bin), "predict", "-f", str(fasta), "-o", str(out_dir)]
    else:
        cmd = [
            psap_python,
            "-c",
            (
                "from psap.cli import main; import sys; "
                f"sys.argv=['psap','predict','-f',r'{fasta}','-o',r'{out_dir}']; main()"
            ),
        ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    preds = list(out_dir.glob("prediction_*.csv"))
    if proc.returncode != 0 or not preds:
        raise RuntimeError(
            f"PSAP failed (code={proc.returncode}):\n{proc.stderr[-3000:]}\n{proc.stdout[-1000:]}"
        )
    pred = preds[0]
    # Prefer matrix for uniprot_id; prediction has protein_name + PSAP_score
    name_to_score: dict[str, float] = {}
    with pred.open() as fh:
        for row in csv.DictReader(fh):
            name = (row.get("protein_name") or "").strip()
            if not name:
                continue
            name_to_score[name] = float(row["PSAP_score"])

    scores: dict[str, float] = {}
    mats = list(out_dir.glob("*_psap_matrix_*.csv"))
    if mats:
        with mats[0].open() as fh:
            for row in csv.DictReader(fh):
                uid = (row.get("uniprot_id") or "").strip()
                pname = (row.get("protein_name") or "").strip()
                if uid and pname in name_to_score:
                    scores[uid] = name_to_score[pname]
                elif uid and uid in name_to_score:
                    scores[uid] = name_to_score[uid]
    # Fallback: protein_name == accession
    for name, sc in name_to_score.items():
        if name not in scores:
            scores[name] = sc
    if len(scores) < max(1, len(batch) // 2):
        raise RuntimeError(f"parsed only {len(scores)}/{len(batch)} PSAP scores from {pred}")
    return scores


def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    if "score_psap" not in cols:
        conn.execute("ALTER TABLE proteins ADD COLUMN score_psap REAL")
    if "comparator_set" not in cols:
        conn.execute("ALTER TABLE proteins ADD COLUMN comparator_set TEXT")


def join_db(conn: sqlite3.Connection, scores: dict[str, float]) -> int:
    ensure_columns(conn)
    updated = 0
    for pid, sc in scores.items():
        cur = conn.execute(
            """
            UPDATE proteins SET
                score_psap = ?,
                comparator_set = COALESCE(comparator_set, 'proteome_native')
            WHERE protein_id = ?
            """,
            (float(sc), pid),
        )
        updated += cur.rowcount
    return updated


def resolve_psap_python(explicit: str | None) -> str:
    if explicit:
        return explicit
    for cand in (
        os.environ.get("PSAP_PYTHON"),
        "/tmp/mamba/envs/psap39/bin/python",
    ):
        if cand and Path(cand).exists():
            return cand
    raise SystemExit(
        "Set --psap-python or PSAP_PYTHON to a Python ≤3.9 env with psap installed "
        "(see docs/EXTERNAL_TOOLS.md)."
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", default=str(DEFAULT_FASTA))
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--checkpoint", default=str(DEFAULT_CKPT))
    ap.add_argument("--batch-size", type=int, default=500)
    ap.add_argument("--limit", type=int, default=0, help="Max new sequences this run")
    ap.add_argument("--psap-python", default=os.environ.get("PSAP_PYTHON"))
    ap.add_argument("--db-ids-only", action="store_true")
    ap.add_argument("--skip-db-join", action="store_true")
    ap.add_argument("--workdir", default="", help="Keep batch workdirs under this path")
    args = ap.parse_args()

    psap_python = resolve_psap_python(args.psap_python)
    ckpt_path = Path(args.checkpoint)
    scores = load_checkpoint(ckpt_path)
    print(f"Checkpoint has {len(scores)} scored IDs; PSAP_PYTHON={psap_python}", flush=True)

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
    keep_root = Path(args.workdir) if args.workdir else None
    if keep_root:
        keep_root.mkdir(parents=True, exist_ok=True)

    for i in range(0, len(pending), args.batch_size):
        batch = pending[i : i + args.batch_size]
        print(f"  batch {i // args.batch_size + 1}: {len(batch)} seqs …", flush=True)
        if keep_root:
            work = keep_root / f"batch_{i // args.batch_size + 1:05d}"
            work.mkdir(parents=True, exist_ok=True)
            part = run_psap_batch(batch, psap_python, work)
        else:
            with tempfile.TemporaryDirectory(prefix="psap_proteome_") as td:
                part = run_psap_batch(batch, psap_python, Path(td))
        scores.update(part)
        save_checkpoint(
            ckpt_path,
            scores,
            {
                "n": len(scores),
                "psap_python": psap_python,
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
    length_map: dict[str, int] = {}
    if scores:
        for pid, seq in iter_fasta(Path(args.fasta)):
            if pid in scores:
                length_map[pid] = len(seq)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["protein_id", "length", "score_psap"], delimiter="\t")
        w.writeheader()
        for pid in sorted(scores):
            w.writerow(
                {
                    "protein_id": pid,
                    "length": length_map.get(pid, ""),
                    "score_psap": round(scores[pid], 8),
                }
            )

    summary = {
        "n": len(scores),
        "out": str(out),
        "mean_score": round(sum(scores.values()) / max(1, len(scores)), 6) if scores else None,
        "elapsed_sec_this_run": round(time.time() - t0, 1),
        "pending_this_run": len(pending),
        "psap_python": psap_python,
        "native": True,
    }
    Path(str(out) + ".summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    if not args.skip_db_join and Path(args.db).exists() and scores:
        conn = sqlite3.connect(str(args.db))
        n = join_db(conn, scores)
        conn.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            ("psap_proteome", json.dumps(summary)),
        )
        conn.commit()
        print(f"Joined {n} rows into {args.db}", flush=True)
        conn.close()


if __name__ == "__main__":
    main()
