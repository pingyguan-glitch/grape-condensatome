#!/usr/bin/env python3
"""Score a protein subset with metapredict and write a joinable TSV.

Default subset = reconstructed gold standards (+ optional UniProt IDs / Tier A sample).
Requires metapredict in the current interpreter, or METAPREDICT_PYTHON /
--metapredict-python pointing at an env that has it.

Example::

    export METAPREDICT_PYTHON=/path/to/meta311/bin/python
    python3 scripts/score_metapredict_subset.py
    # → data/benchmarks/metapredict_subset.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plantllps.metapredict_bridge import resolve_metapredict_scores  # noqa: E402

DATA = ROOT / "data"


def load_gold_seqs() -> list[tuple[str, str, str]]:
    path = DATA / "benchmarks" / "gold_standard_reconstructed.csv"
    if not path.exists():
        # fall back to train rows for known gold IDs
        from scripts.compare_gold_standard import GOLD_NEG, GOLD_POS, load_seqs_from_train

        seqs = load_seqs_from_train()
        rows = []
        for acc, note in GOLD_POS + GOLD_NEG:
            if acc in seqs:
                rows.append((acc, seqs[acc], note))
        return rows
    out = []
    for row in csv.DictReader(path.open()):
        out.append((row["protein_id"], row["sequence"].strip(), row.get("note", "")))
    return out


def load_extra_ids(ids: list[str]) -> list[tuple[str, str, str]]:
    seqs: dict[str, str] = {}
    for name in ("train_formal_v13.csv", "train_formal.csv", "train_formal_v12.csv"):
        p = DATA / name
        if not p.exists():
            continue
        for row in csv.DictReader(p.open()):
            seqs[row["protein_id"]] = row["sequence"].strip()
    # also try proteome FASTA if present
    fa = DATA / "UP000009183.fasta"
    if not fa.exists() and (DATA / "UP000009183.fasta.gz").exists():
        import gzip

        with gzip.open(DATA / "UP000009183.fasta.gz", "rt") as fh:
            pid = None
            buf: list[str] = []
            for line in fh:
                if line.startswith(">"):
                    if pid and pid in ids:
                        seqs[pid] = "".join(buf)
                    header = line[1:].split()[0]
                    pid = header.split("|")[1] if header.startswith("sp|") or "|" in header else header
                    if "|" in header:
                        parts = header.split("|")
                        pid = parts[1] if len(parts) > 1 else parts[0]
                    buf = []
                else:
                    buf.append(line.strip())
            if pid and pid in ids:
                seqs[pid] = "".join(buf)
    rows = []
    for i in ids:
        if i in seqs:
            rows.append((i, seqs[i], "extra"))
        else:
            print(f"warning: no sequence for {i}", file=sys.stderr)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out",
        default=str(DATA / "benchmarks" / "metapredict_subset.tsv"),
    )
    ap.add_argument("--extra-id", action="append", default=[], help="Extra UniProt accession")
    ap.add_argument(
        "--metapredict-python",
        default=os.environ.get("METAPREDICT_PYTHON"),
    )
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument(
        "--summary-json",
        default=str(DATA / "benchmarks" / "metapredict_subset.summary.json"),
    )
    args = ap.parse_args()

    rows = load_gold_seqs()
    if args.extra_id:
        rows.extend(load_extra_ids(args.extra_id))
    # always try VvMAP65-1
    if not any(r[0] == "F6HET4" for r in rows):
        rows.extend(load_extra_ids(["F6HET4"]))

    # dedupe
    seen = set()
    uniq = []
    for r in rows:
        if r[0] in seen:
            continue
        seen.add(r[0])
        uniq.append(r)
    rows = uniq
    print(f"Scoring {len(rows)} sequences with metapredict …", flush=True)
    scores = resolve_metapredict_scores(
        [(r[0], r[1]) for r in rows],
        threshold=args.threshold,
        python_exe=args.metapredict_python,
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["protein_id", "note", "length", "frac_disorder_metapredict"],
            delimiter="\t",
        )
        w.writeheader()
        for pid, seq, note in rows:
            w.writerow(
                {
                    "protein_id": pid,
                    "note": note,
                    "length": len(seq),
                    "frac_disorder_metapredict": round(scores.get(pid, 0.0), 6),
                }
            )
    summary = {
        "n": len(rows),
        "threshold": args.threshold,
        "metapredict_python": args.metapredict_python or sys.executable,
        "out": str(out),
        "mean_frac": round(sum(scores.values()) / max(1, len(scores)), 4),
        "ids": sorted(scores.keys()),
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2))
    print(f"Wrote {out}")
    print(f"Wrote {args.summary_json}")
    print(json.dumps({"n": summary["n"], "mean_frac": summary["mean_frac"]}, indent=2))


if __name__ == "__main__":
    main()
