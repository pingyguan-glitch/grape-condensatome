#!/usr/bin/env python3
"""Run native PSAP on the reconstructed gold set and report AUC vs PlantLLPS.

Requires a Python 3.8–3.9 env with pinned deps (see docs/EXTERNAL_TOOLS.md)::

    micromamba create -n psap39 -c conda-forge python=3.9 pip
    micromamba run -n psap39 pip install 'psap' 'numpy<1.24' 'scipy<1.11'
    export PSAP_PYTHON=$(micromamba run -n psap39 which python)

    python3 scripts/run_psap_gold.py
    # updates data/benchmarks/gold_standard_comparison.json psap_* fields
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plantllps.metrics import roc_auc  # noqa: E402

DATA = ROOT / "data"
BENCH = DATA / "benchmarks"


def write_uniprot_fasta(rows: list[dict], path: Path) -> None:
    with path.open("w") as fh:
        for r in rows:
            note = (r.get("note") or "PROT").replace(" ", "_")
            seq = "".join(c for c in r["sequence"] if c in "ACDEFGHIKLMNPQRSTVWY")
            # PSAP matrix parser expects sp|ACC|NAME headers
            fh.write(f">sp|{r['protein_id']}|{note}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def run_psap(fasta: Path, out_dir: Path, psap_python: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    # Prefer the psap CLI next to the interpreter
    psap_bin = Path(psap_python).resolve().parent / "psap"
    if psap_bin.exists():
        cmd = [str(psap_bin), "predict", "-f", str(fasta), "-o", str(out_dir)]
    else:
        cmd = [
            psap_python,
            "-m",
            "psap.cli",
            "predict",
            "-f",
            str(fasta),
            "-o",
            str(out_dir),
        ]
        # fallback: invoke via python -c
        if not shutil.which("psap", path=str(Path(psap_python).parent)):
            cmd = [
                psap_python,
                "-c",
                (
                    "from psap.cli import main; import sys; "
                    f"sys.argv=['psap','predict','-f',r'{fasta}','-o',r'{out_dir}']; main()"
                ),
            ]
    env = os.environ.copy()
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, check=False)
    pred = out_dir / "prediction_out.csv"
    if proc.returncode != 0 or not pred.exists():
        raise RuntimeError(
            f"PSAP failed (code={proc.returncode}):\n{proc.stderr[-3000:]}\n{proc.stdout[-1000:]}"
        )
    return pred


def parse_psap_prediction(pred: Path, rows: list[dict]) -> dict[str, float]:
    note_to_acc = {}
    for r in rows:
        note_to_acc[r.get("note", "")] = r["protein_id"]
        note_to_acc[(r.get("note") or "").replace(" ", "_")] = r["protein_id"]
    scores: dict[str, float] = {}
    with pred.open() as fh:
        for row in csv.DictReader(fh):
            name = row.get("protein_name") or ""
            acc = note_to_acc.get(name) or note_to_acc.get(name.replace("_", " "))
            if not acc:
                # try uniprot from matrix sibling if present
                continue
            scores[acc] = float(row["PSAP_score"])
    # Prefer matrix file for uniprot_id when available
    for mat in pred.parent.glob("*_psap_matrix_*.csv"):
        with mat.open() as fh:
            for row in csv.DictReader(fh):
                uid = row.get("uniprot_id")
                # matrix may not include PSAP_score; skip if so
                if uid and uid in {r["protein_id"] for r in rows} and uid not in scores:
                    pass
        break
    if len(scores) < len(rows) // 2:
        raise RuntimeError(f"parsed only {len(scores)}/{len(rows)} PSAP scores from {pred}")
    return scores


def auc(pairs: list[tuple[float, int]]) -> float | None:
    if len(pairs) < 4:
        return None
    s = [p[0] for p in pairs]
    y = [p[1] for p in pairs]
    if len(set(y)) < 2:
        return None
    return float(roc_auc(y, s))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--psap-python",
        default=os.environ.get("PSAP_PYTHON"),
        help="Python 3.9 interpreter with psap installed",
    )
    ap.add_argument(
        "--gold-csv",
        default=str(BENCH / "gold_standard_reconstructed.csv"),
    )
    ap.add_argument(
        "--comparison-json",
        default=str(BENCH / "gold_standard_comparison.json"),
    )
    ap.add_argument(
        "--scores-tsv",
        default=str(BENCH / "gold_standard_scores.tsv"),
    )
    ap.add_argument("--keep-workdir", action="store_true")
    args = ap.parse_args()

    if not args.psap_python:
        # try common micromamba path used in this project
        candidate = Path("/tmp/mamba/envs/psap39/bin/python")
        if candidate.exists():
            args.psap_python = str(candidate)
        else:
            raise SystemExit(
                "Set --psap-python or PSAP_PYTHON to a Python 3.9 env with psap installed "
                "(see docs/EXTERNAL_TOOLS.md)."
            )

    gold_path = Path(args.gold_csv)
    if not gold_path.exists():
        print("Gold CSV missing; running compare_gold_standard.py --skip-llphyscore …", flush=True)
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "compare_gold_standard.py"), "--skip-llphyscore"],
            check=True,
        )
    rows = list(csv.DictReader(gold_path.open()))

    td_ctx: tempfile.TemporaryDirectory | None = None
    if args.keep_workdir:
        work = BENCH / "psap_workdir"
        work.mkdir(parents=True, exist_ok=True)
        td_path = work
    else:
        td_ctx = tempfile.TemporaryDirectory(prefix="psap_gold_")
        td_path = Path(td_ctx.name)
    fasta = td_path / "gold.fasta"
    out_dir = td_path / "out"
    write_uniprot_fasta(rows, fasta)
    print(f"Running PSAP via {args.psap_python} …", flush=True)
    pred = run_psap(fasta, out_dir, args.psap_python)
    dest_pred = BENCH / "psap_gold_prediction.csv"
    dest_pred.write_text(pred.read_text())
    # copy matrix if present
    for mat in out_dir.glob("*_psap_matrix_*.csv"):
        (BENCH / "psap_gold_matrix.csv").write_text(mat.read_text())
        break

    psap_scores = parse_psap_prediction(pred, rows)
    pairs = [
        (psap_scores[r["protein_id"]], 1 if r["label"] == "positive" else 0)
        for r in rows
        if r["protein_id"] in psap_scores
    ]
    psap_auc = auc(pairs)

    # merge into scores TSV
    scores_path = Path(args.scores_tsv)
    if scores_path.exists():
        scored = list(csv.DictReader(scores_path.open(), delimiter="\t"))
        fields = list(scored[0].keys()) if scored else []
        if "psap_score" not in fields:
            fields.append("psap_score")
        for r in scored:
            r["psap_score"] = psap_scores.get(r["protein_id"], "")
        with scores_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in scored:
                w.writerow(r)

    # merge into comparison JSON
    report_path = Path(args.comparison_json)
    if report_path.exists():
        report = json.loads(report_path.read_text())
    else:
        report = {"auc": {}}
    report.setdefault("auc", {})["psap"] = psap_auc
    report["psap_status"] = {
        "runnable": True,
        "n_scored": len(psap_scores),
        "psap_python": args.psap_python,
        "prediction_csv": str(dest_pred),
        "auc_gold": psap_auc,
        "cited_auc_gold": 0.692,
        "cited_spearman_map208": 0.368,
        "note": (
            "Native PSAP 1.0.7 (van Mierlo et al.) on reconstructed 13+15 gold set. "
            "Requires Python ≤3.9 + numpy<1.24 + scipy<1.11. "
            "MAP-screen Spearman ρ still blocked (208-protein table not deposited)."
        ),
    }
    # PlantLLPS AUCs for side-by-side if scores TSV present
    if scores_path.exists():
        scored = list(csv.DictReader(scores_path.open(), delimiter="\t"))

        def auc_key(key: str) -> float | None:
            pr = []
            for r in scored:
                if r.get(key) in (None, ""):
                    continue
                try:
                    pr.append((float(r[key]), 1 if r["label"] == "positive" else 0))
                except ValueError:
                    continue
            return auc(pr)

        report["auc"]["plantllps_hybrid"] = auc_key("score_plantllps") or report["auc"].get(
            "plantllps_hybrid"
        )
        report["auc"]["llphyscore_percentile"] = auc_key("llphyscore_percentile") or report[
            "auc"
        ].get("llphyscore_percentile")
    report_path.write_text(json.dumps(report, indent=2))
    print(json.dumps({"psap_auc": psap_auc, "n": len(psap_scores)}, indent=2))
    print(f"Updated {report_path}")
    if td_ctx is not None:
        td_ctx.cleanup()


if __name__ == "__main__":
    main()
