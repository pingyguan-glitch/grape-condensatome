#!/usr/bin/env python3
"""Compare PlantLLPS (+ optional LLPhyScore) on the companion mechanism-matched gold set.

Reconstructs the 13 disorder-driven positives + 15 folded negatives described in the
companion MAP/LLPS manuscript (Table S5 not shipped as a numeric table in the store).
PSAP (van Mierlo et al. 2021) cannot be executed here (PyPI `psap` pins sklearn~=0.21 /
sklearn-json incompatible with Python 3.12); companion-reported PSAP AUC=0.692 is cited.
The 208-protein grapevine MAP screen with PSAP scores is not present in the Project store,
so MAP-screen Spearman ρ cannot be re-computed — companion ρ=0.368 is cited as limitation.

LLPhyScore (standalone 20220627) is run when LLPHYSCORE_HOME points at an unpacked tree.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plantllps.metrics import roc_auc  # noqa: E402
from plantllps.model import ensure_trained_model  # noqa: E402

DATA = ROOT / "data"

# Mechanism-matched gold positives (13): companion plant drivers + animal controls.
# CRY2 intentionally excluded (folded-domain photobody; eval_mechanism_special.csv).
GOLD_POS = [
    ("Q9FLP0", "AtMAP65-1"),
    ("O82804", "ELF3"),
    ("Q8VZR8", "FLOE1"),
    ("Q7DZS7", "RBGD2"),
    ("Q8W234", "SEUSS"),
    ("O04425", "FCA"),
    ("P35637", "FUS"),
    ("O43663", "PRC1"),
    ("Q9ULW0", "TPX2"),
    ("P09651", "hnRNPA1"),
    ("Q9NQI0", "DDX4"),
    ("P10636", "tau"),
    ("Q9NP74", "BuGZ/ZNF207"),
]

GOLD_NEG = [
    ("A0ACT1EZG2", "Vv tubulin alpha"),
    ("P04350", "human tubulin beta"),
    ("P0DP23", "calmodulin"),
    ("P25858", "At GAPDH"),
    ("Q0WL56", "At tubulin alpha"),
    ("Q9FLN4", "ribosomal"),
    ("Q9FNK4", "At OAT enzyme"),
    ("P29514", "tubulin beta-6"),
    ("Q96292", "actin-2"),
    ("P24100", "CDKA;1"),
    ("P25696", "enolase2/LOS2"),
    ("P04778", "LHCB1.3"),
    ("O48549", "ribosomal S6A"),
    ("P49040", "sucrose synthase SUS1"),
    ("Q96528", "catalase-1"),
]

# Companion-reported external baselines (not re-run here for PSAP)
COMPANION_BASELINES = {
    "psap_auc_gold": 0.692,
    "llphyscore_auc_gold": 0.774,
    "psap_spearman_map208": 0.368,
    "llphyscore_spearman_map208": 0.232,
    "companion_composite_auc_gold_insample": 1.00,
}


def load_seqs_from_train() -> dict[str, str]:
    seqs: dict[str, str] = {}
    for name in (
        "train_formal_v13.csv",
        "train_formal_v12.csv",
        "train_formal.csv",
        "training/train_formal_v12.csv",
    ):
        p = DATA / name
        if not p.exists():
            continue
        for row in csv.DictReader(p.open()):
            seqs[row["protein_id"]] = row["sequence"].strip().upper()
    return seqs


def write_gold_csv(path: Path, seqs: dict[str, str]) -> list[dict]:
    rows = []
    missing = []
    for acc, note in GOLD_POS:
        if acc not in seqs:
            missing.append(acc)
            continue
        rows.append(
            {
                "protein_id": acc,
                "label": "positive",
                "note": note,
                "sequence": seqs[acc],
                "in_plantllps_train": None,
            }
        )
    for acc, note in GOLD_NEG:
        if acc not in seqs:
            missing.append(acc)
            continue
        rows.append(
            {
                "protein_id": acc,
                "label": "negative",
                "note": note,
                "sequence": seqs[acc],
                "in_plantllps_train": None,
            }
        )
    if missing:
        raise SystemExit(f"missing sequences for gold IDs: {missing}")

    train_ids = set()
    for name in ("train_formal_v13.csv", "train_formal_v12.csv"):
        p = DATA / name
        if p.exists():
            train_ids |= {r["protein_id"] for r in csv.DictReader(p.open())}
            break
    for r in rows:
        r["in_plantllps_train"] = r["protein_id"] in train_ids

    fields = ["protein_id", "label", "note", "in_plantllps_train", "sequence"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({f: r[f] for f in fields})
    return rows


def score_plantllps(rows: list[dict], model_path: Path) -> list[dict]:
    model = ensure_trained_model(ROOT)
    if model_path.exists():
        from plantllps.model import PlantLLPSModel

        model = PlantLLPSModel.load(model_path)
    out = []
    for r in rows:
        s = model.score(r["sequence"])
        out.append(
            {
                **{k: r[k] for k in ("protein_id", "label", "note", "in_plantllps_train")},
                "score_plantllps": float(s["score_plantllps"]),
                "score_companion": float(s["score_companion"]),
                "score_logistic": float(s["score_logistic"]),
                "map65_tail_like": float(s.get("map65_tail_like", 0)),
                "frac_disorder_proxy": float(s.get("frac_disorder_proxy", 0)),
            }
        )
    return out


def run_llphyscore(rows: list[dict], llphy_home: Path) -> dict[str, float]:
    """Return protein_id -> 8-feature-sum human-proteome percentile."""
    script = llphy_home / "LLPhyScore_standalone.py"
    if not script.exists():
        raise FileNotFoundError(script)
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        fa = td_path / "gold.fasta"
        out_csv = td_path / "llphy_scores.csv"
        with fa.open("w") as fh:
            for r in rows:
                # LLPhyScore fails on non-canonical residues
                seq = "".join(c for c in r["sequence"] if c in "ACDEFGHIKLMNPQRSTVWY")
                fh.write(f">{r['protein_id']}\n{seq}\n")
        cmd = [
            sys.executable,
            str(script),
            "-i",
            str(fa),
            "-m",
            "human+PDB",
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
        if proc.returncode != 0:
            raise RuntimeError(
                f"LLPhyScore failed ({proc.returncode}):\n{proc.stderr[-2000:]}\n{proc.stdout[-2000:]}"
            )
        if not out_csv.exists():
            raise RuntimeError(f"LLPhyScore produced no CSV at {out_csv}")
        # Persist a copy for the repo
        dest = DATA / "benchmarks" / "llphyscore_gold_raw.csv"
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(out_csv.read_text())
        scores: dict[str, float] = {}
        with out_csv.open() as fh:
            for row in csv.DictReader(fh):
                tag = row.get("tag") or row.get("protein_id")
                if not tag:
                    continue
                scores[tag] = float(row["8-feature sum"])
        if len(scores) < len(rows) // 2:
            raise RuntimeError(f"parsed only {len(scores)}/{len(rows)} LLPhyScore rows")
        return scores


def auc_for(scored: list[dict], key: str) -> float | None:
    pairs = [
        (r[key], 1 if r["label"] == "positive" else 0)
        for r in scored
        if r.get(key) is not None and not (
            isinstance(r.get(key), float) and math.isnan(float(r[key]))
        )
    ]
    if len(pairs) < 4:
        return None
    s = [p[0] for p in pairs]
    y = [p[1] for p in pairs]
    if len(set(y)) < 2:
        return None
    return float(roc_auc(y, s))  # type: ignore[arg-type]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--llphyscore-home",
        default=os.environ.get(
            "LLPHYSCORE_HOME",
            "/tmp/LLPhyScore/standalone_package/LLPhyScore",
        ),
    )
    ap.add_argument("--skip-llphyscore", action="store_true")
    ap.add_argument(
        "--model",
        default=str(DATA / "model_formal.json"),
    )
    args = ap.parse_args()

    out_dir = DATA / "benchmarks"
    out_dir.mkdir(parents=True, exist_ok=True)

    seqs = load_seqs_from_train()
    gold_csv = out_dir / "gold_standard_reconstructed.csv"
    rows = write_gold_csv(gold_csv, seqs)

    scored = score_plantllps(rows, Path(args.model))
    llphy_scores: dict[str, float] = {}
    llphy_error = None
    if not args.skip_llphyscore:
        home = Path(args.llphyscore_home)
        try:
            llphy_scores = run_llphyscore(rows, home)
            for r in scored:
                r["llphyscore_percentile"] = llphy_scores.get(r["protein_id"])
        except Exception as e:  # noqa: BLE001 — report limitation honestly
            llphy_error = str(e)
            for r in scored:
                r["llphyscore_percentile"] = None

    for r in scored:
        r["psap_score"] = None  # not runnable in this environment

    # write per-protein table
    tsv = out_dir / "gold_standard_scores.tsv"
    fields = [
        "protein_id",
        "label",
        "note",
        "in_plantllps_train",
        "score_plantllps",
        "score_companion",
        "score_logistic",
        "llphyscore_percentile",
        "psap_score",
        "map65_tail_like",
        "frac_disorder_proxy",
    ]
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore", delimiter="\t")
        w.writeheader()
        for r in scored:
            w.writerow(r)

    n_overlap = sum(1 for r in scored if r["in_plantllps_train"])
    report = {
        "n": len(scored),
        "n_pos": sum(1 for r in scored if r["label"] == "positive"),
        "n_neg": sum(1 for r in scored if r["label"] == "negative"),
        "n_in_plantllps_train": n_overlap,
        "plantllps_train_overlap_warning": (
            "Gold-standard IDs overlap the PlantLLPS training corpus; "
            "PlantLLPS AUCs below are partially in-sample and must not be claimed as "
            "independent transfer benchmarks. LLPhyScore/PSAP are out-of-the-box."
        ),
        "auc": {
            "plantllps_hybrid": auc_for(scored, "score_plantllps"),
            "plantllps_logistic": auc_for(scored, "score_logistic"),
            "companion_composite": auc_for(scored, "score_companion"),
            "llphyscore_percentile": auc_for(scored, "llphyscore_percentile")
            if llphy_scores
            else None,
            "psap": None,
        },
        "companion_reported_baselines": COMPANION_BASELINES,
        "psap_status": {
            "runnable": False,
            "reason": (
                "PyPI package `psap` 1.0.7 requires sklearn~=0.21 and sklearn-json; "
                "sklearn-json fails to import on scikit-learn bundled with Python 3.12 "
                "(cannot import `_gb_losses`). Wheel extracts but model JSON cannot load."
            ),
            "cited_auc_gold": COMPANION_BASELINES["psap_auc_gold"],
            "cited_spearman_map208": COMPANION_BASELINES["psap_spearman_map208"],
        },
        "llphyscore_status": {
            "runnable": bool(llphy_scores) and llphy_error is None,
            "n_scored": len(llphy_scores),
            "error": llphy_error,
            "home": str(args.llphyscore_home),
            "cited_auc_gold": COMPANION_BASELINES["llphyscore_auc_gold"],
        },
        "map_screen_status": {
            "runnable": False,
            "reason": (
                "Companion grapevine MAP screen table (208 proteins with PSAP / "
                "LLPhyScore / composite columns) is not present in the Project store "
                "or GitHub repo; only narrative metrics in the manuscript docx. "
                "Cannot recompute Spearman ρ vs PSAP without that table."
            ),
            "cited_psap_rho": COMPANION_BASELINES["psap_spearman_map208"],
            "cited_llphyscore_rho": COMPANION_BASELINES["llphyscore_spearman_map208"],
        },
        "anchor_ranks": {},
        "files": {"gold_csv": str(gold_csv), "scores_tsv": str(tsv)},
    }

    # ranks for key anchors among scored set
    by_score = sorted(scored, key=lambda r: -r["score_plantllps"])
    for acc in ("Q9FLP0", "F6HET4", "O43663", "P35637"):
        for i, r in enumerate(by_score, 1):
            if r["protein_id"] == acc:
                report["anchor_ranks"][acc] = {
                    "rank_plantllps": i,
                    "n": len(by_score),
                    "score": r["score_plantllps"],
                    "map65_tail_like": r["map65_tail_like"],
                    "llphyscore_percentile": r.get("llphyscore_percentile"),
                }
                break

    # F6HET4 may not be in gold set — still report if in train seqs via model score
    if "F6HET4" not in report["anchor_ranks"] and "F6HET4" in seqs:
        model = ensure_trained_model(ROOT)
        from plantllps.model import PlantLLPSModel

        if Path(args.model).exists():
            model = PlantLLPSModel.load(Path(args.model))
        s = model.score(seqs["F6HET4"])
        report["anchor_ranks"]["F6HET4_not_in_gold"] = {
            "score_plantllps": float(s["score_plantllps"]),
            "map65_tail_like": float(s.get("map65_tail_like", 0)),
            "note": "VvMAP65-1 scored with PlantLLPS but not part of 13+15 gold set",
        }

    out_json = out_dir / "gold_standard_comparison.json"
    out_json.write_text(json.dumps(report, indent=2))
    print(json.dumps(report["auc"], indent=2))
    print(f"overlap_with_train={n_overlap}/{len(scored)}")
    print(f"Wrote {tsv}\nWrote {out_json}")
    if llphy_error:
        print("LLPhyScore error:", llphy_error)


if __name__ == "__main__":
    main()
