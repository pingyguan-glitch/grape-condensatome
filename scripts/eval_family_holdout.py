#!/usr/bin/env python3
"""Family hold-out and Ath/rice transfer AUCs for PlantLLPS (honest transfer).

Retrains logistic on train_formal_v13.csv with leave-one-family-out splits and
cross-species partitions. Reports metrics JSON alongside in-sample gold AUCs
for contrast — do not conflate the two.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plantllps.features import FEATURE_NAMES, extract_features, features_to_vector  # noqa: E402
from plantllps.metrics import roc_auc  # noqa: E402
from plantllps.model import PlantLLPSModel, load_training_csv  # noqa: E402

DEFAULT_TRAIN = ROOT / "data" / "train_formal_v13.csv"
DEFAULT_GOLD = ROOT / "data" / "benchmarks" / "gold_standard_reconstructed.csv"
DEFAULT_OUT = ROOT / "data" / "benchmarks" / "family_holdout_metrics.json"


# Heuristic family tags from source/note (resource-quality, documented)
FAMILY_RULES: List[Tuple[str, re.Pattern]] = [
    ("MAP65", re.compile(r"MAP65|map65", re.I)),
    ("ELF3_thermosensor", re.compile(r"\bELF3\b|thermosensor", re.I)),
    ("FLOE1", re.compile(r"\bFLOE1\b|hydration", re.I)),
    ("RBGD_RGG", re.compile(r"RBGD|RGGB|RGG-repeat", re.I)),
    ("SEUSS_osmotic", re.compile(r"\bSEUSS\b", re.I)),
    ("FCA_FLL", re.compile(r"\bFCA\b|\bFLL2\b", re.I)),
    ("stress_RBP_SG", re.compile(r"RBP47|UBP1|GRP|RBG|RZ1|ECT|YTH|glycine-rich|stress-granule|polyU|DEAD-box", re.I)),
    ("P_body", re.compile(r"\bDCP|VCS|EDC4|TSN1|processing-body|decapping", re.I)),
    ("photobody", re.compile(r"PHY[AB]|COP1|SPA1|PIF4|HY5|TZP|photobody|phytochrome", re.I)),
    ("dehydrin_COR", re.compile(r"dehydrin|COR|ERD1[04]|RD29|RAB18|LTI", re.I)),
    ("WVD2_IQD_MOR1", re.compile(r"WVD2|WDL|IQD|MOR1|WAVE-DAMPENED", re.I)),
    ("MAP70_TPX2_SPR1", re.compile(r"MAP70|TPX2|SPIRAL1|SPR1", re.I)),
    ("mediator_TF", re.compile(r"\bMED\d+|NF-YB|transcription factor", re.I)),
    ("evening_complex", re.compile(r"\bLUX\b|PCL1|evening", re.I)),
    ("LLPSDB_plant", re.compile(r"LLPSDB", re.I)),
    ("animal_control", re.compile(r"animal mechanism|companion_animal|PRC1|FUS|tau|hnRNPA|DDX4|BuGZ|TPX2", re.I)),
    ("hard_folded_neg", re.compile(r"hard_folded|companion_negative|tubulin|calmodulin|GAPDH|histone|ribosomal", re.I)),
    ("hard_disordered_neg", re.compile(r"hard_neg|extensin|AGP|PRP|LRX|prothymosin|PTMA", re.I)),
]


def assign_family(source: str, note: str, kingdom: str, label: str) -> str:
    blob = f"{source} {note}"
    for name, pat in FAMILY_RULES:
        if pat.search(blob):
            return name
    if kingdom == "animal":
        return "animal_control"
    if label.startswith("neg"):
        return "other_negative"
    return "other_positive"


def species_bucket(protein_id: str, note: str, kingdom: str) -> str:
    n = (note or "").lower()
    if kingdom == "animal":
        return "animal"
    if any(k in n for k in ("rice", "oryza", "osmap")):
        return "rice"
    if any(k in n for k in ("maize", "zea", "soybean", "tomato", "apple", "peach", "citrus", "grape", "vitis", "vvmap")):
        return "other_plant"
    # Common Ath UniProt patterns in this corpus
    if protein_id.startswith(("Q9", "Q8", "Q0", "Q1", "O2", "O8", "O0", "O4", "P1", "P2", "P3", "P4", "F4", "A0A")) and "rice" not in n:
        if any(k in n for k in ("arabidopsis", "atmap", "ath", "companion gold", "companion ", "llpsdb", "round_")):
            return "arabidopsis"
        if protein_id.startswith(("Q9", "Q8", "O2", "O8", "O0", "O4", "P1", "P2", "P3", "P4", "F4")):
            return "arabidopsis"
    if "arabidopsis" in n or note.startswith("At") or " At" in note:
        return "arabidopsis"
    return "other_plant"


def build_matrix(ids, seqs, labels, notes, sources, kingdoms):
    X, y, meta = [], [], []
    for i, seq in enumerate(seqs):
        feats = extract_features(seq)
        X.append(features_to_vector(feats, FEATURE_NAMES))
        y.append(labels[i])
        fam = assign_family(sources[i], notes[i], kingdoms[i], "pos" if labels[i] else "neg")
        sp = species_bucket(ids[i], notes[i], kingdoms[i])
        meta.append({"id": ids[i], "family": fam, "species": sp, "note": notes[i], "source": sources[i]})
    return X, y, meta


def load_train(path: Path):
    ids, seqs, labels, notes = [], [], [], []
    sources, kingdoms = [], []
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            ids.append(row["protein_id"])
            seqs.append(row["sequence"].strip())
            lab = row["label"].strip().lower()
            labels.append(1 if lab.startswith("pos") else 0)
            notes.append(row.get("note", ""))
            sources.append(row.get("source", ""))
            kingdoms.append(row.get("kingdom", "plant"))
    return ids, seqs, labels, notes, sources, kingdoms


def fit_and_score(Xtr, ytr, Xte, seed=7):
    m = PlantLLPSModel()
    m.fit_logistic(Xtr, ytr, seed=seed)
    s_log, s_comp, s_hyb = [], [], []
    for row in Xte:
        feats = {n: v for n, v in zip(FEATURE_NAMES, row)}
        sl = m.logistic_score(feats)
        sc = m.companion_score(feats)
        s_log.append(sl)
        s_comp.append(sc)
        s_hyb.append(0.5 * sl + 0.5 * sc)
    return s_log, s_comp, s_hyb


def auc_block(y, s_log, s_comp, s_hyb):
    return {
        "n": len(y),
        "n_pos": int(sum(y)),
        "n_neg": int(len(y) - sum(y)),
        "auc_logistic": roc_auc(y, s_log),
        "auc_companion": roc_auc(y, s_comp),
        "auc_hybrid": roc_auc(y, s_hyb),
    }


def leave_one_family_out(X, y, meta, min_pos=2):
    """Hold out one positive family; score against a shared negative pool.

    Literature families are usually class-pure. Protocol: for each positive
    family with ≥min_pos members, train on all other rows (other positives +
    all negatives), test = held-out family positives ∪ all negatives.
    """
    families = sorted({m["family"] for m in meta})
    neg_idx = [i for i, yi in enumerate(y) if yi == 0]
    results = []
    for fam in families:
        fam_pos = [i for i, m in enumerate(meta) if m["family"] == fam and y[i] == 1]
        fam_neg = [i for i, m in enumerate(meta) if m["family"] == fam and y[i] == 0]
        if len(fam_pos) < min_pos:
            results.append(
                {
                    "held_out_family": fam,
                    "skipped": True,
                    "reason": f"need ≥{min_pos} positives in family",
                    "n_pos_family": len(fam_pos),
                    "n_neg_family": len(fam_neg),
                }
            )
            continue
        tr = [i for i in range(len(y)) if i not in set(fam_pos)]
        te = list(dict.fromkeys(fam_pos + neg_idx))
        ytr = [y[i] for i in tr]
        yte = [y[i] for i in te]
        if sum(ytr) < 2 or sum(1 - yi for yi in ytr) < 2:
            results.append(
                {
                    "held_out_family": fam,
                    "skipped": True,
                    "reason": "train too small / single-class after hold-out",
                    "n_pos_family": len(fam_pos),
                }
            )
            continue
        s_log, s_comp, s_hyb = fit_and_score(
            [X[i] for i in tr], ytr, [X[i] for i in te]
        )
        block = auc_block(yte, s_log, s_comp, s_hyb)
        block.update(
            {
                "held_out_family": fam,
                "skipped": False,
                "protocol": "holdout_pos_family_plus_shared_negatives",
                "n_pos_family": len(fam_pos),
                "test_pos_ids": [meta[i]["id"] for i in fam_pos],
            }
        )
        results.append(block)
    ok = [r for r in results if not r.get("skipped")]

    def mean_key(k):
        vals = [r[k] for r in ok if r.get(k) == r.get(k)]
        return sum(vals) / len(vals) if vals else float("nan")

    return {
        "protocol": "holdout_pos_family_plus_shared_negatives",
        "folds": results,
        "n_families_evaluated": len(ok),
        "mean_auc_logistic": mean_key("auc_logistic"),
        "mean_auc_companion": mean_key("auc_companion"),
        "mean_auc_hybrid": mean_key("auc_hybrid"),
    }


def species_transfer(X, y, meta, train_species: str, test_species: str):
    tr = [i for i, m in enumerate(meta) if m["species"] == train_species]
    te = [i for i, m in enumerate(meta) if m["species"] == test_species]
    # Optionally include animal + negs from train bucket only — keep simple: train all non-test plant+neg
    # Better defined split: train = train_species positives + all negatives not in test species; test = test_species all
    neg_tr = [i for i, m in enumerate(meta) if y[i] == 0 and m["species"] != test_species]
    pos_tr = [i for i in tr if y[i] == 1]
    tr_idx = sorted(set(pos_tr + neg_tr))
    te_idx = te
    yte = [y[i] for i in te_idx]
    ytr = [y[i] for i in tr_idx]
    if sum(yte) < 2 or (len(yte) - sum(yte)) < 1 or sum(ytr) < 2:
        return {
            "train_species": train_species,
            "test_species": test_species,
            "skipped": True,
            "reason": "insufficient class counts",
            "n_train": len(tr_idx),
            "n_test": len(te_idx),
            "n_pos_test": int(sum(yte)),
            "n_neg_test": int(len(yte) - sum(yte)),
        }
    s_log, s_comp, s_hyb = fit_and_score(
        [X[i] for i in tr_idx], ytr, [X[i] for i in te_idx]
    )
    block = auc_block(yte, s_log, s_comp, s_hyb)
    block.update(
        {
            "train_species": train_species,
            "test_species": test_species,
            "skipped": False,
            "n_train": len(tr_idx),
            "test_ids": [meta[i]["id"] for i in te_idx],
        }
    )
    return block


def in_sample_gold(model_path: Path, gold_path: Path) -> dict:
    """Score reconstructed gold with deposited production model (in-sample caveat)."""
    if not model_path.exists() or not gold_path.exists():
        return {"available": False}
    model = PlantLLPSModel.load(model_path)
    y, s_log, s_comp, s_hyb, ids = [], [], [], [], []
    with gold_path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            lab = 1 if row["label"].lower().startswith("pos") else 0
            out = model.score(row["sequence"])
            y.append(lab)
            s_log.append(out["score_logistic"])
            s_comp.append(out["score_companion"])
            s_hyb.append(out["score_plantllps"])
            ids.append(row["protein_id"])
    return {
        "available": True,
        "caveat": "partially_in_sample_gold_ids_overlap_train",
        "n": len(y),
        "auc_logistic": roc_auc(y, s_log),
        "auc_companion": roc_auc(y, s_comp),
        "auc_hybrid": roc_auc(y, s_hyb),
        "ids": ids,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train-csv", type=Path, default=DEFAULT_TRAIN)
    ap.add_argument("--gold-csv", type=Path, default=DEFAULT_GOLD)
    ap.add_argument("--model", type=Path, default=ROOT / "data" / "model_formal.json")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    ids, seqs, labels, notes, sources, kingdoms = load_train(args.train_csv)
    X, y, meta = build_matrix(ids, seqs, labels, notes, sources, kingdoms)

    fam_counts = defaultdict(lambda: {"pos": 0, "neg": 0})
    sp_counts = defaultdict(lambda: {"pos": 0, "neg": 0})
    for m, yi in zip(meta, y):
        fam_counts[m["family"]]["pos" if yi else "neg"] += 1
        sp_counts[m["species"]]["pos" if yi else "neg"] += 1

    lofo = leave_one_family_out(X, y, meta)
    transfers = [
        species_transfer(X, y, meta, "arabidopsis", "rice"),
        species_transfer(X, y, meta, "arabidopsis", "other_plant"),
        species_transfer(X, y, meta, "other_plant", "arabidopsis"),
    ]

    report = {
        "version": 1,
        "train_csv": str(args.train_csv.relative_to(ROOT)),
        "n_train": len(y),
        "n_pos": int(sum(y)),
        "n_neg": int(len(y) - sum(y)),
        "family_counts": dict(fam_counts),
        "species_counts": dict(sp_counts),
        "family_holdout": lofo,
        "species_transfer": transfers,
        "in_sample_gold_production_model": in_sample_gold(args.model, args.gold_csv),
        "honesty": (
            "Family hold-out and Ath/rice transfer AUCs are the appropriate transfer "
            "claims. Production-model AUCs on the reconstructed 13+15 gold set are "
            "partially in-sample (gold IDs overlap train) and must not be advertised "
            "as independent transfer performance. External baselines on that gold set "
            "remain LLPhyScore ≈ 0.769 and native PSAP ≈ 0.718."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2))
    readme = ROOT / "data" / "benchmarks" / "FAMILY_HOLDOUT_README.md"
    readme.write_text(
        f"""# Family hold-out / cross-species transfer metrics

Generated by `scripts/eval_family_holdout.py` → `family_holdout_metrics.json`.

## Headline (v13 train)

- Family hold-out mean hybrid AUC ≈ **{lofo['mean_auc_hybrid']:.3f}** ({lofo['n_families_evaluated']} families evaluated)
- Family hold-out mean logistic AUC ≈ **{lofo['mean_auc_logistic']:.3f}**
- Ath→rice / Ath→other_plant / other_plant→Ath: see JSON `species_transfer`

## vs in-sample gold

Production model on reconstructed 13+15 gold is **partially in-sample** — see
`in_sample_gold_production_model` in the JSON. Use hold-out / transfer numbers
when claiming generalization; use gold AUCs only with the in-sample caveat.
"""
    )
    print(
        f"Family hold-out mean hybrid AUC={lofo['mean_auc_hybrid']:.3f} "
        f"logistic={lofo['mean_auc_logistic']:.3f} "
        f"({lofo['n_families_evaluated']} families) → {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
