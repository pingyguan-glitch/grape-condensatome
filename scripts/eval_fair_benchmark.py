#!/usr/bin/env python3
"""NAR+ fair comparison suite: PlantLLPS vs PSAP / LLPhyScore / PLAAC / ProX-style.

Deposits leakage-controlled AUCs on the same panels:
  * reconstructed 13+15 gold (leave-gold-out retrain for supervised methods)
  * family hold-out (shared negatives)
  * expanded cross-species transfers (Ath↔rice / maize / soybean / crop pool)

Strict leakage control
----------------------
* Production ``model_formal.json`` is **never** used for fair PlantLLPS numbers.
* Supervised methods (PlantLLPS logistic/hybrid, ProX-style logistic) retrain
  on train indices that exclude every test ID for that fold/panel.
* PSAP, LLPhyScore, and PLAAC-approx are out-of-the-box (no plant-label fit).
* In-sample production gold AUCs are reported only under an explicit caveat key.

ProX-style baseline = stdlib logistic on ProX/PU-like features
(IDR/LCR/hydropathy/FCR/PLAAC/aromatic/κ) — not the published LightGBM/PU model.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from plantllps.features import FEATURE_NAMES, extract_features, features_to_vector  # noqa: E402
from plantllps.metrics import roc_auc  # noqa: E402
from plantllps.model import PlantLLPSModel  # noqa: E402
from plantllps.prox_features import (  # noqa: E402
    PROX_FEATURE_NAMES,
    extract_prox_features,
    prox_features_to_vector,
)

# Reuse family tagging from eval_family_holdout (keep in sync)
from eval_family_holdout import assign_family  # noqa: E402

DEFAULT_TRAIN = ROOT / "data" / "train_formal_v13.csv"
DEFAULT_GOLD = ROOT / "data" / "benchmarks" / "gold_standard_reconstructed.csv"
DEFAULT_OUT = ROOT / "data" / "benchmarks" / "fair_comparison_metrics.json"
DEFAULT_SCORES = ROOT / "data" / "benchmarks" / "fair_comparison_scores.tsv"
DEFAULT_CACHE = ROOT / "data" / "benchmarks" / "train_v13_oob_comparator_scores.tsv"

METHODS = [
    "plantllps_hybrid",
    "plantllps_logistic",
    "prox_pu_style",
    "plaac_approx",
    "psap",
    "llphyscore",
]


# --- species tagging (explicit notes; avoid glycine-rich → soybean false hits) ---

_SPECIES_PATTERNS: List[Tuple[str, re.Pattern]] = [
    ("animal", re.compile(r"\banimal\b|companion_animal|\bFUS\b|\btau\b|hnRNPA|DDX4|BuGZ|PRC1", re.I)),
    ("rice", re.compile(r"\brice\b|\boryza\b|\bosmap", re.I)),
    ("maize", re.compile(r"\bmaize\b|\bzea\b", re.I)),
    ("soybean", re.compile(r"\bsoybean\b|\bglycine max\b|\bsoy\b", re.I)),
    ("tomato", re.compile(r"\btomato\b|\bsolanum\b", re.I)),
    ("grape", re.compile(r"\bgrape\b|\bvitis\b|\bvv(?:map|tubulin)?|\bVv[A-Z0-9]", re.I)),
    ("apple_peach_citrus", re.compile(r"\bapple\b|\bpeach\b|\bcitrus\b|\bpear\b", re.I)),
    ("arabidopsis", re.compile(r"arabidopsis|\bath\b|\batmap|companion gold|llpsdb|round_", re.I)),
]


def species_bucket(protein_id: str, note: str, kingdom: str) -> str:
    blob = note or ""
    if kingdom == "animal":
        return "animal"
    for name, pat in _SPECIES_PATTERNS:
        if name == "arabidopsis":
            continue
        if pat.search(blob):
            return name
    if _SPECIES_PATTERNS[-1][1].search(blob):  # arabidopsis
        return "arabidopsis"
    # UniProt Ath-ish accessions common in this corpus when note is thin
    if protein_id.startswith(("Q9", "Q8", "Q0", "Q1", "O2", "O8", "O0", "O4", "P1", "P2", "P3", "P4", "F4")):
        return "arabidopsis"
    if protein_id.startswith("A0A") and "rice" not in blob.lower():
        return "other_plant"
    return "other_plant"


def crop_pool(species: str) -> str:
    """Coarser bucket for expanded transfer panels."""
    if species == "arabidopsis":
        return "arabidopsis"
    if species == "animal":
        return "animal"
    if species in ("rice", "maize"):
        return "cereal"
    if species in ("soybean", "tomato", "grape", "apple_peach_citrus", "other_plant"):
        return "non_ath_crop"
    return species


# --- IO ---

def load_train(path: Path):
    ids, seqs, labels, notes, sources, kingdoms = [], [], [], [], [], []
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


def load_gold(path: Path) -> List[dict]:
    rows = []
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            rows.append(row)
    return rows


def build_matrices(ids, seqs, labels, notes, sources, kingdoms):
    X_pl, X_px, y, meta = [], [], [], []
    plaac_scores = []
    for i, seq in enumerate(seqs):
        feats = extract_features(seq)
        X_pl.append(features_to_vector(feats, FEATURE_NAMES))
        px = extract_prox_features(seq)
        X_px.append(prox_features_to_vector(px, PROX_FEATURE_NAMES))
        plaac_scores.append(float(px["plaac_nllr_approx"]))
        y.append(labels[i])
        fam = assign_family(sources[i], notes[i], kingdoms[i], "pos" if labels[i] else "neg")
        sp = species_bucket(ids[i], notes[i], kingdoms[i])
        meta.append(
            {
                "id": ids[i],
                "family": fam,
                "species": sp,
                "crop_pool": crop_pool(sp),
                "note": notes[i],
                "source": sources[i],
            }
        )
    return X_pl, X_px, y, meta, plaac_scores


# --- supervised scoring ---

def fit_plantllps(Xtr, ytr, seed=7) -> PlantLLPSModel:
    m = PlantLLPSModel()
    m.fit_logistic(Xtr, ytr, seed=seed)
    return m


def fit_prox(Xtr, ytr, seed=7) -> PlantLLPSModel:
    """Reuse PlantLLPSModel SGD logistic on ProX feature matrix."""
    m = PlantLLPSModel()
    m.feature_names = list(PROX_FEATURE_NAMES)
    m.fit_logistic(Xtr, ytr, seed=seed)
    return m


def score_plantllps_rows(model: PlantLLPSModel, Xrows: Sequence[Sequence[float]]):
    s_log, s_hyb = [], []
    for row in Xrows:
        feats = {n: v for n, v in zip(FEATURE_NAMES, row)}
        sl = model.logistic_score(feats)
        sc = model.companion_score(feats)
        s_log.append(sl)
        s_hyb.append(0.5 * sl + 0.5 * sc)
    return s_log, s_hyb


def score_prox_rows(model: PlantLLPSModel, Xrows: Sequence[Sequence[float]]):
    out = []
    for row in Xrows:
        feats = {n: v for n, v in zip(PROX_FEATURE_NAMES, row)}
        out.append(model.logistic_score(feats))
    return out


def safe_auc(y: Sequence[int], scores: Sequence[Optional[float]]) -> Optional[float]:
    pairs = [(s, yi) for s, yi in zip(scores, y) if s is not None and not (isinstance(s, float) and math.isnan(s))]
    if len(pairs) < 4:
        return None
    ys = [p[1] for p in pairs]
    if len(set(ys)) < 2:
        return None
    return float(roc_auc([p[1] for p in pairs], [p[0] for p in pairs]))


def mean_finite(vals: Sequence[Optional[float]]) -> Optional[float]:
    xs = [v for v in vals if v is not None and v == v]
    return sum(xs) / len(xs) if xs else None


# --- OOB comparator runners ---

def run_psap_on_rows(rows: List[dict], psap_python: str) -> Dict[str, float]:
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        fa = td_path / "query.fasta"
        out_dir = td_path / "out"
        out_dir.mkdir()
        with fa.open("w") as fh:
            for r in rows:
                note = (r.get("note") or r["protein_id"]).replace(" ", "_")
                seq = "".join(c for c in r["sequence"] if c in "ACDEFGHIKLMNPQRSTVWY")
                fh.write(f">sp|{r['protein_id']}|{note}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i : i + 60] + "\n")
        psap_bin = Path(psap_python).resolve().parent / "psap"
        if psap_bin.exists():
            cmd = [str(psap_bin), "predict", "-f", str(fa), "-o", str(out_dir)]
        else:
            cmd = [
                psap_python,
                "-c",
                (
                    "from psap.cli import main; import sys; "
                    f"sys.argv=['psap','predict','-f',r'{fa}','-o',r'{out_dir}']; main()"
                ),
            ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        pred = out_dir / "prediction_out.csv"
        if proc.returncode != 0 or not pred.exists():
            raise RuntimeError(
                f"PSAP failed ({proc.returncode}): {proc.stderr[-2000:]}\n{proc.stdout[-1000:]}"
            )
        # Map by protein_name note OR uniprot from matrix
        note_to_acc = {}
        for r in rows:
            note_to_acc[(r.get("note") or r["protein_id"]).replace(" ", "_")] = r["protein_id"]
            note_to_acc[r.get("note") or ""] = r["protein_id"]
            note_to_acc[r["protein_id"]] = r["protein_id"]
        scores: Dict[str, float] = {}
        with pred.open() as fh:
            for row in csv.DictReader(fh):
                name = row.get("protein_name") or ""
                acc = note_to_acc.get(name) or note_to_acc.get(name.replace("_", " "))
                if acc:
                    scores[acc] = float(row["PSAP_score"])
        # Prefer matrix uniprot_id when present
        for mat in out_dir.glob("*_psap_matrix_*.csv"):
            with mat.open() as fh:
                for row in csv.DictReader(fh):
                    uid = row.get("uniprot_id")
                    if uid and uid in {r["protein_id"] for r in rows}:
                        # matrix lacks PSAP_score; keep prediction_out mapping
                        pass
            break
        if len(scores) < max(1, len(rows) // 2):
            raise RuntimeError(f"PSAP parsed {len(scores)}/{len(rows)}")
        return scores


def run_llphyscore_on_rows(
    rows: List[dict], llphy_home: Path, llphy_python: Optional[str] = None
) -> Dict[str, float]:
    script = llphy_home / "LLPhyScore_standalone.py"
    if not script.exists():
        raise FileNotFoundError(script)
    py = llphy_python or sys.executable
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        fa = td_path / "query.fasta"
        out_csv = td_path / "llphy_scores.csv"
        with fa.open("w") as fh:
            for r in rows:
                seq = "".join(c for c in r["sequence"] if c in "ACDEFGHIKLMNPQRSTVWY")
                fh.write(f">{r['protein_id']}\n{seq}\n")
        cmd = [py, str(script), "-i", str(fa), "-m", "human+PDB", "-o", str(out_csv)]
        proc = subprocess.run(
            cmd, cwd=str(llphy_home), capture_output=True, text=True, check=False
        )
        if proc.returncode != 0 or not out_csv.exists():
            raise RuntimeError(
                f"LLPhyScore failed ({proc.returncode}): {proc.stderr[-2000:]}\n{proc.stdout[-1000:]}"
            )
        scores: Dict[str, float] = {}
        with out_csv.open() as fh:
            for row in csv.DictReader(fh):
                tag = row.get("tag") or row.get("protein_id")
                if tag:
                    scores[tag] = float(row["8-feature sum"])
        if len(scores) < max(1, len(rows) // 2):
            raise RuntimeError(f"LLPhyScore parsed {len(scores)}/{len(rows)}")
        return scores


def load_deposited_gold_oob() -> Dict[str, Dict[str, float]]:
    """protein_id -> {psap, llphyscore} from deposited gold TSV."""
    path = ROOT / "data" / "benchmarks" / "gold_standard_scores.tsv"
    out: Dict[str, Dict[str, float]] = {}
    if not path.exists():
        return out
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            d: Dict[str, float] = {}
            if row.get("psap_score"):
                d["psap"] = float(row["psap_score"])
            if row.get("llphyscore_percentile"):
                d["llphyscore"] = float(row["llphyscore_percentile"])
            out[row["protein_id"]] = d
    return out


def ensure_train_oob_cache(
    ids: List[str],
    seqs: List[str],
    notes: List[str],
    cache_path: Path,
    psap_python: Optional[str],
    llphy_home: Optional[Path],
    llphy_python: Optional[str],
    force: bool = False,
) -> Dict[str, Dict[str, Optional[float]]]:
    """Cache PSAP / LLPhyScore for the full v13 train set (for hold-out panels)."""
    cached: Dict[str, Dict[str, Optional[float]]] = {}
    if cache_path.exists() and not force:
        with cache_path.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                cached[row["protein_id"]] = {
                    "psap": float(row["psap"]) if row.get("psap") not in ("", None) else None,
                    "llphyscore": float(row["llphyscore"])
                    if row.get("llphyscore") not in ("", None)
                    else None,
                }
        if len(cached) >= len(ids) * 0.9:
            return cached

    rows = [
        {"protein_id": i, "sequence": s, "note": n}
        for i, s, n in zip(ids, seqs, notes)
    ]
    psap_map: Dict[str, float] = {}
    llphy_map: Dict[str, float] = {}
    status = {"psap": None, "llphyscore": None}

    # Merge deposited gold first
    deposited = load_deposited_gold_oob()
    for pid, d in deposited.items():
        cached.setdefault(pid, {"psap": None, "llphyscore": None})
        if "psap" in d:
            cached[pid]["psap"] = d["psap"]
        if "llphyscore" in d:
            cached[pid]["llphyscore"] = d["llphyscore"]

    missing_psap = [r for r in rows if cached.get(r["protein_id"], {}).get("psap") is None]
    missing_llphy = [r for r in rows if cached.get(r["protein_id"], {}).get("llphyscore") is None]

    if psap_python and missing_psap:
        try:
            # batch all missing
            psap_map = run_psap_on_rows(missing_psap, psap_python)
            status["psap"] = f"scored_{len(psap_map)}"
        except Exception as e:  # noqa: BLE001
            status["psap"] = f"error:{e}"
    elif not psap_python:
        status["psap"] = "skipped_no_PSAP_PYTHON"

    if llphy_home and llphy_home.exists() and missing_llphy:
        try:
            # score in chunks to keep memory bounded
            chunk = 40
            for i in range(0, len(missing_llphy), chunk):
                part = missing_llphy[i : i + chunk]
                llphy_map.update(run_llphyscore_on_rows(part, llphy_home, llphy_python))
            status["llphyscore"] = f"scored_{len(llphy_map)}"
        except Exception as e:  # noqa: BLE001
            status["llphyscore"] = f"error:{e}"
    elif not llphy_home:
        status["llphyscore"] = "skipped_no_LLPHYSCORE_HOME"

    for r in rows:
        pid = r["protein_id"]
        cached.setdefault(pid, {"psap": None, "llphyscore": None})
        if pid in psap_map:
            cached[pid]["psap"] = psap_map[pid]
        if pid in llphy_map:
            cached[pid]["llphyscore"] = llphy_map[pid]

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["protein_id", "psap", "llphyscore"], delimiter="\t")
        w.writeheader()
        for pid in ids:
            d = cached.get(pid, {})
            w.writerow(
                {
                    "protein_id": pid,
                    "psap": "" if d.get("psap") is None else d["psap"],
                    "llphyscore": "" if d.get("llphyscore") is None else d["llphyscore"],
                }
            )
    # stash status on function attribute for report
    ensure_train_oob_cache.status = status  # type: ignore[attr-defined]
    return cached


# --- panel evaluators ---

def evaluate_split(
    name: str,
    tr_idx: List[int],
    te_idx: List[int],
    X_pl,
    X_px,
    y,
    meta,
    plaac_scores,
    oob: Dict[str, Dict[str, Optional[float]]],
    seed: int = 7,
    *,
    allow_shared_negatives: bool = False,
) -> dict:
    ytr = [y[i] for i in tr_idx]
    yte = [y[i] for i in te_idx]
    if sum(yte) < 1 or (len(yte) - sum(yte)) < 1 or sum(ytr) < 2 or (len(ytr) - sum(ytr)) < 2:
        return {
            "panel": name,
            "skipped": True,
            "reason": "insufficient class counts",
            "n_train": len(tr_idx),
            "n_test": len(te_idx),
            "n_pos_test": int(sum(yte)),
            "n_neg_test": int(len(yte) - sum(yte)),
        }

    # Leakage assert: no test ID in train (positives always; negatives unless shared-neg protocol)
    tr_ids = {meta[i]["id"] for i in tr_idx}
    te_ids = [meta[i]["id"] for i in te_idx]
    if allow_shared_negatives:
        leak = [
            meta[i]["id"]
            for i in te_idx
            if y[i] == 1 and meta[i]["id"] in tr_ids
        ]
    else:
        leak = [tid for tid in te_ids if tid in tr_ids]
    if leak:
        raise RuntimeError(f"LEAKAGE in panel {name}: {leak[:5]}")

    m_pl = fit_plantllps([X_pl[i] for i in tr_idx], ytr, seed=seed)
    m_px = fit_prox([X_px[i] for i in tr_idx], ytr, seed=seed)
    s_log, s_hyb = score_plantllps_rows(m_pl, [X_pl[i] for i in te_idx])
    s_px = score_prox_rows(m_px, [X_px[i] for i in te_idx])
    s_plaac = [plaac_scores[i] for i in te_idx]
    s_psap = [oob.get(meta[i]["id"], {}).get("psap") for i in te_idx]
    s_llphy = [oob.get(meta[i]["id"], {}).get("llphyscore") for i in te_idx]

    aucs = {
        "plantllps_hybrid": safe_auc(yte, s_hyb),
        "plantllps_logistic": safe_auc(yte, s_log),
        "prox_pu_style": safe_auc(yte, s_px),
        "plaac_approx": safe_auc(yte, s_plaac),
        "psap": safe_auc(yte, s_psap),
        "llphyscore": safe_auc(yte, s_llphy),
    }
    coverage = {
        "psap": sum(1 for s in s_psap if s is not None),
        "llphyscore": sum(1 for s in s_llphy if s is not None),
        "n_test": len(te_idx),
    }
    return {
        "panel": name,
        "skipped": False,
        "leakage_control": "train_excludes_all_test_ids",
        "n_train": len(tr_idx),
        "n_test": len(te_idx),
        "n_pos_test": int(sum(yte)),
        "n_neg_test": int(len(yte) - sum(yte)),
        "test_ids": te_ids,
        "auc": aucs,
        "oob_coverage": coverage,
        "per_protein": [
            {
                "protein_id": meta[i]["id"],
                "label": y[i],
                "family": meta[i]["family"],
                "species": meta[i]["species"],
                "plantllps_hybrid": s_hyb[j],
                "plantllps_logistic": s_log[j],
                "prox_pu_style": s_px[j],
                "plaac_approx": s_plaac[j],
                "psap": s_psap[j],
                "llphyscore": s_llphy[j],
            }
            for j, i in enumerate(te_idx)
        ],
    }


def gold_leave_out_panel(
    gold_rows: List[dict],
    ids,
    X_pl,
    X_px,
    y,
    meta,
    plaac_scores,
    oob,
    seqs_by_id: Dict[str, str],
):
    gold_ids = [r["protein_id"] for r in gold_rows]
    gold_set = set(gold_ids)
    id_to_idx = {m["id"]: i for i, m in enumerate(meta)}

    # Train = v13 rows not in gold
    tr_idx = [i for i, m in enumerate(meta) if m["id"] not in gold_set]
    # Test = gold; build feature rows (gold seqs may already be in train matrix)
    te_idx = []
    extra_X_pl, extra_X_px, extra_plaac, extra_meta, extra_y = [], [], [], [], []
    for r in gold_rows:
        pid = r["protein_id"]
        lab = 1 if r["label"].lower().startswith("pos") else 0
        if pid in id_to_idx:
            te_idx.append(id_to_idx[pid])
        else:
            seq = r.get("sequence") or seqs_by_id.get(pid, "")
            feats = extract_features(seq)
            px = extract_prox_features(seq)
            idx = len(y) + len(extra_y)
            extra_X_pl.append(features_to_vector(feats, FEATURE_NAMES))
            extra_X_px.append(prox_features_to_vector(px, PROX_FEATURE_NAMES))
            extra_plaac.append(float(px["plaac_nllr_approx"]))
            extra_meta.append(
                {
                    "id": pid,
                    "family": "gold",
                    "species": "gold",
                    "crop_pool": "gold",
                    "note": r.get("note", ""),
                    "source": "gold",
                }
            )
            extra_y.append(lab)
            te_idx.append(idx)

    # Extend matrices temporarily
    X_pl2 = list(X_pl) + extra_X_pl
    X_px2 = list(X_px) + extra_X_px
    y2 = list(y) + extra_y
    meta2 = list(meta) + extra_meta
    plaac2 = list(plaac_scores) + extra_plaac

    block = evaluate_split(
        "gold_leave_out",
        tr_idx,
        te_idx,
        X_pl2,
        X_px2,
        y2,
        meta2,
        plaac2,
        oob,
    )
    block["protocol"] = (
        "Retrain PlantLLPS + ProX-style on v13 excluding all 28 gold IDs; "
        "score reconstructed 13+15 gold. PSAP/LLPhyScore/PLAAC are OOB."
    )
    block["n_gold_overlap_train"] = sum(1 for g in gold_ids if g in id_to_idx)
    return block


def family_holdout_panels(X_pl, X_px, y, meta, plaac_scores, oob, min_pos=2):
    neg_idx = [i for i, yi in enumerate(y) if yi == 0]
    families = sorted({m["family"] for m in meta})
    folds = []
    for fam in families:
        fam_pos = [i for i, m in enumerate(meta) if m["family"] == fam and y[i] == 1]
        if len(fam_pos) < min_pos:
            folds.append(
                {
                    "panel": f"family_holdout:{fam}",
                    "skipped": True,
                    "reason": f"need ≥{min_pos} positives",
                    "n_pos_family": len(fam_pos),
                }
            )
            continue
        te = list(dict.fromkeys(fam_pos + neg_idx))
        tr = [i for i in range(len(y)) if i not in set(fam_pos)]
        # Note: shared negatives appear in both train and test — documented protocol
        # (same as eval_family_holdout.py). Positives are strictly held out.
        block = evaluate_split(
            f"family_holdout:{fam}",
            tr,
            te,
            X_pl,
            X_px,
            y,
            meta,
            plaac_scores,
            oob,
            allow_shared_negatives=True,
        )
        if not block.get("skipped"):
            block["leakage_control"] = (
                "held_out_family_positives_excluded_from_train; "
                "shared_negatives_in_train_and_test (family_holdout protocol)"
            )
            block["n_pos_family"] = len(fam_pos)
            block["held_out_family"] = fam
        folds.append(block)

    ok = [f for f in folds if not f.get("skipped")]
    summary = {
        "protocol": "holdout_pos_family_plus_shared_negatives",
        "n_families_evaluated": len(ok),
        "mean_auc": {
            m: mean_finite([f["auc"].get(m) for f in ok if "auc" in f]) for m in METHODS
        },
        "folds": folds,
    }
    return summary


def species_transfer_panel(
    name: str,
    train_species: Sequence[str],
    test_species: Sequence[str],
    X_pl,
    X_px,
    y,
    meta,
    plaac_scores,
    oob,
    species_key: str = "species",
):
    train_set = set(train_species)
    test_set = set(test_species)
    te_idx = [i for i, m in enumerate(meta) if m[species_key] in test_set]
    # Train: positives from train species + negatives not in test species
    pos_tr = [
        i
        for i, m in enumerate(meta)
        if y[i] == 1 and m[species_key] in train_set
    ]
    neg_tr = [
        i
        for i, m in enumerate(meta)
        if y[i] == 0 and m[species_key] not in test_set
    ]
    tr_idx = sorted(set(pos_tr + neg_tr))
    n_pos_te = sum(1 for i in te_idx if y[i] == 1)
    n_neg_te = sum(1 for i in te_idx if y[i] == 0)
    if len(te_idx) < 4 or n_pos_te < 1 or n_neg_te < 1:
        return {
            "panel": name,
            "skipped": True,
            "reason": "need ≥4 test rows with both classes",
            "n_train": len(tr_idx),
            "n_test": len(te_idx),
            "n_pos_test": n_pos_te,
            "n_neg_test": n_neg_te,
            "train_species": list(train_species),
            "test_species": list(test_species),
        }
    block = evaluate_split(
        name,
        tr_idx,
        te_idx,
        X_pl,
        X_px,
        y,
        meta,
        plaac_scores,
        oob,
    )
    block["train_species"] = list(train_species)
    block["test_species"] = list(test_species)
    return block


def production_insample_gold(model_path: Path, gold_rows: List[dict]) -> dict:
    if not model_path.exists():
        return {"available": False}
    model = PlantLLPSModel.load(model_path)
    y, s_hyb, s_log = [], [], []
    for r in gold_rows:
        lab = 1 if r["label"].lower().startswith("pos") else 0
        out = model.score(r["sequence"])
        y.append(lab)
        s_hyb.append(out["score_plantllps"])
        s_log.append(out["score_logistic"])
    return {
        "available": True,
        "caveat": "partially_in_sample_gold_ids_overlap_train__NOT_a_fair_transfer_claim",
        "n": len(y),
        "auc": {
            "plantllps_hybrid": safe_auc(y, s_hyb),
            "plantllps_logistic": safe_auc(y, s_log),
        },
    }


def write_scores_tsv(path: Path, per_protein: List[dict]) -> None:
    fields = [
        "protein_id",
        "label",
        "family",
        "species",
        "plantllps_hybrid",
        "plantllps_logistic",
        "prox_pu_style",
        "plaac_approx",
        "psap",
        "llphyscore",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in per_protein:
            w.writerow(r)


def write_leaderboard_md(report: dict, path: Path) -> None:
    gold = report["panels"]["gold_leave_out"]
    fam = report["panels"]["family_holdout"]
    transfers = report["panels"]["species_transfer"]

    def fmt(v):
        if v is None or (isinstance(v, float) and (math.isnan(v) or v != v)):
            return "—"
        return f"{v:.3f}"

    lines = [
        "# Fair benchmark leaderboard (NAR+)",
        "",
        "Leakage-controlled comparison of **PlantLLPS** vs native **PSAP**, **LLPhyScore**,",
        "**PLAAC-approx**, and a **ProX/PU-style** feature logistic on the same panels.",
        "Generated by `scripts/eval_fair_benchmark.py` → `data/benchmarks/fair_comparison_metrics.json`.",
        "",
        "## Leakage policy",
        "",
        "- Supervised methods **retrain excluding every test ID** (gold leave-out; family",
        "  positives held out; species-transfer test species held out).",
        "- PSAP / LLPhyScore / PLAAC-approx are **out-of-the-box** (no plant-label fit).",
        "- Production-model gold AUCs are **in-sample** and listed only as a caveat — not on the leaderboard.",
        "- ProX/PU-style = stdlib logistic on IDR/LCR/hydropathy/FCR/PLAAC/aromatic/κ features;",
        "  **not** the published Plant PS proteins LightGBM/PU model.",
        "",
        "## Gold panel (leave-gold-out)",
        "",
        f"n = {gold.get('n_test')} ({gold.get('n_pos_test')} pos / {gold.get('n_neg_test')} neg); "
        f"train n = {gold.get('n_train')}.",
        "",
        "| Method | AUC | Notes |",
        "|--------|----:|-------|",
    ]
    notes = {
        "plantllps_hybrid": "retrained hybrid; fair",
        "plantllps_logistic": "retrained logistic; fair",
        "prox_pu_style": "ProX-like features + logistic; fair",
        "plaac_approx": "stdlib PLAAC-approx OOB",
        "psap": "native PSAP OOB",
        "llphyscore": "LLPhyScore 8-feat %ile OOB",
    }
    order = [
        "plantllps_hybrid",
        "plantllps_logistic",
        "prox_pu_style",
        "llphyscore",
        "psap",
        "plaac_approx",
    ]
    auc_g = gold.get("auc") or {}
    # sort by AUC desc for display ranking
    ranked = sorted(order, key=lambda m: (auc_g.get(m) is None, -(auc_g.get(m) or -1)))
    for m in ranked:
        lines.append(f"| {m} | {fmt(auc_g.get(m))} | {notes[m]} |")

    ins = report.get("caveat_in_sample_production_gold") or {}
    if ins.get("available"):
        lines += [
            "",
            "### Caveat — production model on gold (in-sample; not fair)",
            "",
            f"hybrid {fmt((ins.get('auc') or {}).get('plantllps_hybrid'))}, "
            f"logistic {fmt((ins.get('auc') or {}).get('plantllps_logistic'))}. "
            "Do **not** cite as independent transfer.",
        ]

    lines += [
        "",
        "## Family hold-out (mean over evaluated families)",
        "",
        f"Protocol: hold out ≥2 positives from one family; score vs shared negatives. "
        f"Families evaluated: **{fam.get('n_families_evaluated')}**.",
        "",
        "| Method | Mean AUC |",
        "|--------|---------:|",
    ]
    mean_auc = fam.get("mean_auc") or {}
    fam_ranked = sorted(
        order, key=lambda m: (mean_auc.get(m) is None, -(mean_auc.get(m) or -1))
    )
    for m in fam_ranked:
        lines.append(f"| {m} | {fmt(mean_auc.get(m))} |")

    lines += [
        "",
        "## Cross-species transfer",
        "",
        "| Transfer | n_test (pos/neg) | PlantLLPS hyb | PlantLLPS log | ProX-style | PLAAC | PSAP | LLPhy |",
        "|----------|-----------------:|-------------:|--------------:|-----------:|------:|-----:|------:|",
    ]
    for t in transfers:
        if t.get("skipped"):
            lines.append(
                f"| {t.get('panel')} | skipped ({t.get('reason', '')}) | — | — | — | — | — | — |"
            )
            continue
        a = t.get("auc") or {}
        lines.append(
            f"| {t.get('panel')} | {t.get('n_test')} ({t.get('n_pos_test')}/{t.get('n_neg_test')}) | "
            f"{fmt(a.get('plantllps_hybrid'))} | {fmt(a.get('plantllps_logistic'))} | "
            f"{fmt(a.get('prox_pu_style'))} | {fmt(a.get('plaac_approx'))} | "
            f"{fmt(a.get('psap'))} | {fmt(a.get('llphyscore'))} |"
        )

    lines += [
        "",
        "## Honesty",
        "",
        report.get("honesty", ""),
        "",
        f"Metrics JSON: `data/benchmarks/fair_comparison_metrics.json`  ",
        f"Gold scores TSV: `data/benchmarks/fair_comparison_scores.tsv`  ",
        f"OOB cache: `data/benchmarks/train_v13_oob_comparator_scores.tsv`",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train-csv", type=Path, default=DEFAULT_TRAIN)
    ap.add_argument("--gold-csv", type=Path, default=DEFAULT_GOLD)
    ap.add_argument("--model", type=Path, default=ROOT / "data" / "model_formal.json")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--scores-tsv", type=Path, default=DEFAULT_SCORES)
    ap.add_argument("--oob-cache", type=Path, default=DEFAULT_CACHE)
    ap.add_argument("--psap-python", default=os.environ.get("PSAP_PYTHON"))
    ap.add_argument(
        "--llphyscore-home",
        default=os.environ.get(
            "LLPHYSCORE_HOME",
            "/tmp/LLPhyScore/standalone_package/LLPhyScore",
        ),
    )
    ap.add_argument(
        "--llphyscore-python",
        default=os.environ.get("LLPHYSCORE_PYTHON"),
    )
    ap.add_argument("--skip-oob-run", action="store_true", help="Use cache/deposits only")
    ap.add_argument("--force-oob", action="store_true")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    ids, seqs, labels, notes, sources, kingdoms = load_train(args.train_csv)
    X_pl, X_px, y, meta, plaac_scores = build_matrices(
        ids, seqs, labels, notes, sources, kingdoms
    )
    gold_rows = load_gold(args.gold_csv)
    seqs_by_id = dict(zip(ids, seqs))

    # OOB scores for train (+ gold via deposit merge)
    if args.skip_oob_run:
        oob = {}
        if args.oob_cache.exists():
            with args.oob_cache.open(newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    oob[row["protein_id"]] = {
                        "psap": float(row["psap"]) if row.get("psap") else None,
                        "llphyscore": float(row["llphyscore"]) if row.get("llphyscore") else None,
                    }
        for pid, d in load_deposited_gold_oob().items():
            oob.setdefault(pid, {"psap": None, "llphyscore": None})
            oob[pid].update(d)
        oob_status = {"mode": "cache_or_deposit_only"}
    else:
        llphy_home = Path(args.llphyscore_home) if args.llphyscore_home else None
        oob = ensure_train_oob_cache(
            ids,
            seqs,
            notes,
            args.oob_cache,
            None if args.skip_oob_run else args.psap_python,
            llphy_home,
            args.llphyscore_python,
            force=args.force_oob,
        )
        oob_status = getattr(ensure_train_oob_cache, "status", {})

    # Ensure gold OOB present
    for pid, d in load_deposited_gold_oob().items():
        oob.setdefault(pid, {"psap": None, "llphyscore": None})
        for k, v in d.items():
            if oob[pid].get(k) is None:
                oob[pid][k] = v

    gold_panel = gold_leave_out_panel(
        gold_rows, ids, X_pl, X_px, y, meta, plaac_scores, oob, seqs_by_id
    )
    fam_panel = family_holdout_panels(X_pl, X_px, y, meta, plaac_scores, oob)

    # Expanded species transfers
    transfers = [
        species_transfer_panel(
            "Ath→rice", ["arabidopsis"], ["rice"], X_pl, X_px, y, meta, plaac_scores, oob
        ),
        species_transfer_panel(
            "Ath→maize", ["arabidopsis"], ["maize"], X_pl, X_px, y, meta, plaac_scores, oob
        ),
        species_transfer_panel(
            "Ath→soybean",
            ["arabidopsis"],
            ["soybean"],
            X_pl,
            X_px,
            y,
            meta,
            plaac_scores,
            oob,
        ),
        species_transfer_panel(
            "Ath→cereal(rice+maize)",
            ["arabidopsis"],
            ["rice", "maize"],
            X_pl,
            X_px,
            y,
            meta,
            plaac_scores,
            oob,
        ),
        species_transfer_panel(
            "Ath→non_ath_crop",
            ["arabidopsis"],
            ["rice", "maize", "soybean", "tomato", "grape", "apple_peach_citrus", "other_plant"],
            X_pl,
            X_px,
            y,
            meta,
            plaac_scores,
            oob,
        ),
        species_transfer_panel(
            "non_ath_crop→Ath",
            ["rice", "maize", "soybean", "tomato", "grape", "apple_peach_citrus", "other_plant"],
            ["arabidopsis"],
            X_pl,
            X_px,
            y,
            meta,
            plaac_scores,
            oob,
        ),
        species_transfer_panel(
            "cereal→Ath",
            ["rice", "maize"],
            ["arabidopsis"],
            X_pl,
            X_px,
            y,
            meta,
            plaac_scores,
            oob,
        ),
    ]

    # Counts
    sp_counts = defaultdict(lambda: {"pos": 0, "neg": 0})
    pool_counts = defaultdict(lambda: {"pos": 0, "neg": 0})
    for m, yi in zip(meta, y):
        sp_counts[m["species"]]["pos" if yi else "neg"] += 1
        pool_counts[m["crop_pool"]]["pos" if yi else "neg"] += 1

    report = {
        "version": 1,
        "suite": "nar_plus_fair_comparison",
        "train_csv": str(args.train_csv.relative_to(ROOT)),
        "gold_csv": str(args.gold_csv.relative_to(ROOT)),
        "n_train": len(y),
        "n_pos": int(sum(y)),
        "n_neg": int(len(y) - sum(y)),
        "species_counts": dict(sp_counts),
        "crop_pool_counts": dict(pool_counts),
        "methods": {
            "plantllps_hybrid": "PlantLLPS 0.5·logistic+0.5·companion; leakage-safe retrain",
            "plantllps_logistic": "PlantLLPS SGD logistic on FEATURE_NAMES; leakage-safe retrain",
            "prox_pu_style": (
                "stdlib logistic on ProX/PU-like features "
                f"{PROX_FEATURE_NAMES}; not published LightGBM/PU"
            ),
            "plaac_approx": "stdlib PLAAC-approx NLLR; OOB",
            "psap": "native PSAP 1.0.7; OOB",
            "llphyscore": "LLPhyScore standalone 8-feature-sum human+PDB %ile; OOB",
        },
        "oob_status": oob_status,
        "panels": {
            "gold_leave_out": {
                k: v for k, v in gold_panel.items() if k != "per_protein"
            },
            "family_holdout": fam_panel,
            "species_transfer": [
                {k: v for k, v in t.items() if k != "per_protein"} for t in transfers
            ],
        },
        "caveat_in_sample_production_gold": production_insample_gold(args.model, gold_rows),
        "honesty": (
            "Fair claims use leave-gold-out / family hold-out / species-transfer AUCs only. "
            "Production-model gold AUCs remain partially in-sample and are caveated. "
            "Ath→rice and Ath→maize n are small — treat as exploratory. "
            "ProX/PU-style is a feature-family baseline, not the Plant PS proteins web model. "
            "PLAAC-approx is not the Java PLAAC HMM. Family hold-out shares negatives across "
            "train/test by design (positives strictly held out)."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2), encoding="utf-8")

    if gold_panel.get("per_protein"):
        write_scores_tsv(args.scores_tsv, gold_panel["per_protein"])

    # Family fold TSV (compact)
    fam_tsv = ROOT / "data" / "benchmarks" / "fair_comparison_family.tsv"
    with fam_tsv.open("w", newline="", encoding="utf-8") as fh:
        fields = ["held_out_family", "skipped", "n_pos_test", "n_neg_test"] + METHODS
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for f in fam_panel["folds"]:
            row = {
                "held_out_family": f.get("held_out_family")
                or f.get("panel", "").split(":")[-1],
                "skipped": f.get("skipped", False),
                "n_pos_test": f.get("n_pos_test", f.get("n_pos_family", "")),
                "n_neg_test": f.get("n_neg_test", ""),
            }
            for m in METHODS:
                row[m] = (f.get("auc") or {}).get(m, "")
            w.writerow(row)

    readme = ROOT / "data" / "benchmarks" / "FAIR_COMPARISON_README.md"
    g_auc = gold_panel.get("auc") or {}
    readme.write_text(
        f"""# Fair comparison suite (NAR+)

Generated by `scripts/eval_fair_benchmark.py`.

## Headline (leakage-controlled)

| Panel | PlantLLPS hybrid | ProX-style | LLPhyScore | PSAP | PLAAC-approx |
|-------|-----------------:|-----------:|-----------:|-----:|-------------:|
| Gold leave-out | {g_auc.get('plantllps_hybrid')} | {g_auc.get('prox_pu_style')} | {g_auc.get('llphyscore')} | {g_auc.get('psap')} | {g_auc.get('plaac_approx')} |
| Family hold-out mean | {(fam_panel.get('mean_auc') or {}).get('plantllps_hybrid')} | {(fam_panel.get('mean_auc') or {}).get('prox_pu_style')} | {(fam_panel.get('mean_auc') or {}).get('llphyscore')} | {(fam_panel.get('mean_auc') or {}).get('psap')} | {(fam_panel.get('mean_auc') or {}).get('plaac_approx')} |

See `docs/BENCHMARK_LEADERBOARD.md` for the full ranked table + species transfers.

## Reproduce

```bash
export PSAP_PYTHON=/path/to/psap39/bin/python
export LLPHYSCORE_HOME=/path/to/LLPhyScore/standalone_package/LLPhyScore
export LLPHYSCORE_PYTHON=/path/to/llphy39/bin/python   # optional
python3 scripts/eval_fair_benchmark.py
```

`--skip-oob-run` reuses `train_v13_oob_comparator_scores.tsv` + deposited gold scores.
""",
        encoding="utf-8",
    )

    write_leaderboard_md(report, ROOT / "docs" / "BENCHMARK_LEADERBOARD.md")

    print(
        "Gold leave-out AUCs:",
        json.dumps(g_auc, indent=2),
    )
    print(
        "Family mean AUCs:",
        json.dumps(fam_panel.get("mean_auc"), indent=2),
    )
    print(f"Wrote {args.out}")
    print(f"Wrote {args.scores_tsv}")
    print(f"Wrote {ROOT / 'docs' / 'BENCHMARK_LEADERBOARD.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
