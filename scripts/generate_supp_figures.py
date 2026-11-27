#!/usr/bin/env python3
"""Generate Supplemental Figures S1–S5 (PNG + PDF) from deposited metrics.

Reproducibility
---------------
- Reads deposited JSONs/TSVs/SQLite under ``data/`` (no network).
- Stratified CV ROC (S1) uses the same seed/k as production metrics
  (``seed=7``, ``k=5``, ``data/train_formal_v13.csv``).
- Matplotlib Agg backend; pin ``matplotlib==3.10.8`` or compatible 3.10+/3.11
  (see ``requirements-figures.txt``).

Outputs (default)::

    docs/figures/supp/Figure_S{1..5}.{png,pdf}
    docs/figures/supp/LEGENDS.md

Usage::

    python3 scripts/generate_supp_figures.py
    python3 scripts/generate_supp_figures.py --out-dir docs/figures/supp
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from plantllps.features import FEATURE_NAMES, extract_features, features_to_vector  # noqa: E402
from plantllps.metrics import roc_auc, stratified_kfold_indices  # noqa: E402
from plantllps.model import PlantLLPSModel  # noqa: E402

# Visual language shared with feature_ablation_panel.svg
C_TEAL = "#2a6f6f"
C_TERRACOTTA = "#c45c26"
C_INK = "#1a1a1a"
C_MUTED = "#555555"
C_GRID = "#d8d4cc"
C_BG = "#fbfaf7"
C_TIER = {"A": "#0B6E4F", "B": "#C45C26", "C": "#6B7C85", "D": "#B8B2A6"}
METHOD_COLORS = {
    "plantllps_hybrid": "#0B6E4F",
    "plantllps_logistic": "#2a6f6f",
    "prox_pu_style": "#8B6914",
    "llphyscore": "#3D5A80",
    "psap": "#9A3412",
    "plaac_approx": "#6B7280",
}
METHOD_LABELS = {
    "plantllps_hybrid": "PlantLLPS hybrid",
    "plantllps_logistic": "PlantLLPS logistic",
    "prox_pu_style": "ProX-style",
    "llphyscore": "LLPhyScore",
    "psap": "PSAP",
    "plaac_approx": "PLAAC-approx",
}

SEED = 7
K_FOLDS = 5
DPI = 300


def _style_ax(ax) -> None:
    ax.set_facecolor(C_BG)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#888888")
    ax.spines["bottom"].set_color("#888888")
    ax.tick_params(colors=C_MUTED, labelsize=9)


def _savefig(fig, out_dir: Path, stem: str) -> Tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{stem}.png"
    pdf = out_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=DPI, bbox_inches="tight", facecolor=C_BG)
    fig.savefig(pdf, bbox_inches="tight", facecolor=C_BG)
    plt.close(fig)
    return png, pdf


def _roc_curve(y: Sequence[int], scores: Sequence[float]) -> Tuple[List[float], List[float], float]:
    """Return FPR, TPR, AUC from scores (higher = more positive). Stdlib only."""
    pairs = sorted(zip(scores, y), key=lambda t: -t[0])
    n_pos = sum(1 for yi in y if yi == 1)
    n_neg = len(y) - n_pos
    if n_pos == 0 or n_neg == 0:
        return [0.0, 1.0], [0.0, 1.0], float("nan")
    tp = fp = 0
    fpr = [0.0]
    tpr = [0.0]
    prev = None
    for score, lab in pairs:
        if prev is not None and score != prev:
            fpr.append(fp / n_neg)
            tpr.append(tp / n_pos)
        if lab == 1:
            tp += 1
        else:
            fp += 1
        prev = score
    fpr.append(fp / n_neg)
    tpr.append(tp / n_pos)
    # force endpoints
    if fpr[0] != 0.0 or tpr[0] != 0.0:
        fpr = [0.0] + fpr
        tpr = [0.0] + tpr
    if fpr[-1] != 1.0 or tpr[-1] != 1.0:
        fpr = fpr + [1.0]
        tpr = tpr + [1.0]
    auc = roc_auc(y, scores)
    return fpr, tpr, auc


def _load_train(path: Path):
    ids, seqs, labels = [], [], []
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            ids.append(row["protein_id"])
            seqs.append(row["sequence"].strip())
            lab = row["label"].strip().lower()
            labels.append(1 if lab.startswith("pos") else 0)
    return ids, seqs, labels


def _oof_cv_scores(
    train_csv: Path, *, k: int = K_FOLDS, seed: int = SEED
) -> Dict[str, List[float]]:
    """Out-of-fold logistic / companion / hybrid scores (production CV protocol)."""
    _ids, seqs, y = _load_train(train_csv)
    X = [features_to_vector(extract_features(s), FEATURE_NAMES) for s in seqs]
    oof_log = [float("nan")] * len(y)
    oof_comp = [float("nan")] * len(y)
    oof_hyb = [float("nan")] * len(y)
    fold_auc_log: List[float] = []
    for train_idx, test_idx in stratified_kfold_indices(y, k=k, seed=seed):
        Xtr = [X[i] for i in train_idx]
        ytr = [y[i] for i in train_idx]
        m = PlantLLPSModel()
        m.fit_logistic(Xtr, ytr, seed=seed)
        y_te = [y[i] for i in test_idx]
        s_log = []
        for i in test_idx:
            feats = {n: v for n, v in zip(FEATURE_NAMES, X[i])}
            sl = m.logistic_score(feats)
            sc = m.companion_score(feats)
            sh = 0.5 * sl + 0.5 * sc
            oof_log[i] = sl
            oof_comp[i] = sc
            oof_hyb[i] = sh
            s_log.append(sl)
        fold_auc_log.append(roc_auc(y_te, s_log))
    return {
        "y": y,
        "logistic": oof_log,
        "companion": oof_comp,
        "hybrid": oof_hyb,
        "fold_auc_logistic": fold_auc_log,
    }


# ---------------------------------------------------------------------------
# S1 — stratified CV + gold ROC (in-sample caveat)
# ---------------------------------------------------------------------------


def figure_s1(out_dir: Path) -> Tuple[Path, Path]:
    metrics = json.loads((ROOT / "data" / "model_formal.metrics.v13.json").read_text())
    gold_rows = list(
        csv.DictReader(
            (ROOT / "data" / "benchmarks" / "gold_standard_scores.tsv").open(),
            delimiter="\t",
        )
    )
    oof = _oof_cv_scores(ROOT / "data" / "train_formal_v13.csv")

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.6), facecolor=C_BG)

    # Panel A — OOF ROC
    ax = axes[0]
    _style_ax(ax)
    for key, label, color, ls in (
        ("logistic", "Logistic OOF", C_TEAL, "-"),
        ("hybrid", "Hybrid OOF", C_TERRACOTTA, "-"),
        ("companion", "Companion OOF", "#6B7C85", "--"),
    ):
        fpr, tpr, auc = _roc_curve(oof["y"], oof[key])
        ax.plot(fpr, tpr, color=color, lw=2.0, ls=ls, label=f"{label} (AUC={auc:.3f})")
    ax.plot([0, 1], [0, 1], color="#aaaaaa", lw=1, ls=":")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False positive rate", fontsize=10)
    ax.set_ylabel("True positive rate", fontsize=10)
    ax.set_title(
        f"A  Stratified 5-fold CV (v13, n={len(oof['y'])})\n"
        f"mean logistic AUC={metrics['mean_auc_logistic']:.3f}; "
        f"OOF={metrics['oof_auc_logistic']:.3f}",
        fontsize=11,
        color=C_INK,
        loc="left",
        pad=8,
    )
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    # fold AUC strip
    folds = metrics["fold_auc_logistic"]
    ax.text(
        0.02,
        0.02,
        "fold AUC: " + ", ".join(f"{a:.3f}" for a in folds),
        transform=ax.transAxes,
        fontsize=7,
        color=C_MUTED,
        va="bottom",
    )

    # Panel B — gold ROC (production = partially in-sample)
    ax = axes[1]
    _style_ax(ax)
    y_gold = [1 if r["label"].lower().startswith("pos") else 0 for r in gold_rows]
    for col, label, color, ls in (
        ("score_plantllps", "Hybrid (production)", C_TERRACOTTA, "-"),
        ("score_logistic", "Logistic (production)", C_TEAL, "-"),
        ("score_companion", "Companion", "#6B7C85", "--"),
    ):
        scores = [float(r[col]) for r in gold_rows]
        fpr, tpr, auc = _roc_curve(y_gold, scores)
        ax.plot(fpr, tpr, color=color, lw=2.0, ls=ls, label=f"{label} (AUC={auc:.3f})")
    # OOB comparators when present
    if gold_rows and gold_rows[0].get("psap_score"):
        fpr, tpr, auc = _roc_curve(y_gold, [float(r["psap_score"]) for r in gold_rows])
        ax.plot(fpr, tpr, color=METHOD_COLORS["psap"], lw=1.5, ls="--", label=f"PSAP OOB ({auc:.3f})")
    if gold_rows and gold_rows[0].get("llphyscore_percentile"):
        fpr, tpr, auc = _roc_curve(
            y_gold, [float(r["llphyscore_percentile"]) for r in gold_rows]
        )
        ax.plot(
            fpr,
            tpr,
            color=METHOD_COLORS["llphyscore"],
            lw=1.5,
            ls="--",
            label=f"LLPhy %ile OOB ({auc:.3f})",
        )
    ax.plot([0, 1], [0, 1], color="#aaaaaa", lw=1, ls=":")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False positive rate", fontsize=10)
    ax.set_ylabel("True positive rate", fontsize=10)
    n_in = sum(1 for r in gold_rows if r.get("in_plantllps_train", "").lower() == "true")
    ax.set_title(
        f"B  Gold-standard ROC (n={len(gold_rows)}; {n_in} overlap train)\n"
        "IN-SAMPLE CAVEAT — not independent transfer",
        fontsize=11,
        color=C_TERRACOTTA,
        loc="left",
        pad=8,
    )
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    # caveat badge
    ax.text(
        0.98,
        0.55,
        "⚠ partially\nin-sample",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=9,
        color=C_TERRACOTTA,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#fde8dc", edgecolor=C_TERRACOTTA, lw=1),
    )

    fig.suptitle(
        "Supplemental Figure S1. PlantLLPS v13 stratified CV and gold ROC",
        fontsize=12,
        fontweight="bold",
        color=C_INK,
        y=1.02,
    )
    fig.tight_layout()
    return _savefig(fig, out_dir, "Figure_S1")


# ---------------------------------------------------------------------------
# S2 — fair leave-gold-out / family / Ath→non_ath_crop
# ---------------------------------------------------------------------------


def figure_s2(out_dir: Path) -> Tuple[Path, Path]:
    meta = json.loads((ROOT / "data" / "benchmarks" / "fair_comparison_metrics.json").read_text())
    scores = list(
        csv.DictReader(
            (ROOT / "data" / "benchmarks" / "fair_comparison_scores.tsv").open(),
            delimiter="\t",
        )
    )
    methods = [
        "plantllps_hybrid",
        "plantllps_logistic",
        "prox_pu_style",
        "llphyscore",
        "psap",
        "plaac_approx",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2), facecolor=C_BG)

    # A — leave-gold-out ROC
    ax = axes[0, 0]
    _style_ax(ax)
    y = [int(r["label"]) for r in scores]
    for m in methods:
        raw = [float(r[m]) for r in scores]
        # LLPhy is percentile (higher better); others higher-better already
        fpr, tpr, auc = _roc_curve(y, raw)
        ax.plot(
            fpr,
            tpr,
            color=METHOD_COLORS[m],
            lw=2.0 if m.startswith("plantllps") else 1.5,
            label=f"{METHOD_LABELS[m]} ({auc:.3f})",
        )
    ax.plot([0, 1], [0, 1], color="#aaaaaa", lw=1, ls=":")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    g = meta["panels"]["gold_leave_out"]
    ax.set_title(
        f"A  Leave-gold-out ROC (n={g['n_test']}; train={g['n_train']})",
        fontsize=11,
        loc="left",
        color=C_INK,
    )
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    ax.legend(frameon=False, fontsize=7, loc="lower right")

    def _bar_panel(ax, auc_map: Mapping[str, float], title: str, order: Sequence[str]):
        _style_ax(ax)
        vals = [float(auc_map[m]) for m in order]
        colors = [METHOD_COLORS[m] for m in order]
        ypos = list(range(len(order)))[::-1]
        bars = ax.barh(ypos, vals, color=colors, height=0.65, edgecolor="none")
        ax.set_yticks(ypos)
        ax.set_yticklabels([METHOD_LABELS[m] for m in order], fontsize=8)
        ax.set_xlim(0, 1.05)
        ax.axvline(0.5, color="#aaaaaa", lw=1, ls=":")
        ax.set_xlabel("AUC")
        ax.set_title(title, fontsize=11, loc="left", color=C_INK)
        for bar, v in zip(bars, vals):
            ax.text(
                min(v + 0.02, 1.02),
                bar.get_y() + bar.get_height() / 2,
                f"{v:.3f}",
                va="center",
                fontsize=8,
                color=C_INK,
            )

    order = methods
    _bar_panel(
        axes[0, 1],
        g["auc"],
        "B  Leave-gold-out AUC summary",
        order,
    )

    fam = meta["panels"]["family_holdout"]["mean_auc"]
    _bar_panel(
        axes[1, 0],
        fam,
        f"C  Family hold-out mean AUC ({meta['panels']['family_holdout']['n_families_evaluated']} families)",
        order,
    )

    # Ath→non_ath_crop
    transfer = next(
        (
            item
            for item in meta["panels"]["species_transfer"]
            if item.get("panel") == "Ath→non_ath_crop" and not item.get("skipped")
        ),
        None,
    )

    ax = axes[1, 1]
    if transfer and transfer.get("auc"):
        n_te = transfer.get("n_test", 38)
        _bar_panel(
            ax,
            transfer["auc"],
            f"D  Ath→non_ath_crop transfer (n={n_te})",
            order,
        )
    else:
        _style_ax(ax)
        ax.text(0.5, 0.5, "Ath→non_ath_crop panel missing", ha="center", va="center")

    fig.suptitle(
        "Supplemental Figure S2. Fair benchmark vs ProX-style / PSAP / LLPhy / PLAAC",
        fontsize=12,
        fontweight="bold",
        color=C_INK,
        y=1.01,
    )
    fig.tight_layout()
    return _savefig(fig, out_dir, "Figure_S2")


# ---------------------------------------------------------------------------
# S3 — feature ablation + F6HET4 SHAP grammar
# ---------------------------------------------------------------------------


def figure_s3(out_dir: Path) -> Tuple[Path, Path]:
    abl = json.loads((ROOT / "data" / "benchmarks" / "feature_ablation_metrics.json").read_text())
    shap = json.loads((ROOT / "data" / "explainability" / "F6HET4.json").read_text())

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.0), facecolor=C_BG)

    # A — leave-one-group-out ΔAUC
    ax = axes[0]
    _style_ax(ax)
    rows = sorted(
        abl["leave_one_group_out"],
        key=lambda r: r["delta_mean_auc_logistic"],
    )
    groups = [r["group"] for r in rows]
    deltas = [r["delta_mean_auc_logistic"] for r in rows]
    colors = [C_TERRACOTTA if d < 0 else C_TEAL for d in deltas]
    ypos = list(range(len(groups)))
    ax.barh(ypos, deltas, color=colors, height=0.7, edgecolor="none")
    ax.axvline(0, color="#888888", lw=1)
    ax.set_yticks(ypos)
    ax.set_yticklabels(groups, fontsize=9)
    ax.set_xlabel("Δ mean logistic AUC (ablated − full)")
    full = abl["full_model"]["mean_auc_logistic"]
    ax.set_title(
        f"A  Leave-one-group-out ablation\nfull-model mean logistic AUC = {full:.3f}",
        fontsize=11,
        loc="left",
        color=C_INK,
    )
    for y, d in zip(ypos, deltas):
        ax.text(
            d + (0.0015 if d >= 0 else -0.0015),
            y,
            f"{d:+.3f}",
            va="center",
            ha="left" if d >= 0 else "right",
            fontsize=8,
            color=C_INK,
        )

    # B — F6HET4 grammar SHAP bars
    ax = axes[1]
    _style_ax(ax)
    axes_g = shap["grammar"]["axes"]
    order = ["ser_rich", "charge_block", "idr_disorder"]
    labels = {
        "ser_rich": "Ser-rich / MAP65",
        "charge_block": "Charge-block",
        "idr_disorder": "IDR / disorder",
    }
    vals = [axes_g[k]["logistic_shap_prob"] for k in order]
    colors = [C_TERRACOTTA if v > 0 else C_TEAL for v in vals]
    # IDR negative → muted teal; positive terracotta
    colors = []
    for v in vals:
        if v > 0:
            colors.append(C_TERRACOTTA)
        else:
            colors.append("#6B7C85")
    ypos = list(range(len(order)))[::-1]
    ax.barh(ypos, vals, color=colors, height=0.55, edgecolor="none")
    ax.axvline(0, color="#888888", lw=1)
    ax.set_yticks(ypos)
    ax.set_yticklabels([labels[k] for k in order], fontsize=10)
    ax.set_xlabel("Exact logistic SHAP (probability-local)")
    ax.set_title(
        f"B  F6HET4 (VvMAP65-1) grammar SHAP\n"
        f"Tier {shap['tier']}; PlantLLPS={shap['score_plantllps']:.3f}; "
        f"logistic={shap['score_logistic']:.3f}",
        fontsize=11,
        loc="left",
        color=C_INK,
    )
    for y, v in zip(ypos, vals):
        ax.text(
            v + (0.008 if v >= 0 else -0.008),
            y,
            f"{v:+.3f}",
            va="center",
            ha="left" if v >= 0 else "right",
            fontsize=9,
            color=C_INK,
            fontweight="bold",
        )
    ax.text(
        0.98,
        0.05,
        "Exact linear SHAP in logit space\n(φᵢ = wᵢ·zᵢ); not KernelSHAP",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color=C_MUTED,
    )

    fig.suptitle(
        "Supplemental Figure S3. Feature ablation and F6HET4 SHAP grammar",
        fontsize=12,
        fontweight="bold",
        color=C_INK,
        y=1.02,
    )
    fig.tight_layout()
    return _savefig(fig, out_dir, "Figure_S3")


# ---------------------------------------------------------------------------
# S4 — proteome-wide score distribution + Tier A/B + MAP anchors
# ---------------------------------------------------------------------------


def figure_s4(out_dir: Path) -> Tuple[Path, Path]:
    db = ROOT / "data" / "grape_condensatome.sqlite"
    conn = sqlite3.connect(str(db))
    rows = conn.execute("SELECT protein_id, score, tier, gene_name, protein_name FROM proteins").fetchall()
    conn.close()
    summary = json.loads((ROOT / "data" / "proteome_summary.json").read_text())
    hooks = json.loads((ROOT / "data" / "cross_crop_hooks.json").read_text())

    scores_by_tier: Dict[str, List[float]] = {t: [] for t in "ABCD"}
    score_map: Dict[str, Tuple[float, str, str]] = {}
    for pid, score, tier, gene, pname in rows:
        scores_by_tier.setdefault(tier, []).append(float(score))
        score_map[pid] = (float(score), tier, gene or "")

    anchors = []
    for h in hooks.get("hooks", []):
        for gid in h.get("grape", []):
            if gid in score_map:
                sc, tier, gene = score_map[gid]
                anchors.append(
                    {
                        "id": gid,
                        "score": sc,
                        "tier": tier,
                        "gene": gene,
                        "family": h.get("family", ""),
                        "label": f"{gid} ({h.get('family', '')})",
                    }
                )
    # ensure F6HET4 present
    if not any(a["id"] == "F6HET4" for a in anchors) and "F6HET4" in score_map:
        sc, tier, gene = score_map["F6HET4"]
        anchors.append(
            {
                "id": "F6HET4",
                "score": sc,
                "tier": tier,
                "gene": gene,
                "family": "MAP65-1",
                "label": "F6HET4 (MAP65-1)",
            }
        )

    fig, axes = plt.subplots(
        2, 1, figsize=(10.0, 7.2), facecolor=C_BG, gridspec_kw={"height_ratios": [3, 1.1]}
    )

    ax = axes[0]
    _style_ax(ax)
    bins = [i / 50 for i in range(51)]
    bottom = [0.0] * (len(bins) - 1)
    # stacked hist by tier D→A so A sits on top
    for tier in ("D", "C", "B", "A"):
        vals = scores_by_tier.get(tier, [])
        counts, edges = _histogram(vals, bins)
        ax.bar(
            [(edges[i] + edges[i + 1]) / 2 for i in range(len(counts))],
            counts,
            width=(edges[1] - edges[0]) * 0.95,
            bottom=bottom,
            color=C_TIER[tier],
            edgecolor="none",
            label=f"Tier {tier} (n={len(vals):,})",
            alpha=0.9,
        )
        bottom = [b + c for b, c in zip(bottom, counts)]

    # tier cutoffs
    cut_a = summary["tier_cutoffs"]["score_A"]
    cut_c = summary["tier_cutoffs"]["score_C"]
    ax.axvline(cut_a, color=C_TIER["A"], ls="--", lw=1.2, alpha=0.8)
    ax.axvline(cut_c, color=C_TIER["C"], ls="--", lw=1.2, alpha=0.8)
    ymax = max(bottom) if bottom else 1
    ax.text(cut_a, ymax * 0.92, f" score_A={cut_a:.3f}", color=C_TIER["A"], fontsize=8, va="top")
    ax.text(cut_c, ymax * 0.92, f" score_C={cut_c:.3f}", color=C_MUTED, fontsize=8, va="top")

    # MAP / hook anchors
    for i, a in enumerate(anchors):
        ax.axvline(a["score"], color=C_INK, ls="-", lw=1.0, alpha=0.55)
        ax.annotate(
            a["label"],
            xy=(a["score"], ymax * (0.55 - 0.08 * i)),
            xytext=(a["score"] + 0.04, ymax * (0.72 - 0.12 * i)),
            fontsize=8,
            color=C_INK,
            arrowprops=dict(arrowstyle="->", color=C_INK, lw=0.8),
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#cccccc", lw=0.5),
        )

    ax.set_xlim(0, 1)
    ax.set_xlabel("PlantLLPS hybrid score")
    ax.set_ylabel("Protein count")
    ax.set_title(
        f"A  Proteome-wide score distribution (UP000009183, n={summary['n_scored']:,})",
        fontsize=11,
        loc="left",
        color=C_INK,
    )
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    # B — tier counts bar
    ax = axes[1]
    _style_ax(ax)
    tiers = ["A", "B", "C", "D"]
    counts = [summary["tier_counts"][t] for t in tiers]
    bars = ax.bar(tiers, counts, color=[C_TIER[t] for t in tiers], width=0.6, edgecolor="none")
    ax.set_ylabel("Count")
    ax.set_title("B  Tier counts", fontsize=11, loc="left", color=C_INK)
    for bar, c in zip(bars, counts):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            c + 80,
            f"{c:,}",
            ha="center",
            va="bottom",
            fontsize=9,
            color=C_INK,
        )
    # annotate MAP65-1
    if anchors:
        a0 = next((a for a in anchors if a["id"] == "F6HET4"), anchors[0])
        ax.text(
            0.98,
            0.85,
            f"MAP anchor: {a0['id']} Tier {a0['tier']} score={a0['score']:.3f}",
            transform=ax.transAxes,
            ha="right",
            fontsize=8,
            color=C_TERRACOTTA,
            fontweight="bold",
        )

    fig.suptitle(
        "Supplemental Figure S4. Grapevine condensatome score distribution",
        fontsize=12,
        fontweight="bold",
        color=C_INK,
        y=1.01,
    )
    fig.tight_layout()
    return _savefig(fig, out_dir, "Figure_S4")


def _histogram(vals: Sequence[float], bins: Sequence[float]) -> Tuple[List[float], List[float]]:
    counts = [0.0] * (len(bins) - 1)
    for v in vals:
        if v <= bins[0]:
            counts[0] += 1
            continue
        if v >= bins[-1]:
            counts[-1] += 1
            continue
        # binary search
        lo, hi = 0, len(bins) - 1
        while lo < hi - 1:
            mid = (lo + hi) // 2
            if v >= bins[mid]:
                lo = mid
            else:
                hi = mid
        counts[lo] += 1
    return counts, list(bins)


# ---------------------------------------------------------------------------
# S5 — stress DE enrichment + modules
# ---------------------------------------------------------------------------


def figure_s5(out_dir: Path) -> Tuple[Path, Path]:
    enrich = list(
        csv.DictReader(
            (ROOT / "data" / "stress_modules" / "de_enrichment.tsv").open(),
            delimiter="\t",
        )
    )
    mods = list(
        csv.DictReader(
            (ROOT / "data" / "stress_modules" / "coexpression_modules.tsv").open(),
            delimiter="\t",
        )
    )
    summary = json.loads((ROOT / "data" / "stress_modules" / "summary.json").read_text())

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.4), facecolor=C_BG)

    # A — enrichment fold (all 18 tests)
    ax = axes[0]
    _style_ax(ax)
    # short labels
    labels = []
    folds = []
    fdrs = []
    colors = []
    for r in enrich:
        acc = r["accession"].replace("E-", "")
        lab = f"{acc} {r['contrast_id']} {r['direction']}"
        labels.append(lab)
        folds.append(float(r["enrichment_fold"]))
        fdr = float(r["bh_fdr"])
        fdrs.append(fdr)
        colors.append(C_TERRACOTTA if fdr < 0.05 else C_TEAL)

    ypos = list(range(len(labels)))[::-1]
    ax.barh(ypos, folds, color=colors, height=0.7, edgecolor="none")
    ax.axvline(1.0, color="#888888", lw=1, ls="--")
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.set_xlabel("Enrichment fold (Tier A∪B among DE)")
    ax.set_title(
        f"A  DE × condensatome enrichment (18 tests)\n"
        f"FDR<0.05: {summary['n_enrichment_fdr_lt_0.05']} "
        f"(terracotta); background rate={summary['parameters']['background_tier_ab_rate']:.3f}",
        fontsize=10,
        loc="left",
        color=C_INK,
    )
    # mark significant
    for y, f, fdr in zip(ypos, folds, fdrs):
        if fdr < 0.05:
            ax.text(f + 0.01, y, f"FDR={fdr:.3f}", va="center", fontsize=7, color=C_TERRACOTTA)

    # B — top modules by size
    ax = axes[1]
    _style_ax(ax)
    top = sorted(mods, key=lambda r: -int(r["n_members"]))[:15]
    labels = [
        f"{r['module_id'].replace('CX_stress_mod_', 'M')}  {r['pattern_label'][:42]}"
        for r in top
    ]
    n_a = [int(r["n_tier_a"]) for r in top]
    n_b = [int(r["n_tier_b"]) for r in top]
    ypos = list(range(len(top)))[::-1]
    ax.barh(ypos, n_a, color=C_TIER["A"], height=0.7, label="Tier A", edgecolor="none")
    ax.barh(ypos, n_b, left=n_a, color=C_TIER["B"], height=0.7, label="Tier B", edgecolor="none")
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Module members (Tier A∪B)")
    ax.set_title(
        f"B  Top response-pattern modules (of {summary['n_modules']})\n"
        f"E-MTAB-4324 / E-GEOD-67932 · not WGCNA",
        fontsize=10,
        loc="left",
        color=C_INK,
    )
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    fig.suptitle(
        "Supplemental Figure S5. Stress DE enrichment and condensatome modules",
        fontsize=12,
        fontweight="bold",
        color=C_INK,
        y=1.01,
    )
    fig.tight_layout()
    return _savefig(fig, out_dir, "Figure_S5")


# ---------------------------------------------------------------------------
# Legends
# ---------------------------------------------------------------------------

LEGENDS_MD = """# Supplemental Figure Legends (S1–S5)

Generated by `scripts/generate_supp_figures.py` from deposited v1.4.0 metrics.
Files: `Figure_S{{1..5}}.png` / `.pdf` in this directory.

## Supplemental Figure S1

**PlantLLPS v13 stratified cross-validation and gold-standard ROC.**
**(A)** Out-of-fold ROC curves from stratified five-fold CV on the formal v13
training corpus (*n* = 162; seed = 7), for logistic, hybrid and companion scores.
Mean / OOF logistic AUCs match `data/model_formal.metrics.v13.json`.
**(B)** Production-model ROC on the reconstructed gold panel (*n* = 28).
**Caveat (marked on panel):** gold AUCs are partially in-sample because many gold
IDs overlap the training set — do **not** cite as independent transfer. Prefer
leave-gold-out AUCs in Figure S2 / `docs/BENCHMARK_LEADERBOARD.md`.

## Supplemental Figure S2

**Fair leave-gold-out / family hold-out / Ath→non_ath_crop transfer versus
ProX-style, PSAP, LLPhyScore and PLAAC-approx.**
**(A)** Leave-gold-out ROC (*n*_test = 28; train *n* = 134; supervised methods
retrained excluding every test ID; PSAP/LLPhy/PLAAC out-of-the-box).
**(B)** Leave-gold-out AUC bar summary.
**(C)** Family hold-out mean AUC across 13 evaluated families.
**(D)** Ath→non_ath_crop species-transfer AUC (*n* = 38).
Metrics: `data/benchmarks/fair_comparison_metrics.json`; leaderboard
`docs/BENCHMARK_LEADERBOARD.md`. ProX-style is a feature-family logistic baseline,
not the published Plant PS proteins LightGBM/PU model.

## Supplemental Figure S3

**Feature ablation and F6HET4 Ser-rich / charge-block / IDR SHAP grammar.**
**(A)** Leave-one-group-out Δ mean logistic AUC (ablated − full) on v13 stratified
5-fold CV (`data/benchmarks/feature_ablation_metrics.json`). Largest drop:
`map65_cterm`.
**(B)** Exact logistic SHAP grammar axes for VvMAP65-1 (F6HET4, Tier B): Ser-rich/MAP65,
charge-block and IDR/disorder (`data/explainability/F6HET4.json`). Values are
probability-local linearizations of φᵢ = wᵢ·zᵢ — **not** KernelSHAP/TreeSHAP.

## Supplemental Figure S4

**Grapevine proteome-wide PlantLLPS score distribution (*n* = 25 972).**
**(A)** Stacked histogram of hybrid scores by confidence tier (A–D) with tier score
cutoffs and curated cross-crop MAP anchors (F6HET4 / MAP65-1; WVD2/WDL hooks when
present). Source: `data/grape_condensatome.sqlite` + `data/proteome_summary.json`.
**(B)** Tier counts (A = 1 142, B = 2 728, C = 9 653, D = 12 449).

## Supplemental Figure S5

**Stress DE × condensatome enrichment and response-pattern modules.**
**(A)** Hypergeometric enrichment of Tier A∪B among significant DE genes across 18
contrast×direction tests from Expression Atlas accessions **E-MTAB-4324** and
**E-GEOD-67932**. Terracotta bars: BH FDR < 0.05 (one test: E-GEOD-67932 stage-3
Botrytis down, fold ≈ 1.17).
**(B)** Largest Tier A∪B response-pattern modules (size ≥ 5; top 15 shown).
Descriptive overlap on short Atlas DE vectors — **not** WGCNA.
Source: `data/stress_modules/`.

## Reproducibility

```bash
python3 scripts/generate_supp_figures.py --out-dir docs/figures/supp
```

Pinned plotting stack (optional): see `requirements-figures.txt`
(`matplotlib==3.10.8`, `numpy==2.2.6`). Core PlantLLPS scoring remains stdlib-only.
"""


def write_legends(out_dir: Path) -> Path:
    path = out_dir / "LEGENDS.md"
    path.write_text(LEGENDS_MD)
    return path


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "docs" / "figures" / "supp",
        help="Output directory for PNG/PDF/LEGENDS.md",
    )
    p.add_argument(
        "--only",
        nargs="*",
        choices=["S1", "S2", "S3", "S4", "S5"],
        help="Optional subset of figures to regenerate",
    )
    args = p.parse_args(argv)
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    wanted = set(args.only) if args.only else {"S1", "S2", "S3", "S4", "S5"}
    generators = {
        "S1": figure_s1,
        "S2": figure_s2,
        "S3": figure_s3,
        "S4": figure_s4,
        "S5": figure_s5,
    }
    written: List[str] = []
    for key in ("S1", "S2", "S3", "S4", "S5"):
        if key not in wanted:
            continue
        png, pdf = generators[key](out_dir)
        written.append(f"{png.name} / {pdf.name}")
        print(f"wrote {png} and {pdf}")

    legends = write_legends(out_dir)
    print(f"wrote {legends}")
    print("done:", ", ".join(written))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
