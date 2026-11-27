#!/usr/bin/env python3
"""Systematic feature-ablation panels for PlantLLPS (Figure 3).

Leave-one-group-out and group-only stratified 5-fold CV on train_formal_v13.csv,
matching the production CV seed/path. Family hold-out remains separate
(``scripts/eval_family_holdout.py``).

Deposits::

    data/benchmarks/feature_ablation_metrics.json
    data/benchmarks/feature_ablation_panel.svg
    data/benchmarks/FEATURE_ABLATION_README.md
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plantllps.features import FEATURE_NAMES, extract_features, features_to_vector  # noqa: E402
from plantllps.metrics import roc_auc, stratified_kfold_indices  # noqa: E402
from plantllps.model import PlantLLPSModel  # noqa: E402

DEFAULT_TRAIN = ROOT / "data" / "train_formal_v13.csv"
DEFAULT_OUT = ROOT / "data" / "benchmarks" / "feature_ablation_metrics.json"
DEFAULT_SVG = ROOT / "data" / "benchmarks" / "feature_ablation_panel.svg"
DEFAULT_README = ROOT / "data" / "benchmarks" / "FEATURE_ABLATION_README.md"

# Manuscript-aligned feature groups (Figure 2 / Methods axes).
FEATURE_GROUPS: Dict[str, List[str]] = {
    "disorder": ["frac_disorder_proxy"],
    "low_complexity": ["min_entropy_12", "low_complexity_score"],
    "charge_blockiness": [
        "ncpr",
        "kappa_charge_approx",
        "blockiness_ncpr_sd11",
        "charge_frac",
        "pi_approx",
    ],
    "aromatic": ["frac_aromatic"],
    "ser_arg": ["frac_ser", "frac_arg", "frac_ser_arg"],
    "sp_tp": ["sp_tp_density"],
    "map65_cterm": ["c_term_ser_frac_100", "c_term_ncpr_100", "map65_tail_like"],
    "qn_polar": ["qn_frac"],
    "length": ["length"],
}


def _validate_groups() -> None:
    covered = {n for names in FEATURE_GROUPS.values() for n in names}
    missing = [n for n in FEATURE_NAMES if n not in covered]
    extra = sorted(covered - set(FEATURE_NAMES))
    if missing or extra:
        raise RuntimeError(f"FEATURE_GROUPS mismatch: missing={missing} extra={extra}")


def load_train(path: Path):
    ids, seqs, labels = [], [], []
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            ids.append(row["protein_id"])
            seqs.append(row["sequence"].strip())
            lab = row["label"].strip().lower()
            labels.append(1 if lab.startswith("pos") else 0)
    return ids, seqs, labels


def build_X(seqs: Sequence[str]) -> List[List[float]]:
    return [features_to_vector(extract_features(s), FEATURE_NAMES) for s in seqs]


def mask_features(
    X: Sequence[Sequence[float]],
    keep: Sequence[str] | None = None,
    drop: Sequence[str] | None = None,
) -> Tuple[List[List[float]], List[str]]:
    """Zero-mask features not in ``keep`` (or in ``drop``). Names stay fixed length."""
    names = list(FEATURE_NAMES)
    if keep is not None:
        keep_set = set(keep)
        active = [n for n in names if n in keep_set]
    elif drop is not None:
        drop_set = set(drop)
        active = [n for n in names if n not in drop_set]
    else:
        active = names
    active_set = set(active)
    out: List[List[float]] = []
    for row in X:
        feats = {n: v for n, v in zip(names, row)}
        vec = [float(feats[n]) if n in active_set else 0.0 for n in names]
        out.append(vec)
    return out, active


def _mean(xs: Sequence[float]) -> float:
    ok = [x for x in xs if x == x]
    return sum(ok) / len(ok) if ok else float("nan")


def cross_validate_masked(
    X: Sequence[Sequence[float]],
    y: Sequence[int],
    *,
    active_names: Sequence[str],
    k: int = 5,
    seed: int = 7,
) -> Dict[str, Any]:
    """Stratified k-fold AUC for logistic + companion on a masked feature matrix."""
    folds = stratified_kfold_indices(y, k=k, seed=seed)
    fold_aucs_log: List[float] = []
    fold_aucs_comp: List[float] = []
    fold_aucs_hyb: List[float] = []
    oof_log = [float("nan")] * len(y)
    oof_comp = [float("nan")] * len(y)
    oof_hyb = [float("nan")] * len(y)
    for train_idx, test_idx in folds:
        Xtr = [X[i] for i in train_idx]
        ytr = [y[i] for i in train_idx]
        m = PlantLLPSModel()
        m.fit_logistic(Xtr, ytr, seed=seed)
        y_te = [y[i] for i in test_idx]
        s_log, s_comp, s_hyb = [], [], []
        for i in test_idx:
            feats = {n: v for n, v in zip(FEATURE_NAMES, X[i])}
            # Companion uses the (possibly masked) disorder/LC/blockiness values.
            sl = m.logistic_score(feats)
            sc = m.companion_score(feats)
            sh = 0.5 * sl + 0.5 * sc
            s_log.append(sl)
            s_comp.append(sc)
            s_hyb.append(sh)
            oof_log[i] = sl
            oof_comp[i] = sc
            oof_hyb[i] = sh
        fold_aucs_log.append(roc_auc(y_te, s_log))
        fold_aucs_comp.append(roc_auc(y_te, s_comp))
        fold_aucs_hyb.append(roc_auc(y_te, s_hyb))

    def oof_auc(scores: Sequence[float]) -> float:
        pairs = [(scores[i], y[i]) for i in range(len(y)) if scores[i] == scores[i]]
        return roc_auc([p[1] for p in pairs], [p[0] for p in pairs])

    return {
        "n_active_features": len(active_names),
        "active_features": list(active_names),
        "k": k,
        "fold_auc_logistic": fold_aucs_log,
        "fold_auc_companion": fold_aucs_comp,
        "fold_auc_hybrid": fold_aucs_hyb,
        "mean_auc_logistic": _mean(fold_aucs_log),
        "mean_auc_companion": _mean(fold_aucs_comp),
        "mean_auc_hybrid": _mean(fold_aucs_hyb),
        "oof_auc_logistic": oof_auc(oof_log),
        "oof_auc_companion": oof_auc(oof_comp),
        "oof_auc_hybrid": oof_auc(oof_hyb),
    }


def single_feature_aucs(X: Sequence[Sequence[float]], y: Sequence[int]) -> Dict[str, float]:
    """Univariate ROC AUC for each feature (direction-aware max of raw / inverted)."""
    out: Dict[str, float] = {}
    for j, name in enumerate(FEATURE_NAMES):
        scores = [row[j] for row in X]
        a = roc_auc(y, scores)
        b = roc_auc(y, [-s for s in scores])
        out[name] = max(a, b) if a == a and b == b else (a if a == a else b)
    return out


def svg_bar_panel(
    logo_rows: Sequence[Mapping[str, Any]],
    only_rows: Sequence[Mapping[str, Any]],
    *,
    full_auc: float,
    title: str = "PlantLLPS feature ablation (v13 stratified 5-fold)",
) -> str:
    """Minimal SVG dual-panel bar chart (stdlib; no matplotlib)."""
    W, H = 920, 520
    pad_l, pad_r, pad_t, pad_b = 210, 40, 56, 48
    mid = W // 2
    plot_h = H - pad_t - pad_b

    def bars(rows, x0, x1, ylabel, key="delta_mean_auc_logistic"):
        n = max(1, len(rows))
        gap = 6
        bh = max(10, (plot_h - gap * (n + 1)) / n)
        # scale
        vals = [float(r.get(key, 0.0) or 0.0) for r in rows]
        vmax = max(0.05, max(abs(v) for v in vals) if vals else 0.05)
        cx = x0 + (x1 - x0) * 0.55  # zero line offset for signed deltas
        parts = [
            f'<text x="{(x0+x1)/2:.1f}" y="{pad_t-18}" text-anchor="middle" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="13" font-weight="700">{ylabel}</text>',
            f'<line x1="{cx:.1f}" y1="{pad_t}" x2="{cx:.1f}" y2="{H-pad_b}" '
            f'stroke="#888" stroke-width="1"/>',
        ]
        for i, r in enumerate(rows):
            y = pad_t + gap + i * (bh + gap)
            v = float(r.get(key, 0.0) or 0.0)
            w = (abs(v) / vmax) * ((x1 - x0) * 0.40)
            if key.startswith("delta"):
                # negative delta (drop hurts) → bar left of zero in coral
                if v <= 0:
                    x = cx - w
                    color = "#c45c26"
                else:
                    x = cx
                    color = "#2a6f6f"
                label = f"{v:+.3f}"
            else:
                x = x0 + 8
                w = (max(0.0, v) / max(vmax, 1e-9)) * ((x1 - x0) - 16)
                color = "#2a6f6f"
                label = f"{v:.3f}"
            name = r.get("group") or r.get("label") or ""
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w,0.5):.1f}" height="{bh:.1f}" '
                f'fill="{color}" opacity="0.85"/>'
            )
            parts.append(
                f'<text x="{x0+4:.1f}" y="{y+bh*0.72:.1f}" font-family="Helvetica,Arial,sans-serif" '
                f'font-size="11" fill="#222">{name}</text>'
            )
            parts.append(
                f'<text x="{x1-4:.1f}" y="{y+bh*0.72:.1f}" text-anchor="end" '
                f'font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#333">{label}</text>'
            )
        return "\n".join(parts)

    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        f'<rect width="{W}" height="{H}" fill="#fbfaf7"/>',
        f'<text x="{W/2:.1f}" y="28" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" '
        f'font-size="16" font-weight="700" fill="#1a1a1a">{title}</text>',
        f'<text x="{W/2:.1f}" y="46" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" '
        f'font-size="11" fill="#555">full-model mean logistic AUC = {full_auc:.3f} · '
        f"Δ = ablated − full (negative = group contributes)</text>",
        bars(logo_rows, 10, mid - 10, "Leave-one-group-out ΔAUC (logistic)"),
        bars(
            only_rows,
            mid + 10,
            W - 10,
            "Group-only mean AUC (logistic)",
            key="mean_auc_logistic",
        ),
        "</svg>",
    ]
    return "\n".join(body)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train-csv", type=Path, default=DEFAULT_TRAIN)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--svg", type=Path, default=DEFAULT_SVG)
    ap.add_argument("--readme", type=Path, default=DEFAULT_README)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--folds", type=int, default=5)
    args = ap.parse_args()

    _validate_groups()
    ids, seqs, labels = load_train(args.train_csv)
    X = build_X(seqs)
    y = labels
    k = min(args.folds, sum(y), len(y) - sum(y))

    full_X, full_names = mask_features(X)
    full = cross_validate_masked(full_X, y, active_names=full_names, k=k, seed=args.seed)

    leave_one: List[Dict[str, Any]] = []
    for gname, gfeats in FEATURE_GROUPS.items():
        Xm, active = mask_features(X, drop=gfeats)
        metrics = cross_validate_masked(Xm, y, active_names=active, k=k, seed=args.seed)
        leave_one.append(
            {
                "group": gname,
                "dropped_features": list(gfeats),
                "mean_auc_logistic": metrics["mean_auc_logistic"],
                "mean_auc_companion": metrics["mean_auc_companion"],
                "mean_auc_hybrid": metrics["mean_auc_hybrid"],
                "oof_auc_logistic": metrics["oof_auc_logistic"],
                "delta_mean_auc_logistic": metrics["mean_auc_logistic"] - full["mean_auc_logistic"],
                "delta_mean_auc_hybrid": metrics["mean_auc_hybrid"] - full["mean_auc_hybrid"],
                "cv": metrics,
            }
        )

    group_only: List[Dict[str, Any]] = []
    for gname, gfeats in FEATURE_GROUPS.items():
        Xm, active = mask_features(X, keep=gfeats)
        metrics = cross_validate_masked(Xm, y, active_names=active, k=k, seed=args.seed)
        group_only.append(
            {
                "group": gname,
                "kept_features": list(gfeats),
                "mean_auc_logistic": metrics["mean_auc_logistic"],
                "mean_auc_companion": metrics["mean_auc_companion"],
                "mean_auc_hybrid": metrics["mean_auc_hybrid"],
                "oof_auc_logistic": metrics["oof_auc_logistic"],
                "cv": metrics,
            }
        )

    # Sort panels for readability: largest contribution (most negative Δ) first
    leave_one_sorted = sorted(leave_one, key=lambda r: r["delta_mean_auc_logistic"])
    group_only_sorted = sorted(group_only, key=lambda r: -r["mean_auc_logistic"])

    univ = single_feature_aucs(X, y)

    report = {
        "version": 1,
        "figure": "Figure 3 feature-ablation panels",
        "train_csv": str(args.train_csv.relative_to(ROOT)),
        "n": len(y),
        "n_pos": int(sum(y)),
        "n_neg": int(len(y) - sum(y)),
        "seed": args.seed,
        "k": k,
        "feature_groups": FEATURE_GROUPS,
        "full_model": full,
        "leave_one_group_out": leave_one_sorted,
        "group_only": group_only_sorted,
        "univariate_feature_auc": univ,
        "honesty": (
            "Ablation uses the same stratified 5-fold seed as production CV on the "
            "composition-proxy feature matrix. Zero-masking keeps the logistic input "
            "dimension fixed so fold splits remain comparable. Companion composite only "
            "depends on disorder / low-complexity / blockiness (+ MAP65 bonus terms); "
            "groups outside that set mainly affect logistic AUC. Family hold-out AUCs "
            "remain in family_holdout_metrics.json and are not replaced by these panels."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2))
    args.svg.write_text(
        svg_bar_panel(
            leave_one_sorted,
            group_only_sorted,
            full_auc=full["mean_auc_logistic"],
        )
    )

    # Compact TSV for plotting elsewhere
    tsv_path = args.out.with_suffix(".tsv")
    with tsv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "panel",
                "group",
                "mean_auc_logistic",
                "mean_auc_hybrid",
                "delta_mean_auc_logistic",
                "delta_mean_auc_hybrid",
            ]
        )
        for r in leave_one_sorted:
            w.writerow(
                [
                    "leave_one_group_out",
                    r["group"],
                    f"{r['mean_auc_logistic']:.6f}",
                    f"{r['mean_auc_hybrid']:.6f}",
                    f"{r['delta_mean_auc_logistic']:.6f}",
                    f"{r['delta_mean_auc_hybrid']:.6f}",
                ]
            )
        for r in group_only_sorted:
            w.writerow(
                [
                    "group_only",
                    r["group"],
                    f"{r['mean_auc_logistic']:.6f}",
                    f"{r['mean_auc_hybrid']:.6f}",
                    "",
                    "",
                ]
            )

    top_drop = leave_one_sorted[0] if leave_one_sorted else None
    top_only = group_only_sorted[0] if group_only_sorted else None
    drop_line = (
        f"- Largest leave-one-group-out drop: **{top_drop['group']}** "
        f"(Δlogistic ≈ **{top_drop['delta_mean_auc_logistic']:.3f}**)"
        if top_drop
        else "- Largest leave-one-group-out drop: n/a"
    )
    only_line = (
        f"- Strongest group-only logistic AUC: **{top_only['group']}** "
        f"(≈ **{top_only['mean_auc_logistic']:.3f}**)"
        if top_only
        else "- Strongest group-only logistic AUC: n/a"
    )
    args.readme.write_text(
        "\n".join(
            [
                "# Feature-ablation panels (Figure 3)",
                "",
                "Generated by `scripts/eval_feature_ablation.py`.",
                "",
                f"## Headline (v13 train, seed={args.seed}, {k}-fold)",
                "",
                f"- Full-model mean logistic AUC ≈ **{full['mean_auc_logistic']:.3f}** "
                f"(hybrid ≈ **{full['mean_auc_hybrid']:.3f}**)",
                drop_line,
                only_line,
                "",
                "## Files",
                "",
                "- `feature_ablation_metrics.json` — full CV dumps",
                "- `feature_ablation_metrics.tsv` — compact panel table",
                "- `feature_ablation_panel.svg` — dual-panel bar chart",
                "",
                "## Honesty",
                "",
                "Composition-proxy disorder; same CV protocol as production metrics JSON. Family",
                "hold-out remains the transfer probe (`family_holdout_metrics.json`).",
                "",
            ]
        )
    )

    print(
        f"Full logistic AUC={full['mean_auc_logistic']:.3f} hybrid={full['mean_auc_hybrid']:.3f} "
        f"→ {args.out}",
        file=sys.stderr,
    )
    if top_drop:
        print(
            f"Largest LOGO drop: {top_drop['group']} Δ={top_drop['delta_mean_auc_logistic']:+.3f}",
            file=sys.stderr,
        )
    print(f"SVG → {args.svg}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
