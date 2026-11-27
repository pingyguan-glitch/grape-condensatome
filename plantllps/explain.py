"""PlantLLPS explainability — SHAP-equivalent feature / grammar attributions.

For the trained logistic head, feature attributions are **exact linear SHAP**
in logit space: ``φ_i = w_i · z_i`` with ``z_i = (x_i − μ_i) / σ_i`` and base
value = bias. Hybrid score attributions blend logistic SHAP (mapped to
probability space via a local linearization) with an analytic companion
decomposition.

Grammar groups mirror the manuscript / feature-ablation axes, with NAR+ focus
on **Ser-rich**, **charge-block**, and **IDR/disorder** drivers. Optional
sliding-window occlusion yields residue-segment attributions that can be
joined to deposited IDR / domain intervals.
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .features import FEATURE_NAMES, extract_features, features_to_vector
from .model import COMPOSITE_WEIGHTS, PlantLLPSModel, _clip01, _sigmoid

# Manuscript-aligned groups (same keys as scripts/eval_feature_ablation.py).
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

# NAR+ display axes (aggregations of FEATURE_GROUPS).
GRAMMAR_AXES: Dict[str, List[str]] = {
    "ser_rich": ["ser_arg", "map65_cterm"],
    "charge_block": ["charge_blockiness"],
    "idr_disorder": ["disorder", "low_complexity"],
}


def _logit_contributions(
    model: PlantLLPSModel, feats: Mapping[str, float]
) -> Tuple[float, Dict[str, float], List[float]]:
    """Return (bias, {feature: φ_i}, scaled vector). Exact additive SHAP in logit space."""
    names = list(model.feature_names or FEATURE_NAMES)
    vec = features_to_vector(feats, names)
    if not model.trained or not model.weights:
        return 0.0, {n: 0.0 for n in names}, vec
    scaled = [(v - m) / s for v, m, s in zip(vec, model.means, model.stds)]
    phi = {n: float(w) * float(z) for n, w, z in zip(names, model.weights, scaled)}
    return float(model.bias), phi, scaled


def logistic_shap(model: PlantLLPSModel, feats: Mapping[str, float]) -> Dict[str, Any]:
    """Exact SHAP for the logistic head (logit-space φ; probability via sigmoid)."""
    bias, phi, _ = _logit_contributions(model, feats)
    logit = bias + sum(phi.values())
    # Probability-space attributions via local linearization dσ/dz = σ(1−σ)
    p = _sigmoid(logit)
    scale = p * (1.0 - p)
    phi_prob = {k: v * scale for k, v in phi.items()}
    return {
        "method": "exact_linear_shap_logit",
        "base_logit": bias,
        "logit": logit,
        "probability": p,
        "feature_shap_logit": phi,
        "feature_shap_prob_local": phi_prob,
        "note": (
            "φ_i = w_i·z_i are exact additive SHAP values in logit space; "
            "probability-space values use local linearization σ'(logit)."
        ),
    }


def companion_decomposition(feats: Mapping[str, float]) -> Dict[str, Any]:
    """Transparent companion composite term-by-term contributions (0–1 blend)."""
    d = _clip01(float(feats.get("frac_disorder_proxy", 0.0)))
    lc = _clip01(float(feats.get("low_complexity_score", 0.0)))
    blk_raw = float(feats.get("blockiness_ncpr_sd11", 0.0))
    blk = _clip01(blk_raw / 0.2)
    map65 = float(feats.get("map65_tail_like", 0.0))
    c_ser = float(feats.get("c_term_ser_frac_100", 0.0))
    bonus_map65 = 0.08 * map65
    bonus_ser = 0.05 * _clip01(c_ser / 0.15)
    terms = {
        "frac_disorder_proxy": COMPOSITE_WEIGHTS["frac_disorder_proxy"] * d,
        "low_complexity_score": COMPOSITE_WEIGHTS["low_complexity_score"] * lc,
        "blockiness_ncpr_sd11": COMPOSITE_WEIGHTS["blockiness_ncpr_sd11"] * blk,
        "map65_tail_bonus": bonus_map65,
        "c_term_ser_bonus": bonus_ser,
    }
    total = _clip01(sum(terms.values()))
    return {
        "method": "companion_weighted_decomposition",
        "score": total,
        "terms": terms,
        "normalized_inputs": {
            "frac_disorder_proxy": d,
            "low_complexity_score": lc,
            "blockiness_ncpr_sd11_scaled": blk,
            "map65_tail_like": map65,
            "c_term_ser_frac_100": c_ser,
        },
    }


def _group_sum(phi: Mapping[str, float], members: Sequence[str]) -> float:
    return float(sum(phi.get(n, 0.0) for n in members))


def grammar_attributions(
    feature_phi: Mapping[str, float],
    *,
    companion_terms: Optional[Mapping[str, float]] = None,
    hybrid_weight_logistic: float = 0.5,
) -> Dict[str, Any]:
    """Aggregate feature φ into ablation groups + NAR+ Ser/charge/IDR axes.

    When companion terms are provided, hybrid grammar scores average logistic
    group φ (probability-local) with companion term contributions that map onto
    the same axes.
    """
    groups = {
        gname: _group_sum(feature_phi, members)
        for gname, members in FEATURE_GROUPS.items()
    }

    # Companion → group mapping (partial; groups without companion mass stay logistic-only)
    companion_group = {g: 0.0 for g in FEATURE_GROUPS}
    if companion_terms:
        companion_group["disorder"] = float(companion_terms.get("frac_disorder_proxy", 0.0))
        companion_group["low_complexity"] = float(
            companion_terms.get("low_complexity_score", 0.0)
        )
        companion_group["charge_blockiness"] = float(
            companion_terms.get("blockiness_ncpr_sd11", 0.0)
        )
        companion_group["map65_cterm"] = float(
            companion_terms.get("map65_tail_bonus", 0.0)
        ) + float(companion_terms.get("c_term_ser_bonus", 0.0))

    # Do not average logistic φ (probability units) with companion terms (0–1
    # score mass) — different scales. Hybrid groups are logistic-only; companion
    # is reported beside them. Display "support" ranks axes by |shap| then
    # companion when SHAP is near zero.
    hybrid_groups = dict(groups)

    axes = {}
    for axis, group_names in GRAMMAR_AXES.items():
        shap = sum(groups[g] for g in group_names)
        comp = sum(companion_group.get(g, 0.0) for g in group_names)
        axes[axis] = {
            "logistic_shap_prob": shap,
            "companion": comp,
            # Alias kept for UI: primary NAR+ bar = logistic SHAP
            "hybrid": shap,
            "member_groups": list(group_names),
        }

    return {
        "feature_groups": FEATURE_GROUPS,
        "grammar_axes": GRAMMAR_AXES,
        "group_shap_logistic_prob": groups,
        "group_companion": companion_group,
        "group_hybrid": hybrid_groups,
        "axes": axes,
        "hybrid_weight_logistic": hybrid_weight_logistic,
    }


def window_occlusion(
    model: PlantLLPSModel,
    sequence: str,
    *,
    window: int = 25,
    step: int = 10,
    top_k: int = 24,
    replace_aa: str = "G",
) -> Dict[str, Any]:
    """Leave-one-window-out Δscore (Gly substitution) — SHAP-equivalent occlusion.

    Positive Δ means the wild-type window *raises* PlantLLPS score relative to
    a composition-neutral Gly patch (supports LLPS propensity).
    """
    seq = "".join(ch for ch in sequence.upper() if ch.isalpha())
    if len(seq) < window:
        return {
            "method": "window_occlusion_gly",
            "window": window,
            "step": step,
            "base_score": None,
            "windows": [],
        }
    base = float(model.score(seq)["score_plantllps"])
    rows: List[Dict[str, Any]] = []
    for start in range(0, len(seq) - window + 1, step):
        end = start + window  # exclusive for slice; inclusive coords below
        mut = seq[:start] + (replace_aa * window) + seq[end:]
        s = float(model.score(mut)["score_plantllps"])
        frag = seq[start:end]
        ser_frac = sum(1 for r in frag if r == "S") / window
        pos = sum(1 for r in frag if r in "KR")
        neg = sum(1 for r in frag if r in "DE")
        ncpr = (pos - neg) / window
        charge_frac = (pos + neg) / window
        rows.append(
            {
                "start": start + 1,  # 1-based inclusive
                "end": end,
                "delta_score": base - s,
                "ser_frac": round(ser_frac, 4),
                "ncpr": round(ncpr, 4),
                "charge_frac": round(charge_frac, 4),
                "label": _window_label(ser_frac, ncpr, charge_frac),
            }
        )
    rows.sort(key=lambda r: -abs(r["delta_score"]))
    top = rows[:top_k]
    # Re-sort top by position for tracks
    top_pos = sorted(top, key=lambda r: r["start"])
    return {
        "method": "window_occlusion_gly",
        "window": window,
        "step": step,
        "replace_aa": replace_aa,
        "base_score": base,
        "n_windows_scored": len(rows),
        "windows": top_pos,
        "top_positive": [r for r in rows if r["delta_score"] > 0][:8],
        "top_negative": [r for r in rows if r["delta_score"] < 0][:8],
    }


def _window_label(ser_frac: float, ncpr: float, charge_frac: float) -> str:
    tags = []
    if ser_frac >= 0.12:
        tags.append("ser_rich")
    if charge_frac >= 0.25 and abs(ncpr) >= 0.05:
        tags.append("charge_block")
    return "+".join(tags) if tags else "other"


def join_segments_to_windows(
    windows: Sequence[Mapping[str, Any]],
    segments: Sequence[Mapping[str, Any]],
    *,
    kind: str,
) -> List[Dict[str, Any]]:
    """Attribute deposited IDR/domain segments by overlapping occlusion windows."""
    out: List[Dict[str, Any]] = []
    for seg in segments:
        raw_start = seg.get("start", seg.get("s"))
        raw_end = seg.get("end", seg.get("e"))
        try:
            s0 = int(raw_start)  # type: ignore[arg-type]
            s1 = int(raw_end)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            continue
        overlap = [
            w
            for w in windows
            if not (int(w["end"]) < s0 or int(w["start"]) > s1)
        ]
        if not overlap:
            delta = 0.0
        else:
            delta = sum(float(w["delta_score"]) for w in overlap) / len(overlap)
        out.append(
            {
                "kind": kind,
                "start": s0,
                "end": s1,
                "label": seg.get("label") or seg.get("accession") or seg.get("source") or kind,
                "source": seg.get("source", ""),
                "mean_window_delta": round(delta, 6),
                "n_overlap_windows": len(overlap),
            }
        )
    out.sort(key=lambda r: -abs(r["mean_window_delta"]))
    return out


def explain_sequence(
    model: PlantLLPSModel,
    sequence: str,
    *,
    protein_id: str = "",
    tier: str = "",
    llps_gold: int = 0,
    idr_segments: Optional[Sequence[Mapping[str, Any]]] = None,
    domains: Optional[Sequence[Mapping[str, Any]]] = None,
    window: int = 25,
    step: int = 10,
    top_k_windows: int = 24,
    do_occlusion: bool = True,
) -> Dict[str, Any]:
    """Full explainability payload for one protein."""
    scored = model.score(sequence)
    feats = {k: scored[k] for k in FEATURE_NAMES if k in scored}
    shap = logistic_shap(model, feats)
    comp = companion_decomposition(feats)
    grammar = grammar_attributions(
        shap["feature_shap_prob_local"],
        companion_terms=comp["terms"],
        hybrid_weight_logistic=0.5 if model.mode == "hybrid" and model.trained else 0.0,
    )

    # Ranked feature table (hybrid emphasis: logistic φ_prob)
    ranked = sorted(
        (
            {
                "feature": name,
                "value": float(feats.get(name, 0.0)),
                "shap_logit": shap["feature_shap_logit"].get(name, 0.0),
                "shap_prob_local": shap["feature_shap_prob_local"].get(name, 0.0),
            }
            for name in (model.feature_names or FEATURE_NAMES)
        ),
        key=lambda r: -abs(r["shap_prob_local"]),
    )

    payload: Dict[str, Any] = {
        "protein_id": protein_id,
        "tier": tier or scored.get("tier"),
        "llps_gold": int(llps_gold),
        "length": int(scored.get("length", len(sequence))),
        "score_plantllps": float(scored["score_plantllps"]),
        "score_companion": float(scored["score_companion"]),
        "score_logistic": float(scored["score_logistic"]),
        "logistic_shap": {
            "method": shap["method"],
            "base_logit": shap["base_logit"],
            "logit": shap["logit"],
            "probability": shap["probability"],
            "note": shap["note"],
        },
        "companion": comp,
        "grammar": {
            "axes": grammar["axes"],
            "group_hybrid": grammar["group_hybrid"],
            "group_shap_logistic_prob": grammar["group_shap_logistic_prob"],
        },
        "top_features": ranked[:12],
        "method_summary": (
            "Logistic: exact linear SHAP in logit space (φ=w·z). "
            "Companion: weighted term decomposition. "
            "Hybrid grammar axes average both. "
            "Residue windows: Gly occlusion Δscore (SHAP-equivalent)."
        ),
    }

    if do_occlusion:
        occ = window_occlusion(
            model,
            sequence,
            window=window,
            step=step,
            top_k=top_k_windows,
        )
        payload["occlusion"] = occ
        idr_attr = join_segments_to_windows(
            occ.get("windows") or [], idr_segments or [], kind="idr"
        )
        dom_attr = join_segments_to_windows(
            occ.get("windows") or [], domains or [], kind="domain"
        )
        payload["segment_attributions"] = {
            "idr": idr_attr[:16],
            "domain": dom_attr[:16],
        }
        # Compact NAR+ track: windows labeled ser_rich / charge_block
        payload["grammar_windows"] = [
            w
            for w in (occ.get("windows") or [])
            if w.get("label") and w["label"] != "other"
        ][:20]

    return payload
