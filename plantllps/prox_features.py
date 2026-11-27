"""ProX / Plant-PS–style sequence features (stdlib approximation).

Plant PS proteins / ProX (Zhang et al. 2023 *Mol. Plant*; PU predictor follow-on)
score plant condensates from PLAAC / catGRANULE / IDR / LCR / hydropathy / FCR /
PScore-class inputs under a PU / LightGBM regime. This module extracts a
transparent, stdlib feature vector in that family so a leakage-controlled
logistic baseline can sit on the same train/test splits as PlantLLPS.

Honesty: this is **not** a reimplementation of the published LightGBM/PU model,
catGRANULE, or Java PLAAC. Scores are comparable only as a feature-family
baseline under identical splits and leakage rules.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Dict, List, Mapping, Sequence

from .features import (
    CHARGED,
    _blockiness_ncpr_sd,
    _clean_seq,
    _disorder_proxy_fraction,
    _kappa_approx,
    _ncpr,
    _shannon_entropy,
)
from .plaac_approx import plaac_features

# Kyte–Doolittle hydrophobicity (standard table)
_KD = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

PROX_FEATURE_NAMES: List[str] = [
    "log_length",
    "idr_frac",
    "lcr_score",
    "hydropathy_kd",
    "fcr",
    "ncpr",
    "kappa_charge_approx",
    "aromatic_frac",
    "plaac_qn_frac",
    "plaac_max_window_qn",
    "plaac_nllr_approx",
    "plaac_core_len_norm",
]


def _hydropathy_kd(seq: str) -> float:
    vals = [_KD[r] for r in seq if r in _KD]
    if not vals:
        return 0.0
    return sum(vals) / len(vals)


def _lcr_score(seq: str, win: int = 12) -> float:
    n = len(seq)
    if n == 0:
        return 0.0
    if n >= win:
        ents = [_shannon_entropy(seq[i : i + win]) for i in range(n - win + 1)]
        min_ent = min(ents)
    else:
        min_ent = _shannon_entropy(seq)
    max_ent = math.log2(min(win, 20))
    return 1.0 - (min_ent / max_ent if max_ent > 0 else 0.0)


def extract_prox_features(sequence: str) -> Dict[str, float]:
    """Return ProX/PU-style feature dict for one sequence."""
    seq = _clean_seq(sequence)
    n = len(seq)
    empty = {name: 0.0 for name in PROX_FEATURE_NAMES}
    if n == 0:
        return empty

    counts = Counter(seq)
    aromatic = (counts.get("F", 0) + counts.get("Y", 0) + counts.get("W", 0)) / n
    fcr = sum(1 for r in seq if r in CHARGED) / n
    plaac = plaac_features(seq)
    core_norm = plaac["plaac_core_len"] / n

    return {
        "log_length": math.log(max(1, n)),
        "idr_frac": _disorder_proxy_fraction(seq),
        "lcr_score": _lcr_score(seq),
        "hydropathy_kd": _hydropathy_kd(seq),
        "fcr": fcr,
        "ncpr": _ncpr(seq),
        "kappa_charge_approx": _kappa_approx(seq),
        "aromatic_frac": aromatic,
        "plaac_qn_frac": plaac["plaac_qn_frac"],
        "plaac_max_window_qn": plaac["plaac_max_window_qn"],
        "plaac_nllr_approx": plaac["plaac_nllr_approx"],
        "plaac_core_len_norm": core_norm,
        # Extra diagnostics (not in PROX_FEATURE_NAMES training vector)
        "blockiness_ncpr_sd11": _blockiness_ncpr_sd(seq, 11),
        "plaac_core_len": plaac["plaac_core_len"],
        "length": float(n),
    }


def prox_features_to_vector(
    feats: Mapping[str, float], names: Sequence[str] | None = None
) -> List[float]:
    names = list(names) if names is not None else PROX_FEATURE_NAMES
    return [float(feats.get(n, 0.0)) for n in names]
