"""Sequence features for plant LLPS propensity (composition heuristics).

When metapredict / AlphaFold tracks are unavailable, disorder is proxied by
a sliding-window composition heuristic (high charged/polar/Pro, low
hydrophobic burial). Low-complexity, charge patterning, aromatic and
Ser/Arg grammar follow the companion MAP/LLPS feature logic at reduced
complexity suitable for a stdlib prototype.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Dict, Iterable, List, Mapping, Sequence

AA = "ACDEFGHIKLMNPQRSTVWY"
HYDROPHOBIC = set("AILMFVW")
CHARGED = set("DEKR")
POSITIVE = set("KR")
NEGATIVE = set("DE")
POLAR = set("STNQYC")
AROMATIC = set("FYW")
DISORDER_PROMOTING = set("DESTPQKRNGHA")  # Uversky-ish promoting set (no Cys/Trp)

FEATURE_NAMES: List[str] = [
    "length",
    "frac_disorder_proxy",
    "min_entropy_12",
    "low_complexity_score",
    "ncpr",
    "kappa_charge_approx",
    "blockiness_ncpr_sd11",
    "frac_aromatic",
    "frac_ser",
    "frac_arg",
    "frac_ser_arg",
    "sp_tp_density",
    "qn_frac",
    "charge_frac",
    "pi_approx",
    "c_term_ser_frac_100",
    "c_term_ncpr_100",
    "map65_tail_like",
]


def _clean_seq(seq: str) -> str:
    s = "".join(ch for ch in seq.upper() if ch.isalpha())
    # map nonstandard to X (ignored in most frac counts)
    return "".join(ch if ch in AA else "X" for ch in s)


def _shannon_entropy(window: str) -> float:
    if not window:
        return 0.0
    n = len(window)
    counts = Counter(window)
    return -sum((c / n) * math.log2(c / n) for c in counts.values() if c)


def _disorder_proxy_fraction(seq: str, win: int = 21) -> float:
    """Fraction of residues whose local window is disorder-promoting.

    Heuristic only — not a substitute for metapredict or pLDDT.
    """
    if not seq:
        return 0.0
    half = win // 2
    n = len(seq)
    flags = []
    for i in range(n):
        a = max(0, i - half)
        b = min(n, i + half + 1)
        w = seq[a:b]
        promoting = sum(1 for r in w if r in DISORDER_PROMOTING) / len(w)
        buried = sum(1 for r in w if r in HYDROPHOBIC) / len(w)
        flags.append(1.0 if promoting >= 0.55 and buried <= 0.35 else 0.0)
    return sum(flags) / n


def _ncpr(seq: str) -> float:
    if not seq:
        return 0.0
    pos = sum(1 for r in seq if r in POSITIVE)
    neg = sum(1 for r in seq if r in NEGATIVE)
    return (pos - neg) / len(seq)


def _blockiness_ncpr_sd(seq: str, win: int = 11) -> float:
    if len(seq) < win:
        return 0.0
    vals = [_ncpr(seq[i : i + win]) for i in range(0, len(seq) - win + 1)]
    mean = sum(vals) / len(vals)
    var = sum((v - mean) ** 2 for v in vals) / len(vals)
    return math.sqrt(var)


def _kappa_approx(seq: str) -> float:
    """Rough charge-patterning proxy inspired by Das–Pappu κ.

    Compares observed vs mixed charge asymmetry using window NCPR variance.
    Range roughly [0, 1]; higher ≈ more segregated charges.
    """
    if len(seq) < 30:
        return 0.0
    obs = _blockiness_ncpr_sd(seq, 5)
    # composition-matched null: sort charges to alternate as much as possible
    charged = [r for r in seq if r in CHARGED]
    others = [r for r in seq if r not in CHARGED]
    if len(charged) < 4:
        return 0.0
    # maximally mixed arrangement of +/−
    plus = [r for r in charged if r in POSITIVE]
    minus = [r for r in charged if r in NEGATIVE]
    mixed_charged: List[str] = []
    while plus or minus:
        if plus:
            mixed_charged.append(plus.pop())
        if minus:
            mixed_charged.append(minus.pop())
    # Place charged residues at even spacing to approximate a mixed null.
    step = max(1, len(seq) // max(1, len(mixed_charged)))
    arr2 = ["G"] * len(seq)
    o_i = 0
    c_i = 0
    for i in range(len(seq)):
        if i % step == 0 and c_i < len(mixed_charged):
            arr2[i] = mixed_charged[c_i]
            c_i += 1
        else:
            if o_i < len(others):
                arr2[i] = others[o_i]
                o_i += 1
            elif c_i < len(mixed_charged):
                arr2[i] = mixed_charged[c_i]
                c_i += 1
    mixed_sd = _blockiness_ncpr_sd("".join(arr2), 5)
    # segregated extreme: all + then all -
    seg = "".join(sorted(seq, key=lambda r: (0 if r in POSITIVE else 1 if r in NEGATIVE else 2)))
    seg_sd = _blockiness_ncpr_sd(seg, 5)
    denom = seg_sd - mixed_sd
    if denom <= 1e-9:
        return 0.0
    return max(0.0, min(1.0, (obs - mixed_sd) / denom))


def _approx_pi(seq: str) -> float:
    """Very rough isoelectric-point proxy from net charge composition."""
    # Not a Henderson–Hasselbalch solver; linear map of NCPR for ranking.
    return 7.0 + 12.0 * _ncpr(seq)


def _sp_tp_density(seq: str) -> float:
    if len(seq) < 2:
        return 0.0
    count = 0
    for i in range(len(seq) - 1):
        if seq[i] in "ST" and seq[i + 1] == "P":
            count += 1
    return count / (len(seq) - 1)


def extract_features(sequence: str) -> Dict[str, float]:
    """Return a dict of plant-LLPS-relevant features for one protein sequence."""
    seq = _clean_seq(sequence)
    n = len(seq)
    if n == 0:
        return {name: 0.0 for name in FEATURE_NAMES}

    counts = Counter(seq)
    frac = {aa: counts.get(aa, 0) / n for aa in AA}

    # sliding entropy
    win = 12
    if n >= win:
        ents = [_shannon_entropy(seq[i : i + win]) for i in range(n - win + 1)]
        min_ent = min(ents)
    else:
        min_ent = _shannon_entropy(seq)

    # low-complexity score: invert min entropy (high score = more LC)
    # max entropy for 12-mer of 20 aa ≈ log2(min(12,20))
    max_ent = math.log2(min(win, 20))
    lc_score = 1.0 - (min_ent / max_ent if max_ent > 0 else 0.0)

    frac_ser = frac["S"]
    frac_arg = frac["R"]
    frac_aro = frac["F"] + frac["Y"] + frac["W"]
    frac_sr = frac_ser + frac_arg
    qn = frac["Q"] + frac["N"]
    charge_frac = sum(frac[a] for a in "DEKR")

    cterm = seq[-100:] if n >= 50 else seq
    c_ser = sum(1 for r in cterm if r == "S") / len(cterm)
    c_ncpr = _ncpr(cterm)

    # MAP65-1-like C-terminal grammar: Ser-rich, weakly cationic, aromatic-poor
    map65_like = 0.0
    if (
        c_ser >= 0.08
        and c_ncpr > 0.0
        and c_ncpr < 0.15
        and (sum(1 for r in cterm if r in AROMATIC) / len(cterm)) < 0.05
        and n >= 200
    ):
        map65_like = 1.0
    elif c_ser >= 0.10 and c_ncpr > 0.0 and n >= 150:
        map65_like = 0.5

    feats = {
        "length": float(n),
        "frac_disorder_proxy": _disorder_proxy_fraction(seq),
        "min_entropy_12": min_ent,
        "low_complexity_score": lc_score,
        "ncpr": _ncpr(seq),
        "kappa_charge_approx": _kappa_approx(seq),
        "blockiness_ncpr_sd11": _blockiness_ncpr_sd(seq, 11),
        "frac_aromatic": frac_aro,
        "frac_ser": frac_ser,
        "frac_arg": frac_arg,
        "frac_ser_arg": frac_sr,
        "sp_tp_density": _sp_tp_density(seq),
        "qn_frac": qn,
        "charge_frac": charge_frac,
        "pi_approx": _approx_pi(seq),
        "c_term_ser_frac_100": c_ser,
        "c_term_ncpr_100": c_ncpr,
        "map65_tail_like": map65_like,
    }
    return feats


def features_to_vector(feats: Mapping[str, float], names: Sequence[str] | None = None) -> List[float]:
    names = list(names) if names is not None else FEATURE_NAMES
    return [float(feats.get(n, 0.0)) for n in names]


def batch_extract(records: Iterable[tuple[str, str]]) -> List[Dict[str, float]]:
    out = []
    for pid, seq in records:
        f = extract_features(seq)
        f["protein_id"] = pid  # type: ignore[assignment]
        out.append(f)
    return out
