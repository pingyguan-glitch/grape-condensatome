"""Stdlib PLAAC-like prion / Q–N propensity scores (not the Java PLAAC HMM).

Companion analyses report PLAAC NLLR for MAP65 tails (typically near zero / negative).
This module approximates that signal with a transparent windowed Q/N enrichment score
plus a simple core-length statistic so every proteome row can carry a comparable column
without a Java runtime.

Output fields (joined into SQLite / features):
  * ``plaac_qn_frac`` — whole-sequence (Q+N) fraction
  * ``plaac_max_window_qn`` — max (Q+N) fraction in a 60-aa window
  * ``plaac_nllr_approx`` — log-likelihood-ratio style score vs background 0.1 Q/N
  * ``plaac_core_len`` — longest run of residues in windows with Q/N ≥ 0.30
"""

from __future__ import annotations

import math
from typing import Dict


def plaac_features(sequence: str, *, window: int = 60, background_qn: float = 0.10) -> Dict[str, float]:
    seq = "".join(ch for ch in sequence.upper() if ch.isalpha())
    n = len(seq)
    if n == 0:
        return {
            "plaac_qn_frac": 0.0,
            "plaac_max_window_qn": 0.0,
            "plaac_nllr_approx": 0.0,
            "plaac_core_len": 0.0,
        }
    qn_total = sum(1 for r in seq if r in "QN")
    qn_frac = qn_total / n
    win = min(window, n)
    max_w = 0.0
    hot = [0] * n
    for i in range(0, n - win + 1):
        w = seq[i : i + win]
        f = sum(1 for r in w if r in "QN") / win
        if f > max_w:
            max_w = f
        if f >= 0.30:
            for j in range(i, i + win):
                hot[j] = 1
    # longest contiguous "core"
    best = cur = 0
    for h in hot:
        if h:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    # NLLR-like: compare observed Q/N rate to background (Bernoulli LLR per residue, mean)
    p = min(0.999, max(1e-6, qn_frac))
    b = min(0.999, max(1e-6, background_qn))
    # per-residue mean LLR under observed composition vs background rate
    nllr = float(
        qn_frac * math.log(p / b) + (1.0 - qn_frac) * math.log((1.0 - p) / (1.0 - b))
    )
    # scale to a familiar magnitude; positive ⇒ Q/N-enriched vs background
    nllr_approx = nllr * n / 100.0
    return {
        "plaac_qn_frac": round(qn_frac, 6),
        "plaac_max_window_qn": round(max_w, 6),
        "plaac_nllr_approx": round(nllr_approx, 6),
        "plaac_core_len": float(best),
    }
