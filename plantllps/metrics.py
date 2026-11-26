"""Holdout metrics for PlantLLPS training (stdlib)."""

from __future__ import annotations

import math
from typing import List, Sequence, Tuple


def roc_auc(y_true: Sequence[int], y_score: Sequence[float]) -> float:
    """Mann–Whitney U ROC AUC. Returns nan if only one class present."""
    pairs = sorted(zip(y_score, y_true), key=lambda t: t[0])
    n_pos = sum(1 for y in y_true if y == 1)
    n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    # ranks (average for ties)
    ranks = [0.0] * len(pairs)
    i = 0
    while i < len(pairs):
        j = i
        while j < len(pairs) and pairs[j][0] == pairs[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0  # 1-based ranks
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j
    sum_ranks_pos = sum(r for r, (_, y) in zip(ranks, pairs) if y == 1)
    # AUC = (sum ranks_pos - n_pos*(n_pos+1)/2) / (n_pos*n_neg)
    return (sum_ranks_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def spearman_rho(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman rank correlation (average ranks for ties)."""
    if len(x) != len(y) or len(x) < 3:
        return float("nan")

    def rank(vals: Sequence[float]) -> List[float]:
        indexed = sorted(enumerate(vals), key=lambda t: t[1])
        ranks = [0.0] * len(vals)
        i = 0
        while i < len(indexed):
            j = i
            while j < len(indexed) and indexed[j][1] == indexed[i][1]:
                j += 1
            avg = (i + 1 + j) / 2.0
            for k in range(i, j):
                ranks[indexed[k][0]] = avg
            i = j
        return ranks

    rx, ry = rank(x), rank(y)
    n = len(x)
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    denx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    deny = math.sqrt(sum((b - my) ** 2 for b in ry))
    if denx < 1e-12 or deny < 1e-12:
        return float("nan")
    return num / (denx * deny)


def stratified_kfold_indices(
    y: Sequence[int], k: int = 5, seed: int = 7
) -> List[Tuple[List[int], List[int]]]:
    """Return list of (train_idx, test_idx) stratified folds."""
    import random

    rng = random.Random(seed)
    pos = [i for i, yi in enumerate(y) if yi == 1]
    neg = [i for i, yi in enumerate(y) if yi == 0]
    rng.shuffle(pos)
    rng.shuffle(neg)
    folds: List[Tuple[List[int], List[int]]] = []
    for fold in range(k):
        test = pos[fold::k] + neg[fold::k]
        test_set = set(test)
        train = [i for i in range(len(y)) if i not in test_set]
        if not test or not train:
            continue
        # need both classes in train when possible
        folds.append((train, test))
    return folds
