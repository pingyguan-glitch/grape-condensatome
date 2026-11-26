"""PlantLLPS scoring models.

Ships two scorable modes:
1. companion_composite — transparent weighted blend of disorder / LC / blockiness
   (weights 0.55 / 0.25 / 0.20 from the companion MAP study).
2. logistic — L2 logistic regression (stdlib SGD) trained on LLPSDB plant
   entries + companion gold standards when available.

Holdout metrics are computed via stratified k-fold; plant label scarcity means
AUCs are informative but not definitive. metapredict/pLDDT still pending.
"""

from __future__ import annotations

import json
import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .features import FEATURE_NAMES, extract_features, features_to_vector
from .metrics import roc_auc, spearman_rho, stratified_kfold_indices

# Companion v2 composite feature axes (normalized 0–1 at score time)
COMPOSITE_WEIGHTS = {
    "frac_disorder_proxy": 0.55,
    "low_complexity_score": 0.25,
    "blockiness_ncpr_sd11": 0.20,
}


def _sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def _clip01(x: float) -> float:
    return max(0.0, min(1.0, x))


@dataclass
class PlantLLPSModel:
    """Serializable PlantLLPS model."""

    mode: str = "hybrid"  # companion_composite | logistic_stub | hybrid
    feature_names: List[str] = field(default_factory=lambda: list(FEATURE_NAMES))
    weights: List[float] = field(default_factory=list)
    bias: float = 0.0
    means: List[float] = field(default_factory=list)
    stds: List[float] = field(default_factory=list)
    trained: bool = False
    train_note: str = ""
    metrics: Dict[str, Any] = field(default_factory=dict)
    # Optional proteome-calibrated cut-points (set after proteome scoring)
    tier_cutoffs: Dict[str, float] = field(default_factory=dict)

    def fit_logistic(
        self,
        X: Sequence[Sequence[float]],
        y: Sequence[int],
        lr: float = 0.15,
        epochs: int = 400,
        l2: float = 0.05,
        seed: int = 7,
    ) -> "PlantLLPSModel":
        if not X:
            raise ValueError("empty training matrix")
        rng = random.Random(seed)
        n_feat = len(X[0])
        self.feature_names = list(self.feature_names[:n_feat]) or list(FEATURE_NAMES[:n_feat])
        cols = list(zip(*X))
        self.means = [sum(c) / len(c) for c in cols]
        self.stds = []
        for c, m in zip(cols, self.means):
            var = sum((v - m) ** 2 for v in c) / max(1, len(c) - 1)
            self.stds.append(math.sqrt(var) if var > 1e-12 else 1.0)

        def scale_row(row: Sequence[float]) -> List[float]:
            return [(v - m) / s for v, m, s in zip(row, self.means, self.stds)]

        Xs = [scale_row(r) for r in X]
        w = [rng.uniform(-0.05, 0.05) for _ in range(n_feat)]
        b = 0.0
        idx = list(range(len(Xs)))
        for _ in range(epochs):
            rng.shuffle(idx)
            for i in idx:
                xi = Xs[i]
                yi = y[i]
                z = sum(wj * xj for wj, xj in zip(w, xi)) + b
                p = _sigmoid(z)
                err = p - yi
                for j in range(n_feat):
                    w[j] -= lr * (err * xi[j] + l2 * w[j])
                b -= lr * err
        self.weights = w
        self.bias = b
        self.trained = True
        self.mode = "hybrid"
        return self

    def companion_score(self, feats: Mapping[str, float]) -> float:
        """Transparent 0–1 composite (companion-inspired weights)."""
        # blockiness typically ~0–0.2; rescale for blend
        d = _clip01(float(feats.get("frac_disorder_proxy", 0.0)))
        lc = _clip01(float(feats.get("low_complexity_score", 0.0)))
        blk = _clip01(float(feats.get("blockiness_ncpr_sd11", 0.0)) / 0.2)
        # MAP65 grammar bonus (composition-level rescue → Tier B signal)
        bonus = 0.08 * float(feats.get("map65_tail_like", 0.0))
        bonus += 0.05 * _clip01(float(feats.get("c_term_ser_frac_100", 0.0)) / 0.15)
        score = (
            COMPOSITE_WEIGHTS["frac_disorder_proxy"] * d
            + COMPOSITE_WEIGHTS["low_complexity_score"] * lc
            + COMPOSITE_WEIGHTS["blockiness_ncpr_sd11"] * blk
            + bonus
        )
        return _clip01(score)

    def logistic_score(self, feats: Mapping[str, float]) -> float:
        if not self.trained or not self.weights:
            return self.companion_score(feats)
        vec = features_to_vector(feats, self.feature_names)
        scaled = [(v - m) / s for v, m, s in zip(vec, self.means, self.stds)]
        z = sum(wj * xj for wj, xj in zip(self.weights, scaled)) + self.bias
        return _sigmoid(z)

    def score(self, sequence: str) -> Dict[str, float]:
        feats = extract_features(sequence)
        comp = self.companion_score(feats)
        logi = self.logistic_score(feats) if self.trained else comp
        if self.mode == "companion_composite":
            plant = comp
        elif self.mode == "logistic_stub":
            plant = logi
        else:
            # hybrid: average companion prior with logistic stub when trained
            plant = 0.5 * comp + 0.5 * logi if self.trained else comp
        # Ser/Arg / MAP65-aware lift already partially in companion; expose parts
        out = dict(feats)
        out["score_companion"] = comp
        out["score_logistic"] = logi
        out["score_plantllps"] = plant
        out["tier"] = assign_tier(plant, feats, self.tier_cutoffs or None)
        return out

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2))

    @classmethod
    def load(cls, path: str | Path) -> "PlantLLPSModel":
        from dataclasses import fields as dc_fields

        data = json.loads(Path(path).read_text())
        known = {f.name for f in dc_fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in known})


def assign_tier(
    score: float,
    feats: Mapping[str, float],
    cutoffs: Optional[Mapping[str, float]] = None,
) -> str:
    """Tier labels A–D.

    If proteome `cutoffs` are provided (score_A / score_C thresholds from
    percentiles), use them for A/C/D global strata; MAP65-like grammar still
    forces Tier B when disorder is not extreme (paper design).
    """
    map65 = float(feats.get("map65_tail_like", 0.0))
    disorder = float(feats.get("frac_disorder_proxy", 0.0))
    c_ser = float(feats.get("c_term_ser_frac_100", 0.0))
    if map65 >= 0.5 and disorder < 0.70:
        return "B"
    if score >= 0.45 and c_ser >= 0.10 and disorder < 0.55 and map65 > 0:
        return "B"

    if cutoffs:
        a_cut = float(cutoffs.get("score_A", 0.65))
        c_cut = float(cutoffs.get("score_C", 0.40))
        if score >= a_cut and disorder >= float(cutoffs.get("disorder_A", 0.45)):
            return "A"
        if score >= c_cut:
            return "C"
        return "D"

    if score >= 0.65 and disorder >= 0.45:
        return "A"
    if score >= 0.40:
        return "C"
    return "D"


def load_training_csv(path: str | Path) -> Tuple[List[str], List[str], List[int], List[str]]:
    """Parse training CSV → ids, sequences, labels (1/0), notes."""
    import csv

    ids: List[str] = []
    seqs: List[str] = []
    labels: List[int] = []
    notes: List[str] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            ids.append(row["protein_id"])
            seqs.append(row["sequence"].strip())
            labels.append(1 if row["label"].strip().lower().startswith("pos") else 0)
            notes.append(row.get("note", ""))
    return ids, seqs, labels, notes


def cross_validate(
    X: Sequence[Sequence[float]],
    y: Sequence[int],
    k: int = 5,
    seed: int = 7,
) -> Dict[str, Any]:
    """Stratified k-fold AUC for logistic and companion scores."""
    folds = stratified_kfold_indices(y, k=k, seed=seed)
    fold_aucs_log: List[float] = []
    fold_aucs_comp: List[float] = []
    oof_log = [float("nan")] * len(y)
    oof_comp = [float("nan")] * len(y)
    for train_idx, test_idx in folds:
        Xtr = [X[i] for i in train_idx]
        ytr = [y[i] for i in train_idx]
        m = PlantLLPSModel()
        m.fit_logistic(Xtr, ytr, seed=seed)
        y_te = [y[i] for i in test_idx]
        s_log = []
        s_comp = []
        for i in test_idx:
            # reconstruct feats dict from vector for companion
            feats = {n: v for n, v in zip(FEATURE_NAMES, X[i])}
            sl = m.logistic_score(feats)
            sc = m.companion_score(feats)
            s_log.append(sl)
            s_comp.append(sc)
            oof_log[i] = sl
            oof_comp[i] = sc
        fold_aucs_log.append(roc_auc(y_te, s_log))
        fold_aucs_comp.append(roc_auc(y_te, s_comp))

    def _mean(xs: List[float]) -> float:
        ok = [x for x in xs if x == x]
        return sum(ok) / len(ok) if ok else float("nan")

    # OOF metrics where defined
    oof_pairs_log = [(oof_log[i], y[i]) for i in range(len(y)) if oof_log[i] == oof_log[i]]
    oof_pairs_comp = [(oof_comp[i], y[i]) for i in range(len(y)) if oof_comp[i] == oof_comp[i]]
    oof_auc_log = roc_auc([p[1] for p in oof_pairs_log], [p[0] for p in oof_pairs_log])
    oof_auc_comp = roc_auc([p[1] for p in oof_pairs_comp], [p[0] for p in oof_pairs_comp])
    rho = spearman_rho(
        [p[0] for p in oof_pairs_log],
        [p[0] for p in oof_pairs_comp],
    )
    return {
        "n": len(y),
        "n_pos": int(sum(y)),
        "n_neg": int(len(y) - sum(y)),
        "k": k,
        "fold_auc_logistic": fold_aucs_log,
        "fold_auc_companion": fold_aucs_comp,
        "mean_auc_logistic": _mean(fold_aucs_log),
        "mean_auc_companion": _mean(fold_aucs_comp),
        "oof_auc_logistic": oof_auc_log,
        "oof_auc_companion": oof_auc_comp,
        "spearman_logistic_vs_companion": rho,
    }


def train_from_csv(
    csv_path: str | Path,
    seed: int = 7,
    k_folds: int = 5,
    evaluate: bool = True,
) -> PlantLLPSModel:
    ids, seqs, labels, _ = load_training_csv(csv_path)
    X = [features_to_vector(extract_features(s)) for s in seqs]
    model = PlantLLPSModel()
    metrics: Dict[str, Any] = {}
    if evaluate and len(ids) >= 8:
        metrics = cross_validate(X, labels, k=min(k_folds, sum(labels), len(labels) - sum(labels)), seed=seed)
    model.fit_logistic(X, labels, seed=seed)
    model.metrics = metrics
    model.train_note = (
        f"Trained on {len(ids)} sequences from {Path(csv_path).name} "
        f"(pos={sum(labels)}, neg={len(labels)-sum(labels)}). "
        "Corpus: LLPSDB plant entries (mechanism-filtered) + companion/literature plant "
        "drivers (incl. FLOE1 Q8VZR8) + crop MAP65 orthologs + animal controls + "
        "hard folded plant negatives; 90% identity dedup. "
        "CRY2 held out as mechanism-special. Disorder is composition-proxy (not metapredict)."
    )
    return model


def default_model_path(package_root: Optional[Path] = None) -> Path:
    root = package_root or Path(__file__).resolve().parent.parent
    return root / "data" / "model_formal.json"


def ensure_trained_model(package_root: Optional[Path] = None) -> PlantLLPSModel:
    root = package_root or Path(__file__).resolve().parent.parent
    for name in ("model_formal.json", "model_stub.json"):
        path = root / "data" / name
        if path.exists():
            return PlantLLPSModel.load(path)
    for csv_name in ("train_formal.csv", "train_tiny.csv"):
        train_csv = root / "data" / csv_name
        if train_csv.exists():
            model = train_from_csv(train_csv)
            model.save(root / "data" / "model_formal.json")
            return model
    return PlantLLPSModel(mode="companion_composite")


def percentile_cutoffs(scores: Sequence[float], pct_a: float = 95.0, pct_c: float = 50.0) -> Dict[str, float]:
    """Score cutoffs from proteome percentiles (pct_a = Tier A floor)."""
    if not scores:
        return {"score_A": 0.65, "score_C": 0.40, "disorder_A": 0.45}
    s = sorted(scores)
    def pct(p: float) -> float:
        if not s:
            return 0.0
        i = min(len(s) - 1, max(0, int(round((p / 100.0) * (len(s) - 1)))))
        return s[i]
    return {"score_A": pct(pct_a), "score_C": pct(pct_c), "disorder_A": 0.45}
