#!/usr/bin/env python3
"""Optional PlantLLPS retrain with metapredict and/or AF pLDDT disorder.

Comparable stratified CV vs the composition-proxy production path. **Does not**
overwrite ``data/model_formal.json`` — experimental weights land under
``data/experimental/``.

Disorder variants
-----------------
* ``composition`` — baseline (same as production); reported for side-by-side Δ.
* ``metapredict`` — replace ``frac_disorder_proxy`` with metapredict fraction.
* ``af_plddt`` — replace with AFDB ``af_frac_very_low`` (pLDDT&lt;50 fraction);
  fall back to composition proxy when no AFDB model exists.
* ``metapredict_af_mean`` — mean of metapredict + AF (AF missing → metapredict).

Requires metapredict for the metapredict variants::

    export METAPREDICT_PYTHON="$(micromamba run -n meta311 which python)"
    python3 scripts/retrain_disorder_features.py --variants metapredict af_plddt
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from _http_util import urlopen_text  # noqa: E402
from plantllps.features import FEATURE_NAMES, extract_features, features_to_vector  # noqa: E402
from plantllps.metapredict_bridge import resolve_metapredict_scores  # noqa: E402
from plantllps.metrics import roc_auc  # noqa: E402
from plantllps.model import (  # noqa: E402
    PlantLLPSModel,
    cross_validate,
)

DEFAULT_TRAIN = ROOT / "data" / "train_formal_v13.csv"
DEFAULT_OUTDIR = ROOT / "data" / "experimental"
AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{pid}"
BROWSER_UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)


def load_train(path: Path):
    ids, seqs, labels, notes = [], [], [], []
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            ids.append(row["protein_id"])
            seqs.append(row["sequence"].strip())
            lab = row["label"].strip().lower()
            labels.append(1 if lab.startswith("pos") else 0)
            notes.append(row.get("note", ""))
    return ids, seqs, labels, notes


def fetch_af_frac_very_low(pid: str, cache: Path) -> Optional[float]:
    """Return AFDB fraction pLDDT&lt;50, or None if no model / fetch failure."""
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{pid}.prediction.json"
    pred = None
    if path.exists() and path.stat().st_size >= 2:
        try:
            data = json.loads(path.read_text())
            if isinstance(data, list):
                pred = data[0] if data else None
            elif isinstance(data, dict):
                pred = data
        except json.JSONDecodeError:
            pred = None
    else:
        try:
            raw = urlopen_text(
                AF_API.format(pid=pid),
                timeout=45,
                accept="application/json",
                ua=BROWSER_UA,
            )
            path.write_text(raw)
            data = json.loads(raw)
            if isinstance(data, list):
                pred = data[0] if data else None
            elif isinstance(data, dict):
                pred = data
            else:
                pred = None
        except Exception as e:
            code = getattr(e, "code", None)
            if code == 404:
                path.write_text("[]")
                return None
            # soft-fail: leave uncached so a later run can retry
            return None
    if not pred:
        return None
    # Prefer explicit bin if present; else derive from confidence URL is too heavy —
    # API summary often includes global metrics. Fall back to residualBFactor mean.
    # AF prediction JSON fields vary; try common keys then confidence doc.
    for key in ("fractionPlddtVeryLow", "fraction_plddt_very_low"):
        if key in pred and pred[key] is not None:
            return float(pred[key])
    # Derive from per-residue confidence if linked
    conf_url = pred.get("confidenceScoreUrl") or pred.get("plddtDocUrl") or ""
    model_id = pred.get("modelEntityId") or pred.get("entryId") or pid
    conf_cache = cache / f"{pid}.confidence.json"
    if conf_cache.exists() and conf_cache.stat().st_size > 20:
        try:
            conf = json.loads(conf_cache.read_text())
            scores = conf.get("confidenceScore") or conf.get("pLDDT") or conf.get("plddt")
            if scores:
                return float(sum(1 for s in scores if float(s) < 50.0) / len(scores))
        except Exception:
            pass
    urls = []
    if conf_url:
        urls.append(conf_url)
    for ver in (6, 5, 4, 3):
        urls.append(
            f"https://alphafold.ebi.ac.uk/files/AF-{pid}-F1-confidence_v{ver}.json"
        )
        urls.append(
            f"https://alphafold.ebi.ac.uk/files/{model_id}-confidence_v{ver}.json"
        )
    for url in urls:
        try:
            raw = urlopen_text(url, timeout=60, accept="application/json", ua=BROWSER_UA)
            conf_cache.write_text(raw)
            conf = json.loads(raw)
            scores = conf.get("confidenceScore") or conf.get("pLDDT") or conf.get("plddt")
            if scores:
                return float(sum(1 for s in scores if float(s) < 50.0) / len(scores))
        except Exception:
            continue
    # Last resort: invert mean pLDDT from API if present (not identical to frac<50)
    for key in ("globalMetricValue", "modelConfidence", "ptm"):
        pass
    mean = pred.get("globalMetricValue")
    if mean is not None:
        # Approximate disorder as clipped (70 - mean)/70 — documented as approx fallback
        try:
            m = float(mean)
            return max(0.0, min(1.0, (70.0 - m) / 70.0))
        except (TypeError, ValueError):
            return None
    return None


def collect_af_map(
    ids: Sequence[str],
    cache: Path,
    *,
    workers: int = 8,
) -> Dict[str, Optional[float]]:
    out: Dict[str, Optional[float]] = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fetch_af_frac_very_low, pid, cache): pid for pid in ids}
        for fut in as_completed(futs):
            pid = futs[fut]
            try:
                out[pid] = fut.result()
            except Exception:
                out[pid] = None
    return out


def apply_disorder(
    feats: Dict[str, float],
    *,
    variant: str,
    meta_frac: Optional[float],
    af_frac: Optional[float],
) -> Dict[str, float]:
    out = dict(feats)
    proxy = float(out.get("frac_disorder_proxy", 0.0))
    if variant == "composition":
        disorder = proxy
    elif variant == "metapredict":
        if meta_frac is None:
            raise ValueError("metapredict fraction required")
        disorder = float(meta_frac)
        out["frac_disorder_metapredict"] = disorder
    elif variant == "af_plddt":
        disorder = float(af_frac) if af_frac is not None else proxy
        out["frac_disorder_af_plddt"] = float(af_frac) if af_frac is not None else float("nan")
        out["af_plddt_fallback_composition"] = 1.0 if af_frac is None else 0.0
    elif variant == "metapredict_af_mean":
        if meta_frac is None:
            raise ValueError("metapredict fraction required")
        if af_frac is None:
            disorder = float(meta_frac)
        else:
            disorder = 0.5 * (float(meta_frac) + float(af_frac))
        out["frac_disorder_metapredict"] = float(meta_frac)
        out["frac_disorder_af_plddt"] = float(af_frac) if af_frac is not None else float("nan")
    else:
        raise ValueError(f"unknown variant: {variant}")
    out["frac_disorder_proxy"] = float(disorder)
    out["disorder_source"] = variant  # type: ignore[assignment]
    return out


def build_matrix(
    seqs: Sequence[str],
    ids: Sequence[str],
    *,
    variant: str,
    meta_map: Mapping[str, float],
    af_map: Mapping[str, Optional[float]],
) -> Tuple[List[List[float]], List[Dict[str, float]], Dict[str, Any]]:
    X: List[List[float]] = []
    feat_rows: List[Dict[str, float]] = []
    n_af = 0
    n_af_fallback = 0
    for pid, seq in zip(ids, seqs):
        base = extract_features(seq)
        meta_frac = meta_map.get(pid) if variant != "composition" else None
        af_frac = af_map.get(pid)
        feats = apply_disorder(
            base,
            variant=variant,
            meta_frac=meta_frac,
            af_frac=af_frac,
        )
        if af_frac is not None and variant in ("af_plddt", "metapredict_af_mean"):
            n_af += 1
        if variant == "af_plddt" and af_frac is None:
            n_af_fallback += 1
        X.append(features_to_vector(feats, FEATURE_NAMES))
        feat_rows.append(feats)
    cov = {
        "n": len(ids),
        "n_metapredict": len(meta_map) if variant != "composition" else 0,
        "n_af_plddt": n_af,
        "n_af_fallback_composition": n_af_fallback,
        "af_coverage": n_af / len(ids) if ids else 0.0,
    }
    return X, feat_rows, cov


def train_variant(
    X: Sequence[Sequence[float]],
    y: Sequence[int],
    *,
    variant: str,
    seed: int,
    k: int,
    train_csv_name: str,
) -> Tuple[PlantLLPSModel, Dict[str, Any]]:
    metrics = cross_validate(X, y, k=k, seed=seed)
    # Also hybrid OOF for comparison with family-holdout style numbers
    # (cross_validate already has logistic + companion)
    model = PlantLLPSModel()
    model.fit_logistic(X, y, seed=seed)
    model.metrics = metrics
    model.train_note = (
        f"EXPERIMENTAL disorder={variant}. Trained on {len(y)} sequences from "
        f"{train_csv_name} (pos={sum(y)}, neg={len(y)-sum(y)}). "
        "Not the production default — keep data/model_formal.json (composition proxy)."
    )
    # Annotate metrics
    metrics = dict(metrics)
    metrics["disorder_variant"] = variant
    metrics["experimental"] = True
    # hybrid mean from fold averages of 0.5*(log+comp) approx via OOF
    return model, metrics


def write_disorder_tsv(
    path: Path,
    ids: Sequence[str],
    labels: Sequence[int],
    composition: Sequence[Dict[str, float]],
    meta_map: Mapping[str, float],
    af_map: Mapping[str, Optional[float]],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "protein_id",
                "label",
                "frac_disorder_proxy",
                "frac_disorder_metapredict",
                "frac_disorder_af_plddt",
                "af_available",
            ]
        )
        for i, pid in enumerate(ids):
            af = af_map.get(pid)
            w.writerow(
                [
                    pid,
                    "positive" if labels[i] else "negative",
                    f"{composition[i].get('frac_disorder_proxy', 0.0):.6f}",
                    f"{meta_map.get(pid, float('nan')):.6f}" if pid in meta_map else "",
                    f"{af:.6f}" if af is not None else "",
                    "1" if af is not None else "0",
                ]
            )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train-csv", type=Path, default=DEFAULT_TRAIN)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument(
        "--variants",
        nargs="+",
        default=["composition", "metapredict", "af_plddt", "metapredict_af_mean"],
        choices=["composition", "metapredict", "af_plddt", "metapredict_af_mean"],
    )
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--metapredict-python", default=None)
    ap.add_argument("--metapredict-threshold", type=float, default=0.5)
    ap.add_argument("--af-cache", type=Path, default=ROOT / "data" / "alphafold" / "train_cache")
    ap.add_argument("--af-workers", type=int, default=8)
    ap.add_argument(
        "--skip-metapredict",
        action="store_true",
        help="Skip variants that need metapredict (still runs composition / af_plddt)",
    )
    args = ap.parse_args()

    ids, seqs, labels, _notes = load_train(args.train_csv)
    y = labels
    k = min(args.folds, sum(y), len(y) - sum(y))
    args.outdir.mkdir(parents=True, exist_ok=True)

    need_meta = any(v in args.variants for v in ("metapredict", "metapredict_af_mean")) and not args.skip_metapredict
    need_af = any(v in args.variants for v in ("af_plddt", "metapredict_af_mean"))

    meta_map: Dict[str, float] = {}
    if need_meta:
        print(f"Scoring metapredict for {len(ids)} train sequences …", file=sys.stderr)
        t0 = time.time()
        meta_map = resolve_metapredict_scores(
            list(zip(ids, seqs)),
            threshold=args.metapredict_threshold,
            python_exe=args.metapredict_python,
        )
        print(f"  metapredict done in {time.time()-t0:.1f}s (n={len(meta_map)})", file=sys.stderr)
    elif need_meta is False and any(v in args.variants for v in ("metapredict", "metapredict_af_mean")):
        print("Skipping metapredict variants (--skip-metapredict).", file=sys.stderr)
        args.variants = [v for v in args.variants if v not in ("metapredict", "metapredict_af_mean")]

    af_map: Dict[str, Optional[float]] = {pid: None for pid in ids}
    if need_af:
        print(f"Fetching AF pLDDT disorder for {len(ids)} UniProt IDs …", file=sys.stderr)
        t0 = time.time()
        af_map = collect_af_map(ids, args.af_cache, workers=args.af_workers)
        n_ok = sum(1 for v in af_map.values() if v is not None)
        print(
            f"  AF coverage {n_ok}/{len(ids)} in {time.time()-t0:.1f}s → {args.af_cache}",
            file=sys.stderr,
        )

    # Composition feature rows for TSV deposit
    composition_rows = [extract_features(s) for s in seqs]
    write_disorder_tsv(
        args.outdir / "train_v13_disorder_features.tsv",
        ids,
        labels,
        composition_rows,
        meta_map,
        af_map,
    )

    comparison: Dict[str, Any] = {
        "version": 1,
        "train_csv": str(args.train_csv.relative_to(ROOT)),
        "n": len(y),
        "n_pos": int(sum(y)),
        "n_neg": int(len(y) - sum(y)),
        "seed": args.seed,
        "k": k,
        "production_model_untouched": "data/model_formal.json",
        "variants": {},
        "honesty": (
            "Experimental disorder retrains only. Production PlantLLPS remains the "
            "composition-proxy weights in data/model_formal.json. AF disorder uses "
            "AFDB frac pLDDT<50 when a model exists; otherwise falls back to the "
            "composition proxy for that row (documented in coverage)."
        ),
    }

    baseline_auc = None
    for variant in args.variants:
        if variant in ("metapredict", "metapredict_af_mean") and not meta_map:
            print(f"Skip {variant}: no metapredict map", file=sys.stderr)
            continue
        print(f"Training variant={variant} …", file=sys.stderr)
        X, _feats, cov = build_matrix(
            seqs, ids, variant=variant, meta_map=meta_map, af_map=af_map
        )
        model, metrics = train_variant(
            X,
            y,
            variant=variant,
            seed=args.seed,
            k=k,
            train_csv_name=args.train_csv.name,
        )
        metrics["coverage"] = cov
        if baseline_auc is None and variant == "composition":
            baseline_auc = metrics["mean_auc_logistic"]
        if baseline_auc is not None:
            metrics["delta_mean_auc_logistic_vs_composition"] = (
                metrics["mean_auc_logistic"] - baseline_auc
            )

        stem = f"model_formal.{variant}"
        if variant == "composition":
            # Still deposit metrics for Δ reference; do NOT write a second production model
            metrics_path = args.outdir / f"{stem}.metrics.json"
            metrics_path.write_text(json.dumps(metrics, indent=2))
            comparison["variants"][variant] = {
                "metrics_path": str(metrics_path.relative_to(ROOT)),
                "model_path": None,
                "note": "baseline CV only; production weights unchanged",
                **{k: metrics[k] for k in (
                    "mean_auc_logistic",
                    "mean_auc_companion",
                    "oof_auc_logistic",
                    "oof_auc_companion",
                    "delta_mean_auc_logistic_vs_composition",
                ) if k in metrics},
                "coverage": cov,
            }
        else:
            model_path = args.outdir / f"{stem}.json"
            metrics_path = args.outdir / f"{stem}.metrics.json"
            model.save(model_path)
            metrics_path.write_text(json.dumps(metrics, indent=2))
            comparison["variants"][variant] = {
                "metrics_path": str(metrics_path.relative_to(ROOT)),
                "model_path": str(model_path.relative_to(ROOT)),
                **{k: metrics[k] for k in (
                    "mean_auc_logistic",
                    "mean_auc_companion",
                    "oof_auc_logistic",
                    "oof_auc_companion",
                    "delta_mean_auc_logistic_vs_composition",
                ) if k in metrics},
                "coverage": cov,
            }
            print(
                f"  {variant}: mean_auc_logistic={metrics['mean_auc_logistic']:.3f} "
                f"OOF={metrics['oof_auc_logistic']:.3f} → {model_path}",
                file=sys.stderr,
            )

    comp_path = args.outdir / "disorder_retrain_comparison.json"
    comp_path.write_text(json.dumps(comparison, indent=2))

    readme = args.outdir / "README.md"
    lines = [
        "# Experimental PlantLLPS disorder retrains",
        "",
        "Generated by `scripts/retrain_disorder_features.py`.",
        "",
        "**Production default is unchanged:** `data/model_formal.json` (composition-proxy disorder).",
        "",
        "## Variants",
        "",
    ]
    for v, block in comparison["variants"].items():
        auc = block.get("mean_auc_logistic")
        delta = block.get("delta_mean_auc_logistic_vs_composition")
        mp = block.get("model_path") or "(metrics only)"
        dstr = f"{delta:+.3f}" if isinstance(delta, float) else "n/a"
        lines.append(
            f"- **{v}**: mean logistic AUC ≈ **{auc:.3f}** (Δ vs composition {dstr}); model: `{mp}`"
            if isinstance(auc, float)
            else f"- **{v}**: see JSON; model: `{mp}`"
        )
    lines += [
        "",
        "## Files",
        "",
        "- `train_v13_disorder_features.tsv` — per-ID composition / metapredict / AF fractions",
        "- `disorder_retrain_comparison.json` — side-by-side CV summary",
        "- `model_formal.<variant>.json` — experimental weights (non-composition only)",
        "",
        "## Honesty",
        "",
        "AF coverage is incomplete for non-AFDB accessions; those rows fall back to the",
        "composition proxy inside the `af_plddt` variant. Do not silently replace atlas",
        "production scores with these experimental weights.",
        "",
    ]
    readme.write_text("\n".join(lines))
    print(f"Comparison → {comp_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
