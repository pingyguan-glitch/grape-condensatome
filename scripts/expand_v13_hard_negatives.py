#!/usr/bin/env python3
"""v13: add curated disordered non-LLPS hard negatives; retrain; report delta vs v12.

Hard negatives are plant structural IDPs (extensins, PRPs, AGPs, LRX, linker
histones) plus human prothymosin-α — high composition-proxy disorder, no
condensate / LLPS training labels. Sequences resolved via UniProt REST and
cached under /tmp/plantllps_v13/seqs (or data/hard_neg_v13/).
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plantllps.features import extract_features  # noqa: E402
from plantllps.model import train_from_csv  # noqa: E402

DATA = ROOT / "data"
SEQ_DIRS = [
    ROOT / "data" / "hard_neg_v13" / "seqs",
    Path("/tmp/plantllps_v13/seqs"),
]

# (accession, kingdom, note) — curated; filter again by disorder + 90% dedup
HARD_NEGS: list[tuple[str, str, str]] = [
    ("Q38913", "plant", "EXT1 hydroxyproline-rich extensin; structural cell-wall IDP; no LLPS evidence"),
    ("A0A1P8AM70", "plant", "EXT15 proline-rich glycoprotein; extracellular IDP"),
    ("Q9FZ35", "plant", "PRP1 proline-rich protein; structural IDP"),
    ("Q9LZJ7", "plant", "PRP3 proline-rich protein; structural IDP"),
    ("Q9XIP3", "plant", "PRP-family proline-rich; structural IDP"),
    ("Q9FZA2", "plant", "AGP31 arabinogalactan protein; extracellular glycosylated IDP"),
    ("O22194", "plant", "AGP17 arabinogalactan protein; extracellular IDP"),
    ("Q9ZT16", "plant", "AGP4 arabinogalactan protein; extracellular IDP"),
    ("P26569", "plant", "H1.2 linker histone; disordered chromatin protein; not IDR-driver condensate label"),
    ("F4IPK9", "plant", "HIS1-3 linker histone; disordered chromatin protein"),
    ("Q3ECQ3", "plant", "H1-like linker histone; disordered chromatin protein"),
    ("Q9T0K5", "plant", "LRX3 leucine-rich repeat/extensin; wall-associated IDP"),
    ("O65375", "plant", "LRX1 leucine-rich repeat/extensin; wall-associated IDP"),
    ("F4J0B5", "plant", "EXT16 proline-rich extensin-like; wall IDP"),
    ("Q9STN0", "plant", "EXT22 proline-rich extensin-like; wall IDP"),
    ("A0ACS1AHG4", "plant", "EXT9 proline-rich extensin-like; wall IDP"),
    ("Q9LHF1", "plant", "LRX4 leucine-rich repeat/extensin; wall-associated"),
    ("Q9FS16", "plant", "EXT3 extensin; structural cell-wall IDP"),
    ("P06454", "animal", "human prothymosin-alpha; classic highly disordered non-condensing IDP control"),
]


def load_seq(acc: str) -> str | None:
    for d in SEQ_DIRS:
        p = d / f"{acc}.fasta"
        if not p.exists():
            continue
        lines = p.read_text().splitlines()
        seq = "".join(l.strip() for l in lines if l and not l.startswith(">")).upper()
        if len(seq) >= 50:
            return seq
    return None


def identity(a: str, b: str) -> float:
    if a == b:
        return 1.0
    la, lb = len(a), len(b)
    if min(la, lb) / max(la, lb) < 0.80:
        return 0.0
    if la > lb:
        a, b = b, a
        la, lb = lb, la
    best = 0
    for off in range(lb - la + 1):
        window = b[off : off + la]
        match = sum(1 for x, y in zip(a, window) if x == y)
        if match > best:
            best = match
            if best / la >= 0.98:
                break
    return best / la


def dedup(rows: list[dict], thresh: float = 0.90) -> tuple[list[dict], list[tuple[str, str]]]:
    rows = sorted(
        rows,
        key=lambda r: (
            0 if r["label"] == "positive" else 1,
            0 if r["source"] != "round_v13_hard_neg" else 1,
            -len(r["sequence"]),
            r["protein_id"],
        ),
    )
    kept: list[dict] = []
    dropped: list[tuple[str, str]] = []
    for r in rows:
        dup = None
        for k in kept:
            if k["label"] != r["label"]:
                continue
            if identity(r["sequence"], k["sequence"]) >= thresh:
                dup = k["protein_id"]
                break
        if dup:
            dropped.append((r["protein_id"], dup))
        else:
            kept.append(r)
    return kept, dropped


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = ["protein_id", "label", "kingdom", "source", "note", "sequence"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["label"] != "positive", x["source"], x["protein_id"])):
            w.writerow({f: r[f] for f in fields})


def main() -> None:
    v12_path = DATA / "train_formal_v12.csv"
    if not v12_path.exists():
        # fall back to training/ mirror
        alt = DATA / "training" / "train_formal_v12.csv"
        if not alt.exists():
            raise SystemExit(f"missing {v12_path}")
        v12_path = alt
        (DATA / "train_formal_v12.csv").write_text(v12_path.read_text())
        v12_path = DATA / "train_formal_v12.csv"

    base = list(csv.DictReader(v12_path.open()))
    for r in base:
        r["sequence"] = r["sequence"].strip().upper()
    current = {r["protein_id"]: r for r in base}

    newly: list[dict] = []
    skipped: list[str] = []
    for acc, kingdom, note in HARD_NEGS:
        if acc in current:
            skipped.append(f"{acc}:already_in_set")
            continue
        seq = load_seq(acc)
        if not seq:
            skipped.append(f"{acc}:missing_seq")
            continue
        feats = extract_features(seq)
        dis = float(feats["frac_disorder_proxy"])
        if dis < 0.35:
            skipped.append(f"{acc}:low_disorder={dis:.3f}")
            continue
        newly.append(
            {
                "protein_id": acc,
                "label": "negative",
                "kingdom": kingdom,
                "source": "round_v13_hard_neg",
                "note": f"{note} [disorder_proxy={dis:.3f}]".replace(",", ";"),
                "sequence": seq,
            }
        )

    merged = list(current.values()) + newly
    kept, dropped = dedup(merged, 0.90)
    rows = kept
    csv_path = DATA / "train_formal_v13.csv"
    write_csv(csv_path, rows)
    # mirrors
    (DATA / "train_formal.csv").write_text(csv_path.read_text())
    train_dir = DATA / "training"
    train_dir.mkdir(exist_ok=True)
    (train_dir / "train_formal_v13.csv").write_text(csv_path.read_text())
    (train_dir / "train_formal.csv").write_text(csv_path.read_text())

    # baseline v12 metrics (recompute for fair delta on same code path)
    model_v12 = train_from_csv(DATA / "train_formal_v12.csv", seed=7, k_folds=5, evaluate=True)
    model_v13 = train_from_csv(csv_path, seed=7, k_folds=5, evaluate=True)
    model_v13.train_note = (
        "v13 formal set: v12 + curated disordered non-LLPS hard negatives "
        "(extensin/PRP/AGP/LRX/linker histone + prothymosin-α). "
        "Disorder = composition proxy (not metapredict)."
    )
    model_v13.save(DATA / "model_formal.json")
    (DATA / "model_formal.metrics.json").write_text(json.dumps(model_v13.metrics, indent=2))
    (DATA / "model_formal.metrics.v13.json").write_text(json.dumps(model_v13.metrics, indent=2))
    (train_dir / "model_formal.metrics.v13.json").write_text(json.dumps(model_v13.metrics, indent=2))

    n_pos = sum(1 for r in rows if r["label"] == "positive")
    n_neg = sum(1 for r in rows if r["label"] == "negative")
    added_kept = len([r for r in newly if r["protein_id"] in {x["protein_id"] for x in rows}])

    delta = {
        "version": 13,
        "theme": "Disordered non-LLPS hard negatives (extensin/PRP/AGP/LRX/H1 + PTMA)",
        "n": len(rows),
        "n_pos": n_pos,
        "n_neg": n_neg,
        "added_attempted": len(HARD_NEGS),
        "added_raw": len(newly),
        "added_kept": added_kept,
        "skipped": skipped,
        "dropped_identity": dropped,
        "v12": {
            "n": model_v12.metrics["n"],
            "n_pos": model_v12.metrics["n_pos"],
            "n_neg": model_v12.metrics["n_neg"],
            "mean_auc_logistic": model_v12.metrics["mean_auc_logistic"],
            "mean_auc_companion": model_v12.metrics["mean_auc_companion"],
            "oof_auc_logistic": model_v12.metrics["oof_auc_logistic"],
            "oof_auc_companion": model_v12.metrics.get("oof_auc_companion"),
        },
        "v13": {
            "mean_auc_logistic": model_v13.metrics["mean_auc_logistic"],
            "mean_auc_companion": model_v13.metrics["mean_auc_companion"],
            "oof_auc_logistic": model_v13.metrics["oof_auc_logistic"],
            "oof_auc_companion": model_v13.metrics.get("oof_auc_companion"),
            "spearman": model_v13.metrics.get("spearman_logistic_vs_companion"),
        },
        "delta_mean_auc_logistic": model_v13.metrics["mean_auc_logistic"]
        - model_v12.metrics["mean_auc_logistic"],
        "delta_oof_auc_logistic": model_v13.metrics["oof_auc_logistic"]
        - model_v12.metrics["oof_auc_logistic"],
        "delta_mean_auc_companion": model_v13.metrics["mean_auc_companion"]
        - model_v12.metrics["mean_auc_companion"],
        "hard_neg_ids": [r["protein_id"] for r in newly if r["protein_id"] in {x["protein_id"] for x in rows}],
        "honesty": (
            "Adding hard disordered negatives typically lowers AUC vs folded-heavy v12; "
            "that drop is expected and more informative for submission claims."
        ),
    }
    out = DATA / "train_v13_hard_neg.delta.json"
    out.write_text(json.dumps(delta, indent=2))
    (train_dir / "train_v13_hard_neg.delta.json").write_text(json.dumps(delta, indent=2))

    print(
        f"v13: n={delta['n']} (+{added_kept} hard neg) pos={n_pos} neg={n_neg}\n"
        f"  AUC_L  v12={delta['v12']['mean_auc_logistic']:.3f} → v13={delta['v13']['mean_auc_logistic']:.3f} "
        f"(Δ={delta['delta_mean_auc_logistic']:+.3f})\n"
        f"  OOF_L  v12={delta['v12']['oof_auc_logistic']:.3f} → v13={delta['v13']['oof_auc_logistic']:.3f} "
        f"(Δ={delta['delta_oof_auc_logistic']:+.3f})\n"
        f"  AUC_C  v12={delta['v12']['mean_auc_companion']:.3f} → v13={delta['v13']['mean_auc_companion']:.3f} "
        f"(Δ={delta['delta_mean_auc_companion']:+.3f})",
        flush=True,
    )
    if skipped:
        print("  skipped:", skipped, flush=True)
    if dropped:
        print("  dedup:", dropped, flush=True)
    print(f"Wrote {csv_path}, {out}, data/model_formal.json", flush=True)


if __name__ == "__main__":
    main()
