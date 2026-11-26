#!/usr/bin/env python3
"""Expand train_formal through v3–v12: add batches, 90% dedup, retrain, snapshot."""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plantllps.model import train_from_csv  # noqa: E402

DATA = ROOT / "data"
SEQ_DIRS = [
    Path("/tmp/plantllps_rounds/seqs"),
    Path("/tmp/plantllps_expand/seqs"),
    Path("/tmp/plantllps_formal/seqs"),
]

# Each round: theme + list of (accession, label, kingdom, note)
ROUNDS: dict[int, tuple[str, list[tuple[str, str, str, str]]]] = {
    3: (
        "Stress-granule / P-body plant RBPs + TCA-cycle folded negatives",
        [
            ("Q9SYG4", "positive", "plant", "UBP1A stress-granule polyU-binding"),
            ("Q8VZG7", "positive", "plant", "TSN1 Tudor-SN P-body/SG"),
            ("Q9LTT8", "positive", "plant", "VCS/EDC4 P-body scaffold"),
            ("A0ACS0ZFQ9", "positive", "plant", "DCP1 decapping coactivator"),
            ("Q9LXS6", "negative", "plant", "CSY2 citrate synthase folded"),
            ("Q42560", "negative", "plant", "ACO1 aconitase folded"),
            ("Q94B78", "negative", "plant", "ICDH isocitrate DH folded"),
            ("Q9SUS3", "negative", "plant", "MDH malate DH folded"),
        ],
    ),
    4: (
        "Photobody components + histone structured negatives",
        [
            ("P14713", "positive", "plant", "PHYB phytochrome B photobody"),
            ("P43254", "positive", "plant", "COP1 photomorphogenesis condensate"),
            ("Q9SSA7", "positive", "plant", "PIF4 phytochrome-interacting factor"),
            ("Q9S7N5", "positive", "plant", "SPA1 COP1-SPA complex"),
            ("P59226", "negative", "plant", "histone H3.1/HTR2 structured"),
            ("P59259", "negative", "plant", "histone H3 family structured"),
            ("Q9C5U8", "negative", "plant", "nucleosome-related structured"),
        ],
    ),
    5: (
        "Glycine-rich RBPs / YTH readers + glycolytic folded negatives",
        [
            ("Q03250", "positive", "plant", "GRP7/RBG7 cold-responsive GR-RBP"),
            ("Q03251", "positive", "plant", "GRP8/RBG8 glycine-rich RBP"),
            ("Q9ZQ38", "positive", "plant", "GRP2 glycine-rich RBP"),
            ("Q9LJE5", "positive", "plant", "ECT2 YTH m6A reader"),
            ("Q9LK35", "positive", "plant", "ECT3 YTH reader"),
            ("Q9M7A8", "negative", "plant", "PEP carboxylase PPC1 folded"),
            ("Q9SAU5", "negative", "plant", "fructose-bisphosphate aldolase folded"),
            ("Q9SKT0", "negative", "plant", "glucose-6-phosphate isomerase folded"),
            ("Q9SAU2", "negative", "plant", "aldolase isoform folded"),
        ],
    ),
    6: (
        "Evening-complex / cell-wall structured negatives",
        [
            ("Q9SNB4", "positive", "plant", "LUX/PCL1 evening complex"),
            ("Q06738", "positive", "plant", "RD29A/LTI78 stress-inducible"),
            ("O48946", "negative", "plant", "CESA1 cellulose synthase folded TM"),
            ("Q9FIR9", "negative", "plant", "expansin EXPA1 folded"),
            ("Q9M9A8", "negative", "plant", "cell-wall related structured"),
        ],
    ),
    7: (
        "MAP70 / TPX2 / SPIRAL1 plant MAPs + HSP chaperone folded negatives",
        [
            ("Q9C9X0", "positive", "plant", "AtMAP70.1 microtubule-associated"),
            ("F4I2H7", "positive", "plant", "AtTPX2 microtubule nucleation"),
            ("Q9SJW3", "positive", "plant", "SPIRAL1/SPR1 microtubule"),
            ("Q9LKR3", "negative", "plant", "HSP90.1 chaperone folded"),
            ("Q9STX5", "negative", "plant", "HSP90.2 chaperone folded"),
            ("P22954", "negative", "plant", "HSP70 cytosolic folded"),
            ("Q9LKR0", "negative", "plant", "HSP90 family folded"),
            ("Q9STW6", "negative", "plant", "HSP70-1 folded"),
            ("Q9LF37", "negative", "plant", "CCT chaperonin subunit folded"),
        ],
    ),
    8: (
        "Zinc-finger RBPs (RZ1) + MAPK folded negatives",
        [
            ("Q9LIN3", "positive", "plant", "RZ1A zinc-finger RNA-binding"),
            ("Q9FG31", "positive", "plant", "RZ-family RNA-binding"),
            ("Q9C5X8", "positive", "plant", "RZ1B-like RNA-binding"),
            ("Q9C5X9", "positive", "plant", "RZ1C-like RNA-binding"),
            ("Q39026", "negative", "plant", "MPK6 MAP kinase folded"),
            ("Q39027", "negative", "plant", "MPK4 MAP kinase folded"),
            ("Q9M7A0", "negative", "plant", "MAPK family folded"),
        ],
    ),
    9: (
        "Dehydrins / COR cold proteins + photosynthetic folded negatives",
        [
            ("P42759", "positive", "plant", "ERD10/LTI29 dehydrin"),
            ("P31168", "positive", "plant", "COR47/RD17 dehydrin"),
            ("P30185", "positive", "plant", "RAB18 dehydrin"),
            ("Q9SYI0", "positive", "plant", "COR15A cold-regulated"),
            ("O03042", "negative", "plant", "rbcL RuBisCO large chain folded"),
            ("P56761", "negative", "plant", "ATP synthase CF1 alpha folded"),
            ("P56757", "negative", "plant", "ATP synthase CF1 beta folded"),
            ("Q9SUI4", "negative", "plant", "photosystem I PsaA folded"),
        ],
    ),
    10: (
        "Mediator / HY5 transcription-related + proteostasis folded negatives",
        [
            ("Q7XYY2", "positive", "plant", "MED25/PFT1 mediator subunit"),
            ("Q9ASS3", "positive", "plant", "MED14 mediator subunit"),
            ("O24646", "positive", "plant", "HY5 bZIP photomorphogenesis"),
            ("Q9FZ09", "positive", "plant", "NF-YB transcription factor"),
            ("Q9FXA2", "positive", "plant", "PHYA phytochrome A"),
            ("Q9C5C2", "negative", "plant", "proteasome subunit folded"),
            ("Q9FNN5", "negative", "plant", "ubiquitin-conjugating enzyme folded"),
            ("Q9C5C5", "negative", "plant", "proteasome alpha folded"),
            ("Q9LEW3", "negative", "plant", "metacaspase folded"),
        ],
    ),
    11: (
        "Rice / soybean ortholog sweep + crop folded negatives",
        [
            ("Q10QA5", "positive", "plant", "rice glycine-rich RBP"),
            ("Q6Z7B0", "positive", "plant", "rice DEAD-box RNA helicase"),
            ("Q10CE7", "negative", "plant", "rice actin folded"),
            ("Q5N725", "negative", "plant", "rice histidine kinase folded"),
            ("Q7X8A1", "negative", "plant", "rice metabolic enzyme folded"),
            ("A0A0R0F8K1", "negative", "plant", "soybean GAPDH-like folded"),
        ],
    ),
    12: (
        "WVD2/WDL/IQD/MOR1 plant scaffolds + nuclear transport folded negatives",
        [
            ("Q84ZT9", "positive", "plant", "WVD2 WAVE-DAMPENED 2"),
            ("Q8GYX9", "positive", "plant", "WDL1 WVD2-like"),
            ("Q94FN2", "positive", "plant", "MOR1 microtubule organization"),
            ("Q9SF32", "positive", "plant", "IQD1 CaM-binding scaffold"),
            ("Q9SI19", "positive", "plant", "IQD22 IQ67-domain scaffold"),
            ("Q9C5C1", "negative", "plant", "NUP62 nuclear pore folded"),
            ("Q9C5C0", "negative", "plant", "NUP160 nuclear pore folded"),
            ("Q9LZF6", "negative", "plant", "importin-alpha folded"),
            ("Q9FJH7", "negative", "plant", "exportin folded"),
        ],
    ),
}


def load_seqs() -> dict[str, str]:
    seqs: dict[str, str] = {}
    for d in SEQ_DIRS:
        if not d.exists():
            continue
        for f in d.glob("*.fasta"):
            text = f.read_text()
            if not text.lstrip().startswith(">"):
                continue
            lines = text.splitlines()
            seq = "".join(l.strip() for l in lines if l and not l.startswith(">")).upper()
            if len(seq) >= 50:
                seqs[f.stem] = seq
                # also parse accession from header
                h = lines[0][1:]
                if "|" in h:
                    parts = h.split("|")
                    if len(parts) > 1:
                        seqs[parts[1]] = seq
    # existing training rows
    for name in ("train_formal_v2.csv", "train_formal.csv"):
        p = DATA / name
        if not p.exists():
            continue
        for row in csv.DictReader(p.open()):
            seqs[row["protein_id"]] = row["sequence"].strip().upper()
    return seqs


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
    preferred = {
        "Q9FLP0",
        "F6HET4",
        "O82804",
        "O04425",
        "Q7DZS7",
        "Q8W234",
        "Q8VZR8",
        "Q84TD8",
        "Q84ZT9",
        "Q9SF32",
    }
    rows = sorted(
        rows,
        key=lambda r: (
            0 if r["protein_id"] in preferred else 1,
            0 if r["label"] == "positive" else 1,
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
    seqs = load_seqs()
    # start from v2
    base = list(csv.DictReader((DATA / "train_formal_v2.csv").open()))
    for r in base:
        r["sequence"] = r["sequence"].strip().upper()
    current = {r["protein_id"]: r for r in base}

    trajectory = []
    # record v2 baseline from saved metrics if present
    v2m = json.loads((DATA / "model_formal.metrics.v2.json").read_text())
    trajectory.append(
        {
            "version": 2,
            "theme": "baseline (prior expansion)",
            "n": v2m["n"],
            "n_pos": v2m["n_pos"],
            "n_neg": v2m["n_neg"],
            "added_attempted": 0,
            "added_kept": 0,
            "skipped_missing_seq": [],
            "dropped_identity": [],
            "mean_auc_logistic": v2m["mean_auc_logistic"],
            "mean_auc_companion": v2m["mean_auc_companion"],
            "oof_auc_logistic": v2m["oof_auc_logistic"],
            "oof_auc_companion": v2m.get("oof_auc_companion"),
            "spearman": v2m.get("spearman_logistic_vs_companion"),
        }
    )

    for ver in range(3, 13):
        theme, batch = ROUNDS[ver]
        attempted = 0
        skipped = []
        newly = []
        for acc, label, kingdom, note in batch:
            attempted += 1
            if acc in current:
                skipped.append(f"{acc}:already_in_set")
                continue
            seq = seqs.get(acc)
            if not seq or len(seq) < 50:
                skipped.append(f"{acc}:missing_seq")
                continue
            newly.append(
                {
                    "protein_id": acc,
                    "label": label,
                    "kingdom": kingdom,
                    "source": f"round_v{ver}",
                    "note": note.replace(",", ";"),
                    "sequence": seq,
                }
            )

        merged = list(current.values()) + newly
        kept, dropped = dedup(merged, 0.90)
        # only keep drops that were in newly or collisions
        current = {r["protein_id"]: r for r in kept}
        rows = list(current.values())
        csv_path = DATA / f"train_formal_v{ver}.csv"
        write_csv(csv_path, rows)

        model = train_from_csv(csv_path, seed=7, k_folds=5, evaluate=True)
        metrics_path = DATA / f"model_formal.metrics.v{ver}.json"
        metrics_path.write_text(json.dumps(model.metrics, indent=2))
        # also update canonical model on last round only later

        n_pos = sum(1 for r in rows if r["label"] == "positive")
        n_neg = sum(1 for r in rows if r["label"] == "negative")
        added_kept = sum(1 for r in newly if r["protein_id"] in current)
        # newly that survived dedup
        added_kept = len([r for r in newly if r["protein_id"] in current])

        entry = {
            "version": ver,
            "theme": theme,
            "n": len(rows),
            "n_pos": n_pos,
            "n_neg": n_neg,
            "added_attempted": attempted,
            "added_raw": len(newly),
            "added_kept": added_kept,
            "skipped_missing_seq": skipped,
            "dropped_identity": dropped,
            "mean_auc_logistic": model.metrics.get("mean_auc_logistic"),
            "mean_auc_companion": model.metrics.get("mean_auc_companion"),
            "oof_auc_logistic": model.metrics.get("oof_auc_logistic"),
            "oof_auc_companion": model.metrics.get("oof_auc_companion"),
            "spearman": model.metrics.get("spearman_logistic_vs_companion"),
        }
        trajectory.append(entry)
        print(
            f"v{ver}: n={entry['n']} (+{entry['added_kept']} kept/{entry['added_raw']} fetched) "
            f"pos={n_pos} neg={n_neg} "
            f"AUC_L={entry['mean_auc_logistic']:.3f} AUC_C={entry['mean_auc_companion']:.3f} "
            f"OOF_L={entry['oof_auc_logistic']:.3f}",
            flush=True,
        )
        if skipped:
            print(f"  skipped: {skipped}", flush=True)
        if dropped:
            print(f"  dedup dropped: {dropped[:8]}{'...' if len(dropped)>8 else ''}", flush=True)

    # finalize canonical files from v12
    v12_csv = DATA / "train_formal_v12.csv"
    (DATA / "train_formal.csv").write_text(v12_csv.read_text())
    model = train_from_csv(v12_csv, seed=7, k_folds=5, evaluate=True)
    model.save(DATA / "model_formal.json")
    (DATA / "model_formal.metrics.json").write_text(json.dumps(model.metrics, indent=2))
    (DATA / "model_formal.metrics.v12.json").write_text(json.dumps(model.metrics, indent=2))

    traj_path = DATA / "train_rounds_v2_to_v12.trajectory.json"
    traj_path.write_text(json.dumps(trajectory, indent=2))
    print(f"Wrote trajectory → {traj_path}")
    print(f"Final n={model.metrics['n']} pos={model.metrics['n_pos']} neg={model.metrics['n_neg']}")


if __name__ == "__main__":
    main()
