"""Command-line interface: train, score FASTA → TSV, proteome → SQLite."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import List, Optional

from . import __version__
from .db import connect, set_meta, summarize, upsert_protein
from .fasta import parse_fasta
from .features import FEATURE_NAMES
from .model import (
    PlantLLPSModel,
    assign_tier,
    ensure_trained_model,
    percentile_cutoffs,
    train_from_csv,
)

PACKAGE_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TRAIN = PACKAGE_ROOT / "data" / "train_formal.csv"
DEFAULT_TRAIN_FALLBACK = PACKAGE_ROOT / "data" / "train_tiny.csv"
DEFAULT_SEED = PACKAGE_ROOT / "data" / "seed_grapevine.fasta"
DEFAULT_DB = PACKAGE_ROOT / "data" / "grape_condensatome.sqlite"
DEFAULT_MODEL = PACKAGE_ROOT / "data" / "model_formal.json"
DEFAULT_PROTEOME = PACKAGE_ROOT / "data" / "UP000009183.fasta.gz"

TSV_COLUMNS: List[str] = (
    ["protein_id", "score_plantllps", "score_companion", "score_logistic", "tier", "percentile"]
    + FEATURE_NAMES
)


def _resolve_train_csv(path: Optional[str]) -> Path:
    if path:
        return Path(path)
    if DEFAULT_TRAIN.exists():
        return DEFAULT_TRAIN
    return DEFAULT_TRAIN_FALLBACK


def cmd_train(args: argparse.Namespace) -> int:
    csv_path = _resolve_train_csv(args.train_csv)
    out = Path(args.model_out)
    model = train_from_csv(csv_path, seed=args.seed, k_folds=args.folds, evaluate=True)
    model.save(out)
    metrics_path = out.with_suffix(".metrics.json")
    metrics_path.write_text(json.dumps(model.metrics, indent=2))
    print(f"Wrote model → {out}", file=sys.stderr)
    print(model.train_note, file=sys.stderr)
    if model.metrics:
        m = model.metrics
        print(
            f"Holdout (stratified {m.get('k')}-fold): "
            f"AUC_logistic={m.get('mean_auc_logistic'):.3f} "
            f"AUC_companion={m.get('mean_auc_companion'):.3f} "
            f"OOF_AUC_logistic={m.get('oof_auc_logistic'):.3f} "
            f"ρ(logi,comp)={m.get('spearman_logistic_vs_companion'):.3f} "
            f"n={m.get('n')} (pos={m.get('n_pos')}, neg={m.get('n_neg')})",
            file=sys.stderr,
        )
        print(f"Metrics JSON → {metrics_path}", file=sys.stderr)
    print(
        "Caveat: plant positives remain scarce; AUC is on folded-heavy negatives "
        "and composition-proxy disorder — not a final claim vs PSAP/LLPhyScore.",
        file=sys.stderr,
    )
    return 0


def cmd_score(args: argparse.Namespace) -> int:
    fasta = Path(args.fasta)
    out_tsv = Path(args.out_tsv)
    model = ensure_trained_model(PACKAGE_ROOT)
    if args.model:
        model = PlantLLPSModel.load(args.model)
    if args.mode:
        model.mode = args.mode

    records = parse_fasta(fasta)
    if not records:
        print(f"No sequences in {fasta}", file=sys.stderr)
        return 1

    min_len = int(args.min_length)
    records = [(p, s) for p, s in records if len(s) >= min_len]
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for pid, seq in records:
        scored = model.score(seq)
        row = {"protein_id": pid, "percentile": "", **scored}
        rows.append(row)

    # percentiles within this batch (rank-based)
    scores = [r["score_plantllps"] for r in rows]
    order = sorted(range(len(rows)), key=lambda i: scores[i])
    ranks = [0.0] * len(rows)
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and scores[order[j]] == scores[order[i]]:
            j += 1
        avg = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = avg
        i = j
    for idx, r in enumerate(rows):
        r["percentile"] = round(100.0 * ranks[idx] / len(rows), 3)

    with open(out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in TSV_COLUMNS})

    print(f"Scored {len(rows)} sequences → {out_tsv}", file=sys.stderr)

    if args.db:
        _write_db(Path(args.db), records, rows, source=str(fasta.name), status=args.db_status)
    return 0


def _write_db(db_path: Path, records, rows, source: str, status: str) -> None:
    conn = connect(db_path)
    set_meta(conn, "tool_version", __version__)
    set_meta(conn, "status", status)
    set_meta(
        conn,
        "disclaimer",
        "Scores use composition-proxy disorder; not metapredict/pLDDT. "
        "Training set is LLPSDB plant (filtered) + companion gold standards; still small.",
    )
    for (pid, seq), row in zip(records, rows):
        upsert_protein(
            conn,
            protein_id=pid,
            sequence=seq,
            score_row=row,
            annotations=None,
            source=source,
        )
    conn.commit()
    summary = summarize(conn)
    conn.close()
    print(f"SQLite → {db_path}", file=sys.stderr)
    print(f"  proteins={summary['n_proteins']} tiers={summary['tiers']}", file=sys.stderr)


def cmd_init_db(args: argparse.Namespace) -> int:
    args.fasta = args.fasta or str(DEFAULT_SEED)
    args.out_tsv = args.out_tsv or str(PACKAGE_ROOT / "data" / "seed_scores.tsv")
    args.db = args.db or str(PACKAGE_ROOT / "data" / "grape_condensatome.prototype.sqlite")
    args.model = getattr(args, "model", None)
    args.mode = getattr(args, "mode", None)
    args.min_length = getattr(args, "min_length", 80)
    args.db_status = "prototype_seed_only"
    if not DEFAULT_MODEL.exists():
        train_csv = _resolve_train_csv(None)
        if train_csv.exists():
            m = train_from_csv(train_csv)
            m.save(DEFAULT_MODEL)
    return cmd_score(args)


def cmd_score_proteome(args: argparse.Namespace) -> int:
    """Score Vitis reference proteome and write TSV + SQLite with calibrated tiers."""
    fasta = Path(args.fasta)
    if not fasta.exists():
        # try gz sibling
        gz = Path(str(fasta) + ".gz") if not str(fasta).endswith(".gz") else fasta
        alt = PACKAGE_ROOT / "data" / "UP000009183.fasta.gz"
        if gz.exists():
            fasta = gz
        elif alt.exists():
            fasta = alt
        else:
            print(f"Proteome FASTA not found: {args.fasta}", file=sys.stderr)
            return 1

    # Decompress if needed to a temp-ish path under data/
    if str(fasta).endswith(".gz"):
        import gzip

        plain = PACKAGE_ROOT / "data" / "UP000009183.fasta"
        if not plain.exists() or plain.stat().st_mtime < fasta.stat().st_mtime:
            print(f"Decompressing {fasta} → {plain}", file=sys.stderr)
            with gzip.open(fasta, "rt") as src, open(plain, "w") as dst:
                while True:
                    chunk = src.read(1 << 20)
                    if not chunk:
                        break
                    dst.write(chunk)
        fasta = plain

    model = PlantLLPSModel.load(args.model) if args.model else ensure_trained_model(PACKAGE_ROOT)
    if args.mode:
        model.mode = args.mode

    print(f"Parsing {fasta} …", file=sys.stderr)
    records = parse_fasta(fasta)
    min_len = int(args.min_length)
    records = [(p, s) for p, s in records if len(s) >= min_len]
    # drop sequences with too many nonstandard
    clean = []
    for p, s in records:
        aa = sum(1 for c in s.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
        if aa / max(1, len(s)) >= 0.95:
            clean.append((p, s))
    records = clean
    print(f"Scoring {len(records)} proteins (min_length={min_len}) …", file=sys.stderr)

    t0 = time.time()
    rows = []
    for i, (pid, seq) in enumerate(records, 1):
        scored = model.score(seq)
        rows.append({"protein_id": pid, "percentile": "", **scored})
        if i % 2000 == 0:
            elapsed = time.time() - t0
            rate = i / max(elapsed, 1e-6)
            eta = (len(records) - i) / max(rate, 1e-6)
            print(f"  {i}/{len(records)} ({rate:.0f}/s, ETA {eta/60:.1f} min)", file=sys.stderr)

    # Calibrate tiers on proteome percentiles (O(n log n) rank → percentile)
    scores = [r["score_plantllps"] for r in rows]
    cutoffs = percentile_cutoffs(scores, pct_a=float(args.tier_a_pct), pct_c=float(args.tier_c_pct))
    model.tier_cutoffs = cutoffs
    order = sorted(range(len(rows)), key=lambda i: scores[i])
    ranks = [0.0] * len(rows)
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and scores[order[j]] == scores[order[i]]:
            j += 1
        # average 1-based rank
        avg = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = avg
        i = j
    n = len(rows)
    for idx, r in enumerate(rows):
        feats = {k: r[k] for k in FEATURE_NAMES if k in r}
        r["tier"] = assign_tier(r["score_plantllps"], feats, cutoffs)
        r["percentile"] = round(100.0 * ranks[idx] / n, 3)

    # persist cutoffs on model
    if args.model_out:
        model.save(args.model_out)
    else:
        model.save(DEFAULT_MODEL)

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in TSV_COLUMNS})

    db_path = Path(args.db)
    conn = connect(db_path)
    set_meta(conn, "tool_version", __version__)
    set_meta(conn, "status", "proteome_scored_v1")
    set_meta(conn, "proteome", "UniProt UP000009183 Vitis vinifera")
    set_meta(conn, "n_scored", str(len(rows)))
    set_meta(conn, "tier_cutoffs", json.dumps(cutoffs))
    set_meta(
        conn,
        "disclaimer",
        "Whole-proteome prototype scores with composition-proxy disorder. "
        "Not a finished atlas with Expression Atlas / metapredict.",
    )
    for (pid, seq), row in zip(records, rows):
        upsert_protein(
            conn,
            protein_id=pid,
            sequence=seq,
            score_row=row,
            annotations=None,
            source="UP000009183",
        )
    conn.commit()
    tier_counts = Counter(r["tier"] for r in rows)
    top = sorted(rows, key=lambda r: -r["score_plantllps"])[:20]
    summary_path = Path(args.summary_json)
    summary = {
        "n_scored": len(rows),
        "min_length": min_len,
        "tier_counts": dict(tier_counts),
        "tier_cutoffs": cutoffs,
        "elapsed_sec": round(time.time() - t0, 1),
        "top20": [
            {
                "protein_id": r["protein_id"],
                "score": round(r["score_plantllps"], 4),
                "tier": r["tier"],
                "percentile": r["percentile"],
                "disorder_proxy": round(r.get("frac_disorder_proxy", 0), 4),
                "map65_tail_like": r.get("map65_tail_like", 0),
            }
            for r in top
        ],
        "tsv": str(out_tsv),
        "sqlite": str(db_path),
    }
    summary_path.write_text(json.dumps(summary, indent=2))
    conn.close()

    print(f"TSV → {out_tsv}", file=sys.stderr)
    print(f"SQLite → {db_path}", file=sys.stderr)
    print(f"Summary → {summary_path}", file=sys.stderr)
    print(f"Tier counts: {dict(tier_counts)}", file=sys.stderr)
    print(f"Cutoffs: {cutoffs}", file=sys.stderr)
    print("Top 5:", file=sys.stderr)
    for r in top[:5]:
        print(
            f"  {r['protein_id']}\t{r['score_plantllps']:.3f}\ttier={r['tier']}\tpct={r['percentile']}",
            file=sys.stderr,
        )
    return 0


def cmd_serve(args: argparse.Namespace) -> int:
    from .web import run_server

    run_server(
        host=args.host,
        port=int(args.port),
        db_path=Path(args.db),
        model_path=Path(args.model) if args.model else None,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plantllps",
        description="PlantLLPS: plant-aware LLPS features + GrapeCondensatome SQLite",
    )
    p.add_argument("--version", action="version", version=f"plantllps {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    t = sub.add_parser("train", help="Train logistic model + report stratified holdout AUC")
    t.add_argument("--train-csv", default=None, help="Default: data/train_formal.csv")
    t.add_argument("--model-out", default=str(DEFAULT_MODEL))
    t.add_argument("--seed", type=int, default=7)
    t.add_argument("--folds", type=int, default=5)
    t.set_defaults(func=cmd_train)

    s = sub.add_parser("score", help="Score FASTA → TSV (+ optional SQLite)")
    s.add_argument("fasta", help="Input FASTA")
    s.add_argument("-o", "--out-tsv", required=True, help="Output TSV path")
    s.add_argument("--db", help="Optional SQLite path to upsert scores")
    s.add_argument("--model", help="Optional model JSON path")
    s.add_argument("--mode", choices=["hybrid", "companion_composite", "logistic_stub"])
    s.add_argument("--min-length", type=int, default=1)
    s.add_argument("--db-status", default="custom_score")
    s.set_defaults(func=cmd_score)

    i = sub.add_parser("init-db", help="Score seed grapevine FASTA → prototype SQLite")
    i.add_argument("--fasta", default=str(DEFAULT_SEED))
    i.add_argument("--out-tsv", default=str(PACKAGE_ROOT / "data" / "seed_scores.tsv"))
    i.add_argument("--db", default=str(PACKAGE_ROOT / "data" / "grape_condensatome.prototype.sqlite"))
    i.add_argument("--model", default=None)
    i.add_argument("--mode", default=None)
    i.add_argument("--min-length", type=int, default=80)
    i.set_defaults(func=cmd_init_db)

    sp = sub.add_parser("score-proteome", help="Score UniProt Vitis proteome → TSV + SQLite")
    sp.add_argument(
        "--fasta",
        default=str(DEFAULT_PROTEOME),
        help="Proteome FASTA or .fasta.gz (default: data/UP000009183.fasta.gz)",
    )
    sp.add_argument(
        "-o",
        "--out-tsv",
        default=str(PACKAGE_ROOT / "data" / "grape_condensatome.proteome.tsv"),
    )
    sp.add_argument("--db", default=str(DEFAULT_DB))
    sp.add_argument("--model", default=None)
    sp.add_argument("--model-out", default=str(DEFAULT_MODEL))
    sp.add_argument("--mode", default=None)
    sp.add_argument("--min-length", type=int, default=80)
    sp.add_argument("--tier-a-pct", type=float, default=95.0, help="Percentile floor for Tier A")
    sp.add_argument("--tier-c-pct", type=float, default=50.0, help="Percentile floor for Tier C")
    sp.add_argument(
        "--summary-json",
        default=str(PACKAGE_ROOT / "data" / "proteome_summary.json"),
    )
    sp.set_defaults(func=cmd_score_proteome)

    sv = sub.add_parser(
        "serve",
        help="Local GrapeCondensatome web UI (stdlib; not public hosting)",
    )
    sv.add_argument("--host", default="0.0.0.0")
    sv.add_argument("--port", type=int, default=8765)
    sv.add_argument("--db", default=str(DEFAULT_DB))
    sv.add_argument("--model", default=str(DEFAULT_MODEL))
    sv.set_defaults(func=cmd_serve)

    return p


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
