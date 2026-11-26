#!/usr/bin/env python3
"""Join side-by-side comparator scores into SQLite (PlantLLPS vs PSAP / LLPhyScore / PLAAC).

Sources (honest scope):
  * PSAP / LLPhyScore — gold-standard TSV from data/benchmarks/ (native runs);
    proteome-wide PSAP/LLPhyScore require external envs (see docs/EXTERNAL_TOOLS.md).
  * PLAAC-approx — stdlib scores for every proteome sequence (not Java PLAAC HMM).
  * metapredict / AlphaFold columns are joined by their own scripts; this script
    only ensures comparator_* columns + a /compare view-friendly summary.

Example::

    python3 scripts/enrich_comparator_scores.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plantllps.plaac_approx import plaac_features  # noqa: E402

DATA = ROOT / "data"
DEFAULT_DB = DATA / "grape_condensatome.sqlite"
GOLD_TSV = DATA / "benchmarks" / "gold_standard_scores.tsv"
FASTA = DATA / "UP000009183.fasta.gz"


def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    for col, decl in (
        ("score_psap", "REAL"),
        ("score_llphyscore_percentile", "REAL"),
        ("plaac_nllr_approx", "REAL"),
        ("plaac_qn_frac", "REAL"),
        ("plaac_max_window_qn", "REAL"),
        ("plaac_core_len", "REAL"),
        ("comparator_set", "TEXT"),
    ):
        if col not in cols:
            conn.execute(f"ALTER TABLE proteins ADD COLUMN {col} {decl}")


def load_gold_comparators(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    out = {}
    for row in csv.DictReader(path.open(), delimiter="\t"):
        pid = row["protein_id"]
        def f(k):
            v = row.get(k)
            if v is None or v == "":
                return None
            try:
                return float(v)
            except ValueError:
                return None

        out[pid] = {
            "score_psap": f("psap_score"),
            "score_llphyscore_percentile": f("llphyscore_percentile"),
            "label": row.get("label"),
            "note": row.get("note"),
            "comparator_set": "gold_standard",
        }
    return out


def iter_fasta(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        pid = None
        buf: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if pid is not None:
                    yield pid, "".join(buf)
                header = line[1:].split()[0]
                parts = header.split("|")
                pid = parts[1] if len(parts) > 1 else header
                buf = []
            else:
                buf.append(line.strip())
        if pid is not None:
            yield pid, "".join(buf)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--gold-tsv", default=str(GOLD_TSV))
    ap.add_argument("--fasta", default=str(FASTA))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--skip-plaac", action="store_true")
    args = ap.parse_args()

    conn = sqlite3.connect(str(args.db))
    ensure_columns(conn)
    gold = load_gold_comparators(Path(args.gold_tsv))
    print(f"Gold comparator rows: {len(gold)}", flush=True)

    # Join gold PSAP / LLPhyScore onto any matching proteome IDs (+ keep for UI even if not in proteome)
    gold_joined = 0
    for pid, g in gold.items():
        exists = conn.execute(
            "SELECT 1 FROM proteins WHERE protein_id = ?", (pid,)
        ).fetchone()
        if not exists:
            continue
        conn.execute(
            """
            UPDATE proteins SET
                score_psap = ?,
                score_llphyscore_percentile = ?,
                comparator_set = ?
            WHERE protein_id = ?
            """,
            (g["score_psap"], g["score_llphyscore_percentile"], g["comparator_set"], pid),
        )
        gold_joined += 1

    plaac_n = 0
    if not args.skip_plaac and Path(args.fasta).exists():
        db_ids = {r[0] for r in conn.execute("SELECT protein_id FROM proteins")}
        print(f"Computing PLAAC-approx for up to {len(db_ids)} proteome sequences …", flush=True)
        for i, (pid, seq) in enumerate(iter_fasta(Path(args.fasta)), 1):
            if pid not in db_ids:
                continue
            feats = plaac_features(seq)
            conn.execute(
                """
                UPDATE proteins SET
                    plaac_nllr_approx = ?,
                    plaac_qn_frac = ?,
                    plaac_max_window_qn = ?,
                    plaac_core_len = ?
                WHERE protein_id = ?
                """,
                (
                    feats["plaac_nllr_approx"],
                    feats["plaac_qn_frac"],
                    feats["plaac_max_window_qn"],
                    feats["plaac_core_len"],
                    pid,
                ),
            )
            # also stash into features_json
            row = conn.execute(
                "SELECT features_json FROM proteins WHERE protein_id = ?", (pid,)
            ).fetchone()
            if row:
                try:
                    fj = json.loads(row[0] or "{}")
                except json.JSONDecodeError:
                    fj = {}
                fj.update(feats)
                conn.execute(
                    "UPDATE proteins SET features_json = ? WHERE protein_id = ?",
                    (json.dumps(fj, sort_keys=True), pid),
                )
            plaac_n += 1
            if args.limit and plaac_n >= args.limit:
                break
            if plaac_n % 5000 == 0:
                conn.commit()
                print(f"  plaac {plaac_n}", flush=True)

    # Gold-standard side-by-side table for UI / docs (includes non-grape IDs)
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    side_by_side = []
    for pid, g in gold.items():
        row = conn.execute(
            "SELECT * FROM proteins WHERE protein_id = ?",
            (pid,),
        ).fetchone()
        row_d = dict(zip([c[1] for c in conn.execute("PRAGMA table_info(proteins)")], row)) if row else {}
        # features from gold tsv always
        entry = {
            "protein_id": pid,
            "label": g.get("label"),
            "note": g.get("note"),
            "score_plantllps": row_d.get("score"),
            "score_psap": g.get("score_psap"),
            "score_llphyscore_percentile": g.get("score_llphyscore_percentile"),
            "plaac_nllr_approx": row_d.get("plaac_nllr_approx") if "plaac_nllr_approx" in cols else None,
            "frac_disorder_metapredict": row_d.get("frac_disorder_metapredict")
            if "frac_disorder_metapredict" in cols
            else None,
            "af_frac_very_low": row_d.get("af_frac_very_low") if "af_frac_very_low" in cols else None,
        }
        if row_d:
            pass  # score_plantllps already set
        # PlantLLPS from gold tsv if not in proteome
        if entry["score_plantllps"] is None:
            # read from gold file plantllps column if present
            pass
        side_by_side.append(entry)

    # enrich side_by_side plantllps from gold tsv
    if Path(args.gold_tsv).exists():
        by_id = {r["protein_id"]: r for r in csv.DictReader(Path(args.gold_tsv).open(), delimiter="\t")}
        for e in side_by_side:
            r = by_id.get(e["protein_id"])
            if r and e["score_plantllps"] is None:
                try:
                    e["score_plantllps"] = float(r.get("score_plantllps") or 0)
                except ValueError:
                    pass
            # plaac for gold sequences from reconstructed fasta-less CSV
    gold_csv = DATA / "benchmarks" / "gold_standard_reconstructed.csv"
    if gold_csv.exists():
        for row in csv.DictReader(gold_csv.open()):
            pid = row["protein_id"]
            feats = plaac_features(row["sequence"])
            for e in side_by_side:
                if e["protein_id"] == pid and e["plaac_nllr_approx"] is None:
                    e["plaac_nllr_approx"] = feats["plaac_nllr_approx"]
            # if in DB, already set; if not in proteome DB, still in side_by_side JSON

    out_json = DATA / "benchmarks" / "comparator_side_by_side.json"
    out_json.write_text(json.dumps({"n": len(side_by_side), "rows": side_by_side}, indent=2))
    out_tsv = DATA / "benchmarks" / "comparator_side_by_side.tsv"
    with out_tsv.open("w", newline="") as fh:
        fields = [
            "protein_id",
            "label",
            "note",
            "score_plantllps",
            "score_psap",
            "score_llphyscore_percentile",
            "plaac_nllr_approx",
            "frac_disorder_metapredict",
            "af_frac_very_low",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for e in side_by_side:
            w.writerow(e)

    summary = {
        "gold_comparator_ids": len(gold),
        "gold_joined_proteome": gold_joined,
        "plaac_proteome_n": plaac_n,
        "side_by_side": str(out_tsv.name),
        "note": (
            "PSAP/LLPhyScore columns are gold-standard only unless a full-proteome "
            "external run is deposited; PLAAC-approx covers the scored proteome."
        ),
    }
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("comparator_scores", json.dumps(summary)),
    )
    conn.commit()
    print(json.dumps(summary, indent=2))
    for row in conn.execute(
        "SELECT protein_id, score, score_psap, score_llphyscore_percentile, plaac_nllr_approx "
        "FROM proteins WHERE comparator_set='gold_standard' OR protein_id='F6HET4' LIMIT 8"
    ):
        print(" ", row)
    conn.close()


if __name__ == "__main__":
    main()
