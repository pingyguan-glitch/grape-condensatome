#!/usr/bin/env python3
"""Write empty tomato condensatome SQLite schema + blocker JSON when proteome download fails.

Only use after score_second_crop_proteome.py --crop tomato has been attempted and blocked.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data" / "second_crop"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    ap.add_argument(
        "--blocker",
        default="UniProt proteome download blocked or failed after retries",
    )
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    db_path = out / "tomato_condensatome.sqlite"
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(str(db_path))
    conn.executescript(
        """
        CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE proteins (
            protein_id TEXT PRIMARY KEY,
            length INTEGER NOT NULL,
            score REAL NOT NULL,
            score_companion REAL,
            score_logistic REAL,
            tier TEXT NOT NULL,
            percentile REAL,
            gene_name TEXT,
            protein_name TEXT,
            grape_ortholog TEXT,
            grape_ortholog_gene TEXT,
            orthodb_id TEXT
        );
        CREATE INDEX idx_sc_score ON proteins(score DESC);
        CREATE INDEX idx_sc_tier ON proteins(tier);
        """
    )
    meta = {
        "crop": "tomato",
        "species": "Solanum lycopersicum",
        "proteome": "UP000004994",
        "taxon": 4081,
        "n_scored": 0,
        "stub": True,
        "blocker": args.blocker,
        "model": "model_formal.json",
        "disorder": "composition_proxy",
        "honesty": (
            "Empty schema stub only — tomato proteome was not scored. "
            "Re-run scripts/score_second_crop_proteome.py --crop tomato when download works."
        ),
    }
    for k, v in meta.items():
        conn.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?)",
            (k, json.dumps(v) if not isinstance(v, str) else v),
        )
    conn.commit()
    conn.close()

    summary = {
        **meta,
        "tiers": {"A": 0, "B": 0, "C": 0, "D": 0},
        "tsv": None,
        "sqlite": db_path.name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    (out / "tomato_condensatome.summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    join = {
        "crop": "tomato",
        "species": "Solanum lycopersicum",
        "proteome": "UP000004994",
        "n_grape": 0,
        "n_matched": 0,
        "n_matched_with_crop_score": 0,
        "n_hook_overlay": 0,
        "n_crop_scored": 0,
        "crop_db": db_path.name,
        "stub": True,
        "blocker": args.blocker,
        "columns": [
            "ortholog_sly",
            "ortholog_sly_gene",
            "ortholog_sly_score",
            "ortholog_sly_tier",
        ],
    }
    (out / "tomato_grape_ortholog_join.json").write_text(json.dumps(join, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
