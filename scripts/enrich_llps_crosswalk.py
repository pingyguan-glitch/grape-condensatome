#!/usr/bin/env python3
"""Cross-annotate PlantLLPS gold + grape Tier A/B/C candidates vs external LLPS DBs.

Inputs (from ``scripts/fetch_llps_external_ids.py``):
  data/llps_crosswalk/ids/*.txt

Outputs:
  data/llps_crosswalk/gold_crosswalk.tsv
  data/llps_crosswalk/gold_crosswalk.json
  data/llps_crosswalk/candidates_crosswalk.tsv.gz
  data/llps_crosswalk/summary.json
  SQLite columns on data/grape_condensatome.sqlite:
    ext_llps_dbs, ext_llps_unique, ext_llps_json

Match rules (documented in data/llps_crosswalk/README.md):
  - direct UniProt ID in PhaSePro / LLPSDB / DrLLPS-experimental / CD-CODE
  - Arabidopsis ortholog UniProt in those same sets
  - grape gene_name (VIT_…) in DrLLPS Vitis ortholog-expansion gene list

``ext_llps_unique=1`` means no match under any rule above (unique to
GrapeCondensatome / PlantLLPS relative to these four catalogs).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IDS_DIR = ROOT / "data" / "llps_crosswalk" / "ids"
OUT_DIR = ROOT / "data" / "llps_crosswalk"
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"

DB_KEYS = ("phasepro", "llpsdb", "drllps", "cdcode")


def _load_id_file(path: Path) -> set[str]:
    out: set[str] = set()
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}; run scripts/fetch_llps_external_ids.py first")
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        out.add(line.split()[0])
    return out


def load_external_sets() -> dict[str, set[str]]:
    return {
        "phasepro": _load_id_file(IDS_DIR / "phasepro_uniprot.txt"),
        "llpsdb": _load_id_file(IDS_DIR / "llpsdb_uniprot.txt"),
        "drllps_experimental": _load_id_file(IDS_DIR / "drllps_experimental_uniprot.txt"),
        "cdcode": _load_id_file(IDS_DIR / "cdcode_uniprot.txt"),
        "drllps_vitis_genes": _load_id_file(IDS_DIR / "drllps_vitis_genes.txt"),
    }


def match_uniprot(uid: str | None, sets: dict[str, set[str]]) -> dict[str, str]:
    """Return {db_key: match_type} for a UniProt accession."""
    hits: dict[str, str] = {}
    if not uid:
        return hits
    u = str(uid).strip().split("-")[0]
    if u in sets["phasepro"]:
        hits["phasepro"] = "direct"
    if u in sets["llpsdb"]:
        hits["llpsdb"] = "direct"
    if u in sets["drllps_experimental"]:
        hits["drllps"] = "direct"
    if u in sets["cdcode"]:
        hits["cdcode"] = "direct"
    return hits


def merge_hits(*parts: dict[str, str]) -> dict[str, str]:
    """Prefer direct over ortholog / expansion labels when merging."""
    priority = {"direct": 3, "ortholog_ath": 2, "drllps_vitis_gene": 1}
    out: dict[str, str] = {}
    for part in parts:
        for db, how in part.items():
            if db not in out or priority.get(how, 0) > priority.get(out[db], 0):
                out[db] = how
    return out


def annotate_protein(
    *,
    protein_id: str,
    gene_name: str | None,
    ortholog_ath: str | None,
    sets: dict[str, set[str]],
) -> dict:
    direct = match_uniprot(protein_id, sets)
    via_ath: dict[str, str] = {}
    if ortholog_ath:
        for db, how in match_uniprot(ortholog_ath, sets).items():
            via_ath[db] = "ortholog_ath"
    via_viv: dict[str, str] = {}
    gene = (gene_name or "").strip().split()[0] if gene_name else ""
    if gene and gene in sets["drllps_vitis_genes"]:
        via_viv["drllps"] = "drllps_vitis_gene"

    hits = merge_hits(direct, via_ath, via_viv)
    dbs = sorted(hits.keys())
    return {
        "protein_id": protein_id,
        "ext_llps_dbs": ";".join(dbs),
        "ext_llps_unique": 0 if dbs else 1,
        "n_ext_dbs": len(dbs),
        "match": hits,
        "match_gene": gene or None,
        "match_ortholog_ath": (ortholog_ath or None),
    }


def ensure_columns(con: sqlite3.Connection) -> None:
    cols = {r[1] for r in con.execute("PRAGMA table_info(proteins)")}
    alters = []
    if "ext_llps_dbs" not in cols:
        alters.append("ALTER TABLE proteins ADD COLUMN ext_llps_dbs TEXT")
    if "ext_llps_unique" not in cols:
        alters.append("ALTER TABLE proteins ADD COLUMN ext_llps_unique INTEGER DEFAULT 1")
    if "ext_llps_json" not in cols:
        alters.append("ALTER TABLE proteins ADD COLUMN ext_llps_json TEXT")
    for sql in alters:
        con.execute(sql)
    con.commit()


def load_gold_panel() -> list[dict]:
    """PlantLLPS gold-standard panel: train ∪ curated ∪ reconstructed gold."""
    rows: dict[str, dict] = {}

    def upsert(pid: str, **kwargs: object) -> None:
        pid = pid.strip()
        if not pid:
            return
        cur = rows.setdefault(
            pid,
            {
                "protein_id": pid,
                "in_train": "no",
                "train_label": "",
                "in_curated_gold": "no",
                "in_reconstructed_gold": "no",
                "species": "",
                "gene_symbol": "",
                "grape_ortholog": "",
                "sources": [],
            },
        )
        for k, v in kwargs.items():
            if k == "sources" and isinstance(v, str):
                if v not in cur["sources"]:
                    cur["sources"].append(v)
            elif k in ("in_train", "in_curated_gold", "in_reconstructed_gold") and v == "yes":
                cur[k] = "yes"
            elif k in ("species", "gene_symbol", "grape_ortholog", "train_label") and v:
                if not cur.get(k):
                    cur[k] = v

    train_path = ROOT / "data" / "train_formal_v13.csv"
    with train_path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            upsert(
                r["protein_id"],
                in_train="yes",
                train_label=r.get("label") or "",
                sources="train_formal_v13",
            )

    curated_path = ROOT / "data" / "curated" / "llps_experimental_gold.tsv"
    with curated_path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            upsert(
                r["protein_id"],
                in_curated_gold="yes",
                species=r.get("species") or "",
                gene_symbol=r.get("gene_symbol") or "",
                grape_ortholog=r.get("grape_ortholog") or "",
                sources="curated_llps_experimental_gold",
            )
            if r.get("grape_ortholog"):
                upsert(
                    r["grape_ortholog"],
                    in_curated_gold="yes",
                    species="Vitis vinifera",
                    gene_symbol=r.get("gene_symbol") or "",
                    grape_ortholog=r["grape_ortholog"],
                    sources="curated_grape_ortholog",
                )

    gold_path = ROOT / "data" / "benchmarks" / "gold_standard_reconstructed.csv"
    with gold_path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            upsert(
                r["protein_id"],
                in_reconstructed_gold="yes",
                train_label=r.get("label") or "",
                sources="gold_standard_reconstructed",
            )

    out = []
    for pid, row in sorted(rows.items()):
        row["sources"] = ";".join(row["sources"])
        out.append(row)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument(
        "--tiers",
        default="A,B,C",
        help="Browsable candidate tiers to deposit (default A,B,C)",
    )
    args = ap.parse_args()

    sets = load_external_sets()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- Gold panel ---
    gold_rows = load_gold_panel()
    gold_out = []
    gold_unique = gold_shared = 0
    for g in gold_rows:
        # For non-grape gold, ortholog_ath is N/A; use grape_ortholog as secondary ID check
        hits = merge_hits(
            match_uniprot(g["protein_id"], sets),
            {
                db: "grape_ortholog"
                for db in match_uniprot(g.get("grape_ortholog") or "", sets)
            },
        )
        # Also check DrLLPS Vitis if this is a grape ID with gene later via DB
        dbs = sorted(hits.keys())
        unique = 0 if dbs else 1
        if unique:
            gold_unique += 1
        else:
            gold_shared += 1
        gold_out.append(
            {
                **g,
                "ext_llps_dbs": ";".join(dbs),
                "ext_llps_unique": unique,
                "n_ext_dbs": len(dbs),
                "in_phasepro": 1 if "phasepro" in hits else 0,
                "in_llpsdb": 1 if "llpsdb" in hits else 0,
                "in_drllps": 1 if "drllps" in hits else 0,
                "in_cdcode": 1 if "cdcode" in hits else 0,
                "match_json": json.dumps(hits, sort_keys=True),
            }
        )

    gold_tsv = OUT_DIR / "gold_crosswalk.tsv"
    fields = [
        "protein_id",
        "species",
        "gene_symbol",
        "grape_ortholog",
        "in_train",
        "train_label",
        "in_curated_gold",
        "in_reconstructed_gold",
        "sources",
        "ext_llps_dbs",
        "ext_llps_unique",
        "n_ext_dbs",
        "in_phasepro",
        "in_llpsdb",
        "in_drllps",
        "in_cdcode",
        "match_json",
    ]
    with gold_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in gold_out:
            w.writerow(row)

    gold_json = {
        "n_panel": len(gold_out),
        "n_unique": gold_unique,
        "n_shared": gold_shared,
        "note": (
            "PlantLLPS gold-standard panel = train_formal_v13 ∪ curated experimental gold "
            "∪ gold_standard_reconstructed (unique UniProt IDs)."
        ),
        "rows": gold_out,
    }
    (OUT_DIR / "gold_crosswalk.json").write_text(
        json.dumps(gold_json, indent=2) + "\n", encoding="utf-8"
    )

    # --- Proteome SQLite + Tier A/B/C candidates ---
    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row
    ensure_columns(con)

    tier_keep = {t.strip().upper() for t in str(args.tiers).split(",") if t.strip()}
    cand_path = OUT_DIR / "candidates_crosswalk.tsv.gz"
    cand_fields = [
        "protein_id",
        "gene_name",
        "tier",
        "score",
        "ortholog_ath",
        "ext_llps_dbs",
        "ext_llps_unique",
        "n_ext_dbs",
        "in_phasepro",
        "in_llpsdb",
        "in_drllps",
        "in_cdcode",
        "match_json",
    ]

    summary_tiers: dict[str, dict] = {
        t: {"n": 0, "unique": 0, "shared": 0, "by_db": Counter()} for t in "ABCD"
    }
    n_updated = 0
    cand_rows = []

    # Reset then update all proteins (D also gets flags for protein pages)
    con.execute("UPDATE proteins SET ext_llps_dbs = NULL, ext_llps_unique = 1, ext_llps_json = NULL")
    for row in con.execute(
        "SELECT protein_id, gene_name, ortholog_ath, tier, score FROM proteins"
    ):
        ann = annotate_protein(
            protein_id=row["protein_id"],
            gene_name=row["gene_name"],
            ortholog_ath=row["ortholog_ath"],
            sets=sets,
        )
        payload = {
            "dbs": ann["match"],
            "gene": ann["match_gene"],
            "ortholog_ath": ann["match_ortholog_ath"],
        }
        con.execute(
            "UPDATE proteins SET ext_llps_dbs = ?, ext_llps_unique = ?, ext_llps_json = ? "
            "WHERE protein_id = ?",
            (
                ann["ext_llps_dbs"] or None,
                ann["ext_llps_unique"],
                json.dumps(payload, sort_keys=True),
                row["protein_id"],
            ),
        )
        n_updated += 1
        tier = (row["tier"] or "?").upper()
        if tier in summary_tiers:
            summary_tiers[tier]["n"] += 1
            if ann["ext_llps_unique"]:
                summary_tiers[tier]["unique"] += 1
            else:
                summary_tiers[tier]["shared"] += 1
                for db in ann["match"]:
                    summary_tiers[tier]["by_db"][db] += 1

        if tier in tier_keep:
            cand_rows.append(
                {
                    "protein_id": row["protein_id"],
                    "gene_name": row["gene_name"] or "",
                    "tier": tier,
                    "score": row["score"],
                    "ortholog_ath": row["ortholog_ath"] or "",
                    "ext_llps_dbs": ann["ext_llps_dbs"],
                    "ext_llps_unique": ann["ext_llps_unique"],
                    "n_ext_dbs": ann["n_ext_dbs"],
                    "in_phasepro": 1 if "phasepro" in ann["match"] else 0,
                    "in_llpsdb": 1 if "llpsdb" in ann["match"] else 0,
                    "in_drllps": 1 if "drllps" in ann["match"] else 0,
                    "in_cdcode": 1 if "cdcode" in ann["match"] else 0,
                    "match_json": json.dumps(ann["match"], sort_keys=True),
                }
            )

    con.commit()

    # Sort candidates by score desc
    cand_rows.sort(key=lambda r: (-float(r["score"] or 0), r["protein_id"]))
    with gzip.open(cand_path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cand_fields, delimiter="\t")
        w.writeheader()
        for r in cand_rows:
            w.writerow(r)

    def pack_tier(t: str) -> dict:
        s = summary_tiers[t]
        return {
            "n": s["n"],
            "unique": s["unique"],
            "shared": s["shared"],
            "by_db": dict(s["by_db"]),
        }

    abc_unique = sum(summary_tiers[t]["unique"] for t in "ABC")
    abc_shared = sum(summary_tiers[t]["shared"] for t in "ABC")
    ab_unique = sum(summary_tiers[t]["unique"] for t in "AB")
    ab_shared = sum(summary_tiers[t]["shared"] for t in "AB")

    summary = {
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "tool_version": (ROOT / "VERSION").read_text(encoding="utf-8").strip(),
        "external_id_counts": {
            "phasepro": len(sets["phasepro"]),
            "llpsdb": len(sets["llpsdb"]),
            "drllps_experimental": len(sets["drllps_experimental"]),
            "cdcode": len(sets["cdcode"]),
            "drllps_vitis_genes": len(sets["drllps_vitis_genes"]),
        },
        "gold_panel": {
            "n": len(gold_out),
            "unique": gold_unique,
            "shared": gold_shared,
            "definition": "train_formal_v13 ∪ curated_llps_experimental_gold ∪ gold_standard_reconstructed",
        },
        "proteome": {
            "n_updated": n_updated,
            "tiers": {t: pack_tier(t) for t in "ABCD"},
            "tier_ABC": {"n": abc_unique + abc_shared, "unique": abc_unique, "shared": abc_shared},
            "tier_AB": {"n": ab_unique + ab_shared, "unique": ab_unique, "shared": ab_shared},
        },
        "match_rules": [
            "direct UniProt in PhaSePro / LLPSDB / DrLLPS-experimental / CD-CODE",
            "ortholog_ath UniProt in those sets",
            "gene_name in DrLLPS Vitis ortholog-expansion gene list → counts as DrLLPS",
        ],
        "honesty": (
            "DrLLPS Vitis matches are ortholog expansions of animal/yeast scaffolds, not "
            "plant-experimental LLPS assays. 'Unique' means absent from these four catalogs "
            "under the match rules — not experimental proof of novelty."
        ),
        "files": {
            "gold_tsv": "data/llps_crosswalk/gold_crosswalk.tsv",
            "gold_json": "data/llps_crosswalk/gold_crosswalk.json",
            "candidates_tsv_gz": "data/llps_crosswalk/candidates_crosswalk.tsv.gz",
            "ids": "data/llps_crosswalk/ids/",
        },
    }

    # Serialize Counters
    for t in summary["proteome"]["tiers"].values():
        t["by_db"] = dict(t["by_db"])

    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("llps_crosswalk", json.dumps(summary)),
    )
    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("tool_version", (ROOT / "VERSION").read_text(encoding="utf-8").strip()),
    )
    con.commit()
    con.close()

    print(json.dumps({
        "gold": summary["gold_panel"],
        "tier_AB": summary["proteome"]["tier_AB"],
        "tier_ABC": summary["proteome"]["tier_ABC"],
        "tiers": summary["proteome"]["tiers"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
