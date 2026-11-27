#!/usr/bin/env python3
"""Join curated experimentally anchored LLPS gold into grape SQLite + JSON.

Reads data/curated/llps_experimental_gold.tsv, resolves grape orthologs via
explicit column + OrthoDB reverse lookup, and marks matching proteome rows.
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TSV = ROOT / "data" / "curated" / "llps_experimental_gold.tsv"
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
DEFAULT_JSON = ROOT / "data" / "curated" / "llps_experimental_gold.json"


def ensure_columns(con: sqlite3.Connection) -> None:
    cols = {r[1] for r in con.execute("PRAGMA table_info(proteins)")}
    alters = []
    if "llps_gold" not in cols:
        alters.append("ALTER TABLE proteins ADD COLUMN llps_gold INTEGER DEFAULT 0")
    if "llps_gold_json" not in cols:
        alters.append("ALTER TABLE proteins ADD COLUMN llps_gold_json TEXT")
    for sql in alters:
        con.execute(sql)
    con.commit()


def load_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def reverse_ortholog_map(con: sqlite3.Connection) -> dict[str, list[str]]:
    """Map Ath/rice UniProt → grape protein_ids present in atlas."""
    out: dict[str, list[str]] = {}
    for pid, ath, osa in con.execute(
        "SELECT protein_id, ortholog_ath, ortholog_osa FROM proteins"
    ):
        for oid in (ath, osa):
            if not oid:
                continue
            out.setdefault(str(oid).strip().upper(), []).append(pid)
    return out


def load_hooks_grape_ids() -> dict[str, list[str]]:
    """family / Ath ID → grape IDs from cross_crop_hooks.json (only IDs that may exist)."""
    path = ROOT / "data" / "cross_crop_hooks.json"
    if not path.exists():
        return {}
    data = json.loads(path.read_text())
    out: dict[str, list[str]] = {}
    for h in data.get("hooks", []):
        grapes = [g for g in (h.get("grape") or []) if g]
        fam = (h.get("family") or "").strip()
        if fam and grapes:
            out.setdefault(fam.upper(), []).extend(grapes)
        for ath in h.get("arabidopsis") or []:
            out.setdefault(str(ath).upper(), []).extend(grapes)
    return out


# Short / ambiguous symbols must not substring-match (e.g. FUS ⊂ "fusing").
_NAME_MATCH_ALLOW = {
    "MAP65-1",
    "MAP65",
    "ELF3",
    "FLOE1",
    "RBGD2",
    "SEUSS",
    "FCA",
    "FLL2",
    "WVD2",
    "WDL",
    "PHYB",
    "COP1",
    "SPA1",
    "TPX2",
}


def name_symbol_matches(con: sqlite3.Connection, symbol: str) -> list[str]:
    """Conservative grape hits by protein_name / gene_name for allow-listed symbols."""
    sym = (symbol or "").strip()
    if sym.upper() not in {s.upper() for s in _NAME_MATCH_ALLOW} and not any(
        sym.upper().startswith(a.upper()) for a in ("MAP65", "WVD2", "WDL")
    ):
        return []
    if len(sym) < 3:
        return []
    # Prefer token-ish patterns: "ELF3", "ELF3-like", "MAP65-1", not bare 3-letter traps
    like = f"%{sym}%"
    rows = con.execute(
        """SELECT protein_id, protein_name, gene_name FROM proteins
           WHERE protein_name LIKE ? COLLATE NOCASE
              OR gene_name LIKE ? COLLATE NOCASE
           LIMIT 40""",
        (like, like),
    ).fetchall()
    out = []
    su = sym.upper()
    for pid, pname, gname in rows:
        blob = f"{pname or ''} {gname or ''}".upper()
        # reject accidental substring hits unless symbol appears as a token prefix
        if su in blob.replace("-", " ").replace("_", " ").split() or f"{su}-" in blob or f"{su} " in blob or blob.startswith(su):
            out.append(pid)
        elif su.startswith("MAP65") and "MICROTUBULE-ASSOCIATED PROTEIN" in blob and "65" in blob:
            out.append(pid)
        elif su in {"WVD2", "WDL"} and "WVD2" in blob:
            out.append(pid)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", type=Path, default=DEFAULT_TSV)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--json-out", type=Path, default=DEFAULT_JSON)
    args = ap.parse_args()

    rows = load_tsv(args.tsv)
    con = sqlite3.connect(args.db)
    ensure_columns(con)
    con.execute("UPDATE proteins SET llps_gold = 0, llps_gold_json = NULL")
    ortho = reverse_ortholog_map(con)
    hooks = load_hooks_grape_ids()
    present = {
        r[0] for r in con.execute("SELECT protein_id FROM proteins")
    }

    grape_hits: dict[str, list[dict]] = {}
    for r in rows:
        pid = (r.get("protein_id") or "").strip()
        grape_field = (r.get("grape_ortholog") or "").strip()
        symbol = (r.get("gene_symbol") or "").strip()
        targets = set()
        for g in grape_field.replace(";", ",").split(","):
            g = g.strip()
            if g and g in present:
                targets.add(g)
        # Direct grape ID in curated table
        if pid in present:
            targets.add(pid)
        for g in ortho.get(pid.upper(), []):
            if g in present:
                targets.add(g)
        for g in hooks.get(pid.upper(), []):
            if g in present:
                targets.add(g)
        for g in hooks.get(symbol.upper(), []):
            if g in present:
                targets.add(g)
        # Name rescue for allow-listed plant symbols when OrthoDB/hooks miss
        if symbol:
            for g in name_symbol_matches(con, symbol):
                targets.add(g)
        payload = {
            "protein_id": pid,
            "species": r.get("species", ""),
            "gene_symbol": symbol,
            "condensate_type": r.get("condensate_type", ""),
            "assay_evidence": r.get("assay_evidence", ""),
            "citation": r.get("citation", ""),
            "label": r.get("label", ""),
            "notes": r.get("notes", ""),
        }
        for g in targets:
            grape_hits.setdefault(g, []).append(payload)

    n_marked = 0
    for gpid, entries in grape_hits.items():
        # Prefer positive / excluded labels for badge text
        labels = {e["label"] for e in entries}
        flag = 1 if any(x.startswith("pos") for x in labels) else (
            2 if any("excluded" in x for x in labels) else 1
        )
        con.execute(
            "UPDATE proteins SET llps_gold = ?, llps_gold_json = ? WHERE protein_id = ? COLLATE NOCASE",
            (flag, json.dumps(entries), gpid),
        )
        if con.total_changes:
            n_marked += 1
    # total_changes is cumulative; recount
    n_marked = con.execute("SELECT COUNT(*) FROM proteins WHERE llps_gold > 0").fetchone()[0]

    summary = {
        "n_curated_rows": len(rows),
        "n_grape_proteins_marked": n_marked,
        "n_positives": sum(1 for r in rows if str(r.get("label", "")).startswith("pos")),
        "n_excluded": sum(1 for r in rows if "excluded" in str(r.get("label", ""))),
        "source_tsv": str(args.tsv.relative_to(ROOT)),
        "honesty": (
            "Curated overlay from companion gold + modest literature/LLPSDB rows "
            "already in project docs/train; not an exhaustive LLPSDB dump."
        ),
        "rows": rows,
        "grape_join": {k: v for k, v in sorted(grape_hits.items())},
    }
    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (
            "llps_curated_gold",
            json.dumps(
                {
                    "n_curated_rows": summary["n_curated_rows"],
                    "n_grape_proteins_marked": summary["n_grape_proteins_marked"],
                    "n_positives": summary["n_positives"],
                    "n_excluded": summary["n_excluded"],
                    "source_tsv": summary["source_tsv"],
                    "honesty": summary["honesty"],
                }
            ),
        ),
    )
    con.commit()
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(summary, indent=2))
    print(
        f"Curated LLPS gold: {len(rows)} rows → marked {n_marked} grape proteins; "
        f"JSON → {args.json_out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
