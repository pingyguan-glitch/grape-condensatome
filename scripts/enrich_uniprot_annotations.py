#!/usr/bin/env python3
"""Download UniProt gene / GO / keyword annotations for UP000009183 and join into SQLite.

Writes annotations JSON into proteins.annotations and adds optional columns
gene_name, protein_name, go_terms, keywords when absent.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
DEFAULT_CACHE = ROOT / "data" / "uniprot_up000009183_annotations.tsv.gz"

FIELDS = "accession,gene_names,protein_name,go_id,keyword,organism_name"
STREAM_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?"
    + urllib.parse.urlencode(
        {
            "compressed": "true",
            "format": "tsv",
            "query": "(proteome:UP000009183)",
            "fields": FIELDS,
        }
    )
)


def download(cache: Path, force: bool = False) -> Path:
    if cache.exists() and not force and cache.stat().st_size > 1000:
        print(f"Using cached annotations → {cache}", flush=True)
        return cache
    cache.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading UniProt annotations…\n  {STREAM_URL}", flush=True)
    req = urllib.request.Request(
        STREAM_URL,
        headers={"User-Agent": "GrapeCondensatome/0.3 (research; pyguan@cau.edu.cn)"},
    )
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=300) as resp, open(cache, "wb") as out:
                while True:
                    chunk = resp.read(1024 * 256)
                    if not chunk:
                        break
                    out.write(chunk)
            print(f"Wrote {cache} ({cache.stat().st_size} bytes)", flush=True)
            return cache
        except (urllib.error.URLError, TimeoutError) as e:
            wait = 4 * (2**attempt)
            print(f"download failed ({e}); retry in {wait}s", flush=True)
            time.sleep(wait)
    raise SystemExit("UniProt download failed after retries")


def parse_tsv_gz(path: Path) -> dict[str, dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    ann: dict[str, dict] = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # UniProt stream headers vary slightly by field names
        # Expected-ish: Entry, Gene Names, Protein names, Gene Ontology IDs, Keywords, Organism
        idx = {h: i for i, h in enumerate(header)}

        def col(*names: str) -> int | None:
            for n in names:
                if n in idx:
                    return idx[n]
            return None

        i_acc = col("Entry", "accession")
        i_gene = col("Gene Names", "gene_names", "Gene Names (primary)")
        i_prot = col("Protein names", "protein_name", "Protein Name")
        i_go = col("Gene Ontology IDs", "go_id", "Gene Ontology (GO)")
        i_kw = col("Keywords", "keyword")
        if i_acc is None:
            raise SystemExit(f"unexpected UniProt TSV header: {header}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= i_acc:
                continue
            acc = parts[i_acc].strip()
            if not acc:
                continue
            gene_raw = parts[i_gene].strip() if i_gene is not None and i_gene < len(parts) else ""
            # Gene Names can be space-separated; keep primary token
            gene_primary = gene_raw.split()[0] if gene_raw else ""
            protein_name = parts[i_prot].strip() if i_prot is not None and i_prot < len(parts) else ""
            go_raw = parts[i_go].strip() if i_go is not None and i_go < len(parts) else ""
            kw_raw = parts[i_kw].strip() if i_kw is not None and i_kw < len(parts) else ""
            go_terms = [g.strip() for g in go_raw.replace(";", " ").split() if g.strip().startswith("GO:")]
            if not go_terms and go_raw:
                go_terms = [g.strip() for g in go_raw.split(";") if g.strip()]
            keywords = [k.strip() for k in kw_raw.split(";") if k.strip()]
            # Truncate long protein names for UI
            if len(protein_name) > 240:
                protein_name = protein_name[:237] + "..."
            ann[acc] = {
                "gene_name": gene_primary,
                "gene_names": gene_raw,
                "protein_name": protein_name,
                "go_terms": go_terms[:40],
                "keywords": keywords[:40],
            }
    return ann


def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    for col, decl in (
        ("gene_name", "TEXT"),
        ("protein_name", "TEXT"),
        ("go_terms", "TEXT"),
        ("keywords", "TEXT"),
    ):
        if col not in cols:
            conn.execute(f"ALTER TABLE proteins ADD COLUMN {col} {decl}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--cache", default=str(DEFAULT_CACHE))
    ap.add_argument("--force-download", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="Debug: annotate at most N rows")
    args = ap.parse_args()

    cache = download(Path(args.cache), force=args.force_download)
    ann = parse_tsv_gz(cache)
    print(f"Parsed {len(ann)} UniProt annotation rows", flush=True)

    db = Path(args.db)
    conn = sqlite3.connect(str(db))
    ensure_columns(conn)
    ids = [r[0] for r in conn.execute("SELECT protein_id FROM proteins")]
    if args.limit:
        ids = ids[: args.limit]
    matched = 0
    updated = 0
    for pid in ids:
        a = ann.get(pid)
        if not a:
            continue
        matched += 1
        payload = json.dumps(a, sort_keys=True)
        conn.execute(
            """
            UPDATE proteins SET
                annotations = ?,
                gene_name = ?,
                protein_name = ?,
                go_terms = ?,
                keywords = ?
            WHERE protein_id = ?
            """,
            (
                payload,
                a.get("gene_name") or None,
                a.get("protein_name") or None,
                ";".join(a.get("go_terms") or []) or None,
                ";".join(a.get("keywords") or []) or None,
                pid,
            ),
        )
        updated += 1
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("annotations", f"uniprot_up000009183 gene/GO/keywords n_matched={matched}"),
    )
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("annotations_cache", str(cache.name)),
    )
    conn.commit()
    n_gene = conn.execute(
        "SELECT COUNT(*) FROM proteins WHERE IFNULL(gene_name,'') != ''"
    ).fetchone()[0]
    print(
        f"DB {db}: matched={matched}/{len(ids)} updated={updated} with_gene_name={n_gene}",
        flush=True,
    )
    # sample
    for row in conn.execute(
        "SELECT protein_id, gene_name, substr(protein_name,1,60) FROM proteins "
        "WHERE gene_name IS NOT NULL LIMIT 5"
    ):
        print(" ", row)
    conn.close()


if __name__ == "__main__":
    main()
