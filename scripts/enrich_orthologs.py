#!/usr/bin/env python3
"""Map grapevine proteins to Arabidopsis / rice via shared OrthoDB groups + curated hooks.

1. Stream UniProt OrthoDB xrefs for UP000009183 (grape), UP000006548 (Ath),
   UP000059680 (rice japonica).
2. Join on OrthoDB group id → store primary Ath/rice accessions + gene names.
3. Overlay curated cross-crop condensatome hooks (MAP65 / TPX2 / WVD2 / …).

Example::

    python3 scripts/enrich_orthologs.py
"""

from __future__ import annotations

import argparse
import gzip
import json
import sqlite3
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from _http_util import DEFAULT_UA, urlopen_bytes  # noqa: E402

DATA = ROOT / "data"
DEFAULT_DB = DATA / "grape_condensatome.sqlite"
CACHE = DATA / "orthologs"
HOOKS = DATA / "cross_crop_hooks.json"

PROTEOMES = {
    "grape": ("UP000009183", 29760),
    "ath": ("UP000006548", 3702),
    "osa": ("UP000059680", 39947),
}


def stream_orthodb(proteome: str, dest: Path, force: bool = False) -> Path:
    if dest.exists() and not force and dest.stat().st_size > 1000:
        print(f"Using cache {dest}", flush=True)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    q = urllib.parse.urlencode(
        {
            "compressed": "true",
            "format": "tsv",
            "query": f"(proteome:{proteome})",
            "fields": "accession,gene_names,protein_name,xref_orthodb",
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/stream?{q}"
    print(f"Downloading OrthoDB stream {proteome} → {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": DEFAULT_UA})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=600) as resp, open(dest, "wb") as out:
                while True:
                    chunk = resp.read(1024 * 256)
                    if not chunk:
                        break
                    out.write(chunk)
            print(f"  wrote {dest.stat().st_size} bytes", flush=True)
            return dest
        except Exception as e:
            wait = 4 * (2**attempt)
            print(f"  failed ({e}); retry in {wait}s", flush=True)
            time.sleep(wait)
    raise SystemExit(f"download failed for {proteome}")


def parse_stream(path: Path) -> dict[str, dict]:
    """accession → {gene, protein_name, orthodb: [ids]}"""
    opener = gzip.open if str(path).endswith(".gz") else open
    out: dict[str, dict] = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}

        def col(*names: str) -> int | None:
            for n in names:
                if n in idx:
                    return idx[n]
            return None

        i_acc = col("Entry", "accession")
        i_gene = col("Gene Names", "gene_names")
        i_prot = col("Protein names", "protein_name")
        i_odb = col("OrthoDB", "xref_orthodb")
        if i_acc is None:
            raise SystemExit(f"bad header in {path}: {header}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= i_acc:
                continue
            acc = parts[i_acc].strip()
            gene_raw = parts[i_gene].strip() if i_gene is not None and i_gene < len(parts) else ""
            gene = gene_raw.split()[0] if gene_raw else ""
            pname = parts[i_prot].strip() if i_prot is not None and i_prot < len(parts) else ""
            odb_raw = parts[i_odb].strip() if i_odb is not None and i_odb < len(parts) else ""
            odbs = [x.strip().rstrip(";") for x in odb_raw.replace(";", " ").split() if x.strip()]
            out[acc] = {"gene": gene, "protein_name": pname[:200], "orthodb": odbs}
    return out


def invert_odb(rows: dict[str, dict]) -> dict[str, list[str]]:
    """orthodb id → list of accessions (stable order)."""
    inv: dict[str, list[str]] = defaultdict(list)
    for acc, info in rows.items():
        for oid in info["orthodb"]:
            inv[oid].append(acc)
    return inv


def pick_best(accs: list[str], rows: dict[str, dict], limit: int = 3) -> list[dict]:
    ranked = sorted(
        accs,
        key=lambda a: (
            0 if rows[a]["gene"] else 1,
            0 if rows[a]["gene"].startswith(("At", "Os", "LOC_Os")) else 1,
            a,
        ),
    )
    out = []
    for a in ranked[:limit]:
        out.append(
            {
                "accession": a,
                "gene": rows[a]["gene"],
                "protein_name": rows[a]["protein_name"],
            }
        )
    return out


DEFAULT_HOOKS = {
    "version": 1,
    "note": (
        "Curated cross-crop condensatome anchors (companion MAP/LLPS + literature). "
        "Used to seed UI links even when OrthoDB groups are broad/missing."
    ),
    "hooks": [
        {
            "family": "MAP65-1",
            "grape": ["F6HET4"],
            "arabidopsis": ["Q9FLP0"],
            "rice": ["Q0D4Q4", "Q7X8B4"],
            "other_crops": {
                "soybean": [],
                "maize": [],
                "tomato": [],
            },
            "note": "Proven LLPS driver (Sahu et al. 2023); Ser-rich C-terminal grammar",
        },
        {
            "family": "TPX2-like",
            "grape": [],
            "arabidopsis": ["Q9C5Y0"],
            "rice": [],
            "note": "Animal paradigm extension; often Tier A in grape MAP screen",
        },
        {
            "family": "ELF3",
            "grape": [],
            "arabidopsis": ["O82804"],
            "rice": [],
            "note": "Prion-like thermosensor",
        },
        {
            "family": "FLOE1",
            "grape": [],
            "arabidopsis": ["Q8VZR8"],
            "rice": [],
            "note": "Hydration-dependent condensate",
        },
        {
            "family": "WVD2/WDL",
            "grape": ["A0ACT1LGI9"],
            "arabidopsis": ["Q9LZ18"],
            "rice": [],
            "note": "Plant-specific disordered MAP family; top of extended grape screen",
        },
    ],
}


def ensure_hooks(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(DEFAULT_HOOKS, indent=2))
    return DEFAULT_HOOKS


def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    for col, decl in (
        ("ortholog_ath", "TEXT"),
        ("ortholog_ath_gene", "TEXT"),
        ("ortholog_osa", "TEXT"),
        ("ortholog_osa_gene", "TEXT"),
        ("orthodb_id", "TEXT"),
        ("ortholog_json", "TEXT"),
        ("cross_crop_family", "TEXT"),
    ):
        if col not in cols:
            conn.execute(f"ALTER TABLE proteins ADD COLUMN {col} {decl}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--cache-dir", default=str(CACHE))
    ap.add_argument("--hooks", default=str(HOOKS))
    ap.add_argument("--force-download", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    cache = Path(args.cache_dir)
    hooks = ensure_hooks(Path(args.hooks))

    grape_path = stream_orthodb(
        PROTEOMES["grape"][0], cache / "grape_orthodb.tsv.gz", force=args.force_download
    )
    ath_path = stream_orthodb(
        PROTEOMES["ath"][0], cache / "ath_orthodb.tsv.gz", force=args.force_download
    )
    osa_path = stream_orthodb(
        PROTEOMES["osa"][0], cache / "osa_orthodb.tsv.gz", force=args.force_download
    )

    grape = parse_stream(grape_path)
    ath = parse_stream(ath_path)
    osa = parse_stream(osa_path)
    ath_inv = invert_odb(ath)
    osa_inv = invert_odb(osa)
    print(
        f"Parsed grape={len(grape)} ath={len(ath)} osa={len(osa)} "
        f"ath_odb_groups={len(ath_inv)} osa_odb_groups={len(osa_inv)}",
        flush=True,
    )

    # family hooks → grape id
    hook_by_grape: dict[str, dict] = {}
    for h in hooks.get("hooks") or []:
        for gid in h.get("grape") or []:
            hook_by_grape[gid] = h

    conn = sqlite3.connect(str(args.db))
    ensure_columns(conn)
    ids = [r[0] for r in conn.execute("SELECT protein_id FROM proteins")]
    if args.limit:
        ids = ids[: args.limit]

    matched_ath = matched_osa = matched_hook = 0
    for pid in ids:
        info = grape.get(pid) or {"orthodb": [], "gene": "", "protein_name": ""}
        odbs = info.get("orthodb") or []
        ath_hits: list[dict] = []
        osa_hits: list[dict] = []
        primary_odb = odbs[0] if odbs else None
        for oid in odbs:
            if oid in ath_inv and not ath_hits:
                ath_hits = pick_best(ath_inv[oid], ath)
            if oid in osa_inv and not osa_hits:
                osa_hits = pick_best(osa_inv[oid], osa)
            if ath_hits and osa_hits:
                break

        hook = hook_by_grape.get(pid)
        family = None
        if hook:
            matched_hook += 1
            family = hook.get("family")
            # overlay curated accessions if OrthoDB missed them
            if not ath_hits and hook.get("arabidopsis"):
                ath_hits = [
                    {
                        "accession": a,
                        "gene": (ath.get(a) or {}).get("gene", ""),
                        "protein_name": (ath.get(a) or {}).get("protein_name", ""),
                        "source": "curated_hook",
                    }
                    for a in hook["arabidopsis"]
                ]
            if not osa_hits and hook.get("rice"):
                osa_hits = [
                    {
                        "accession": a,
                        "gene": (osa.get(a) or {}).get("gene", ""),
                        "protein_name": (osa.get(a) or {}).get("protein_name", ""),
                        "source": "curated_hook",
                    }
                    for a in hook["rice"]
                ]

        if ath_hits:
            matched_ath += 1
        if osa_hits:
            matched_osa += 1

        payload = {
            "orthodb": odbs[:5],
            "arabidopsis": ath_hits,
            "rice": osa_hits,
            "cross_crop_family": family,
            "hook_note": (hook or {}).get("note"),
        }
        conn.execute(
            """
            UPDATE proteins SET
                ortholog_ath = ?,
                ortholog_ath_gene = ?,
                ortholog_osa = ?,
                ortholog_osa_gene = ?,
                orthodb_id = ?,
                ortholog_json = ?,
                cross_crop_family = ?
            WHERE protein_id = ?
            """,
            (
                ath_hits[0]["accession"] if ath_hits else None,
                ath_hits[0].get("gene") if ath_hits else None,
                osa_hits[0]["accession"] if osa_hits else None,
                osa_hits[0].get("gene") if osa_hits else None,
                primary_odb,
                json.dumps(payload, sort_keys=True),
                family,
                pid,
            ),
        )

    summary = {
        "n_proteins": len(ids),
        "n_matched_ath": matched_ath,
        "n_matched_osa": matched_osa,
        "n_hook_overlay": matched_hook,
        "proteomes": {k: v[0] for k, v in PROTEOMES.items()},
        "hooks_file": str(Path(args.hooks).name),
    }
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("orthologs", json.dumps(summary)),
    )
    # also store hooks blob for UI
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("cross_crop_hooks", json.dumps(hooks)),
    )
    conn.commit()
    print(json.dumps(summary, indent=2))
    for row in conn.execute(
        "SELECT protein_id, ortholog_ath, ortholog_ath_gene, ortholog_osa, ortholog_osa_gene, cross_crop_family "
        "FROM proteins WHERE protein_id='F6HET4' OR cross_crop_family IS NOT NULL LIMIT 8"
    ):
        print(" ", row)
    conn.close()


if __name__ == "__main__":
    main()
