#!/usr/bin/env python3
"""Join a scored second-crop condensatome to grape via shared OrthoDB groups.

1. Stream UniProt OrthoDB xrefs for grape + second crop (default soybean).
2. Join OrthoDB groups → grape proteins get ortholog_gma / ortholog_sly columns
   + second-crop PlantLLPS score/tier when the crop SQLite is present.
3. Update crop SQLite with reverse grape ortholog pointers.
4. Overlay curated cross_crop_hooks.json other_crops IDs when present.

Example::

    python3 scripts/enrich_second_crop_orthologs.py
    python3 scripts/enrich_second_crop_orthologs.py --crop tomato
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
from _http_util import DEFAULT_UA  # noqa: E402

DATA = ROOT / "data"
DEFAULT_DB = DATA / "grape_condensatome.sqlite"
CACHE = DATA / "orthologs"
HOOKS = DATA / "cross_crop_hooks.json"
SECOND = DATA / "second_crop"

CROPS = {
    "soybean": {
        "proteome": "UP000008827",
        "taxon": 3847,
        "species": "Glycine max",
        "slug": "soybean",
        "col": "gma",
        "label": "Soybean",
        "gene_prefixes": ("GLYMA", "Glyma", "LOC"),
    },
    "tomato": {
        "proteome": "UP000004994",
        "taxon": 4081,
        "species": "Solanum lycopersicum",
        "slug": "tomato",
        "col": "sly",
        "label": "Tomato",
        "gene_prefixes": ("Solyc", "LOC"),
    },
}

GRAPE_PROTEOME = "UP000009183"


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
    inv: dict[str, list[str]] = defaultdict(list)
    for acc, info in rows.items():
        for oid in info["orthodb"]:
            inv[oid].append(acc)
    return inv


def pick_best(accs: list[str], rows: dict[str, dict], prefixes: tuple[str, ...], limit: int = 3) -> list[dict]:
    ranked = sorted(
        accs,
        key=lambda a: (
            0 if rows[a]["gene"] else 1,
            0 if any(rows[a]["gene"].startswith(p) for p in prefixes) else 1,
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


def ensure_grape_columns(conn: sqlite3.Connection, col: str) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    for name, decl in (
        (f"ortholog_{col}", "TEXT"),
        (f"ortholog_{col}_gene", "TEXT"),
        (f"ortholog_{col}_score", "REAL"),
        (f"ortholog_{col}_tier", "TEXT"),
    ):
        if name not in cols:
            conn.execute(f"ALTER TABLE proteins ADD COLUMN {name} {decl}")


def load_crop_scores(db_path: Path) -> dict[str, dict]:
    if not db_path.exists():
        return {}
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    out = {}
    for r in conn.execute(
        "SELECT protein_id, score, tier, length, percentile FROM proteins"
    ):
        out[r["protein_id"]] = dict(r)
    conn.close()
    return out


def ensure_hooks(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {"version": 1, "hooks": []}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--crop", choices=sorted(CROPS), default="soybean")
    ap.add_argument("--db", default=str(DEFAULT_DB), help="Grape condensatome SQLite")
    ap.add_argument("--crop-db", default="", help="Second-crop slim SQLite (default under data/second_crop/)")
    ap.add_argument("--cache-dir", default=str(CACHE))
    ap.add_argument("--hooks", default=str(HOOKS))
    ap.add_argument("--force-download", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    crop = CROPS[args.crop]
    col = crop["col"]
    cache = Path(args.cache_dir)
    hooks = ensure_hooks(Path(args.hooks))
    crop_db_path = Path(args.crop_db) if args.crop_db else SECOND / f"{crop['slug']}_condensatome.sqlite"
    crop_scores = load_crop_scores(crop_db_path)
    print(f"Loaded {len(crop_scores)} scored {crop['slug']} proteins from {crop_db_path}", flush=True)

    grape_path = stream_orthodb(
        GRAPE_PROTEOME, cache / "grape_orthodb.tsv.gz", force=args.force_download
    )
    crop_odb_path = stream_orthodb(
        crop["proteome"],
        cache / f"{crop['slug']}_orthodb.tsv.gz",
        force=args.force_download,
    )
    grape = parse_stream(grape_path)
    crop_rows = parse_stream(crop_odb_path)
    crop_inv = invert_odb(crop_rows)
    print(
        f"Parsed grape={len(grape)} {crop['slug']}={len(crop_rows)} "
        f"odb_groups={len(crop_inv)}",
        flush=True,
    )

    # curated other_crops overlays keyed by grape id
    hook_by_grape: dict[str, dict] = {}
    for h in hooks.get("hooks") or []:
        for gid in h.get("grape") or []:
            hook_by_grape[gid] = h

    conn = sqlite3.connect(str(args.db))
    ensure_grape_columns(conn, col)
    ids = [r[0] for r in conn.execute("SELECT protein_id FROM proteins")]
    if args.limit:
        ids = ids[: args.limit]

    matched = matched_scored = matched_hook = 0
    grape_to_crop: dict[str, list[dict]] = {}

    for pid in ids:
        info = grape.get(pid) or {"orthodb": [], "gene": "", "protein_name": ""}
        odbs = info.get("orthodb") or []
        hits: list[dict] = []
        primary_odb = odbs[0] if odbs else None
        for oid in odbs:
            if oid in crop_inv:
                hits = pick_best(crop_inv[oid], crop_rows, crop["gene_prefixes"])
                primary_odb = oid
                break

        hook = hook_by_grape.get(pid)
        other = ((hook or {}).get("other_crops") or {}).get(crop["slug"]) or []
        if not hits and other:
            matched_hook += 1
            hits = [
                {
                    "accession": a,
                    "gene": (crop_rows.get(a) or {}).get("gene", ""),
                    "protein_name": (crop_rows.get(a) or {}).get("protein_name", ""),
                    "source": "curated_hook",
                }
                for a in other
            ]

        score = tier = None
        if hits:
            matched += 1
            sc = crop_scores.get(hits[0]["accession"])
            if sc:
                matched_scored += 1
                score = sc["score"]
                tier = sc["tier"]
            for h in hits:
                h["score"] = (crop_scores.get(h["accession"]) or {}).get("score")
                h["tier"] = (crop_scores.get(h["accession"]) or {}).get("tier")
            grape_to_crop[pid] = hits

        # Keep thin ortholog_{col}* columns as the crop join source of truth.
        # Do NOT expand ortholog_json with full second-crop hit lists (GitHub ~100 MiB
        # SQLite limit); Ath/rice compact entries stay as deposited by enrich_orthologs.py.
        conn.execute(
            f"""
            UPDATE proteins SET
                ortholog_{col} = ?,
                ortholog_{col}_gene = ?,
                ortholog_{col}_score = ?,
                ortholog_{col}_tier = ?
            WHERE protein_id = ?
            """,
            (
                hits[0]["accession"] if hits else None,
                hits[0].get("gene") if hits else None,
                score,
                tier,
                pid,
            ),
        )
        if hook and hook.get("family"):
            conn.execute(
                """
                UPDATE proteins SET cross_crop_family = COALESCE(cross_crop_family, ?)
                WHERE protein_id = ?
                """,
                (hook.get("family"), pid),
            )

    summary = {
        "crop": crop["slug"],
        "species": crop["species"],
        "proteome": crop["proteome"],
        "n_grape": len(ids),
        "n_matched": matched,
        "n_matched_with_crop_score": matched_scored,
        "n_hook_overlay": matched_hook,
        "n_crop_scored": len(crop_scores),
        "crop_db": str(crop_db_path.name) if crop_db_path.exists() else None,
        "columns": [
            f"ortholog_{col}",
            f"ortholog_{col}_gene",
            f"ortholog_{col}_score",
            f"ortholog_{col}_tier",
        ],
    }
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (f"second_crop_{crop['slug']}", json.dumps(summary)),
    )
    # keep hooks blob fresh for UI
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("cross_crop_hooks", json.dumps(hooks)),
    )
    conn.commit()

    # reverse pointers on crop DB + gene/protein names from OrthoDB stream
    if crop_db_path.exists() and crop_scores:
        cconn = sqlite3.connect(str(crop_db_path))
        cols = {r[1] for r in cconn.execute("PRAGMA table_info(proteins)")}
        for name, decl in (
            ("gene_name", "TEXT"),
            ("protein_name", "TEXT"),
            ("grape_ortholog", "TEXT"),
            ("grape_ortholog_gene", "TEXT"),
            ("orthodb_id", "TEXT"),
            ("score", "REAL"),
        ):
            if name not in cols and name != "score":
                cconn.execute(f"ALTER TABLE proteins ADD COLUMN {name} {decl}")
        # invert grape→crop to crop→grape (first grape wins)
        crop_to_grape: dict[str, tuple[str, str, str | None]] = {}
        for gpid, hits in grape_to_crop.items():
            ggene = (grape.get(gpid) or {}).get("gene", "")
            for h in hits:
                cid = h["accession"]
                if cid not in crop_to_grape:
                    odb = None
                    for oid in (grape.get(gpid) or {}).get("orthodb") or []:
                        if cid in crop_inv.get(oid, []):
                            odb = oid
                            break
                    crop_to_grape[cid] = (gpid, ggene, odb)
        for acc, info in crop_rows.items():
            gpid, ggene, odb = crop_to_grape.get(acc, (None, None, None))
            if gpid is None and info.get("orthodb"):
                # leave grape null; still fill names
                pass
            cconn.execute(
                """
                UPDATE proteins SET
                    gene_name = COALESCE(?, gene_name),
                    protein_name = COALESCE(?, protein_name),
                    grape_ortholog = COALESCE(?, grape_ortholog),
                    grape_ortholog_gene = COALESCE(?, grape_ortholog_gene),
                    orthodb_id = COALESCE(?, orthodb_id)
                WHERE protein_id = ?
                """,
                (
                    info.get("gene") or None,
                    info.get("protein_name") or None,
                    gpid,
                    ggene or None,
                    odb,
                    acc,
                ),
            )
        cconn.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            ("grape_join", json.dumps({"n_with_grape_ortholog": len(crop_to_grape), **summary})),
        )
        cconn.commit()
        cconn.close()

    # write join summary sidecar
    SECOND.mkdir(parents=True, exist_ok=True)
    join_path = SECOND / f"{crop['slug']}_grape_ortholog_join.json"
    join_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    for row in conn.execute(
        f"SELECT protein_id, ortholog_{col}, ortholog_{col}_gene, ortholog_{col}_score, "
        f"ortholog_{col}_tier, cross_crop_family FROM proteins "
        f"WHERE protein_id='F6HET4' OR cross_crop_family IS NOT NULL LIMIT 8"
    ):
        print(" ", row, flush=True)
    conn.close()


if __name__ == "__main__":
    main()
