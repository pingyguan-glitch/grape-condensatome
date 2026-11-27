#!/usr/bin/env python3
"""Join fruit-crop condensatomes to grape via shared OrthoDB groups.

Unlike enrich_second_crop_orthologs.py (soy/tomato), this script **does not**
add per-species columns to grape_condensatome.sqlite (keeps grape DB ≤100 MiB).
It only:

1. Streams UniProt OrthoDB xrefs for grape + each fruit crop.
2. Writes reverse grape_ortholog / gene / orthodb_id onto fruit SQLite.
3. Optionally copies Arabidopsis accession from grape when the matched grape
   row has ortholog_ath (Ath hook chain: fruit → grape → Ath).
4. Writes join JSON under data/fruit_crops/.

Example::

    python3 scripts/enrich_fruit_crop_orthologs.py --all
    python3 scripts/enrich_fruit_crop_orthologs.py --species apple,peach
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
DEFAULT_MANIFEST = DATA / "fruit_crops" / "species_manifest.yaml"
DEFAULT_GRAPE = DATA / "grape_condensatome.sqlite"
CACHE = DATA / "orthologs"
FRUIT = DATA / "fruit_crops"
GRAPE_PROTEOME = "UP000009183"


def load_manifest(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in {".yaml", ".yml"}:
        import yaml

        return yaml.safe_load(text)
    return json.loads(text)


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
    raise RuntimeError(f"download failed for {proteome}")


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


def load_grape_ath(conn: sqlite3.Connection) -> dict[str, tuple[str | None, str | None]]:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    if "ortholog_ath" not in cols:
        return {}
    out = {}
    for r in conn.execute(
        "SELECT protein_id, ortholog_ath, ortholog_ath_gene FROM proteins"
    ):
        out[r[0]] = (r[1], r[2])
    return out


def enrich_one(
    crop: dict,
    *,
    grape: dict[str, dict],
    grape_inv: dict[str, list[str]],
    grape_ath: dict[str, tuple[str | None, str | None]],
    crop_db_path: Path,
    cache: Path,
    force: bool,
) -> dict:
    slug = crop["slug"]
    prefixes = tuple(crop.get("gene_prefixes") or ("LOC",))
    if not crop_db_path.exists():
        return {
            "slug": slug,
            "ok": False,
            "error": f"missing crop DB {crop_db_path}",
        }

    try:
        crop_odb_path = stream_orthodb(
            crop["proteome"], cache / f"{slug}_orthodb.tsv.gz", force=force
        )
    except Exception as e:
        return {"slug": slug, "ok": False, "error": f"OrthoDB download: {e}"}

    crop_rows = parse_stream(crop_odb_path)
    crop_inv = invert_odb(crop_rows)
    print(
        f"[{slug}] OrthoDB rows={len(crop_rows)} groups={len(crop_inv)}",
        flush=True,
    )

    # Map crop accession → best grape hit via shared OrthoDB group
    crop_to_grape: dict[str, tuple[str, str, str | None]] = {}
    for oid, crop_accs in crop_inv.items():
        grape_accs = grape_inv.get(oid) or []
        if not grape_accs:
            continue
        g_hits = pick_best(grape_accs, grape, ("VIT_", "LOC"), limit=1)
        if not g_hits:
            continue
        gpid = g_hits[0]["accession"]
        ggene = grape[gpid]["gene"]
        for cid in crop_accs:
            if cid not in crop_to_grape:
                crop_to_grape[cid] = (gpid, ggene, oid)

    cconn = sqlite3.connect(str(crop_db_path))
    cols = {r[1] for r in cconn.execute("PRAGMA table_info(proteins)")}
    for name, decl in (
        ("gene_name", "TEXT"),
        ("protein_name", "TEXT"),
        ("grape_ortholog", "TEXT"),
        ("grape_ortholog_gene", "TEXT"),
        ("orthodb_id", "TEXT"),
        ("ath_ortholog", "TEXT"),
        ("ath_ortholog_gene", "TEXT"),
    ):
        if name not in cols:
            cconn.execute(f"ALTER TABLE proteins ADD COLUMN {name} {decl}")

    n_grape = n_ath = 0
    for acc, info in crop_rows.items():
        gpid, ggene, odb = crop_to_grape.get(acc, (None, None, None))
        ath_id = ath_gene = None
        if gpid:
            n_grape += 1
            ath_id, ath_gene = grape_ath.get(gpid, (None, None))
            if ath_id:
                n_ath += 1
        cconn.execute(
            """
            UPDATE proteins SET
                gene_name = COALESCE(?, gene_name),
                protein_name = COALESCE(?, protein_name),
                grape_ortholog = COALESCE(?, grape_ortholog),
                grape_ortholog_gene = COALESCE(?, grape_ortholog_gene),
                orthodb_id = COALESCE(?, orthodb_id),
                ath_ortholog = COALESCE(?, ath_ortholog),
                ath_ortholog_gene = COALESCE(?, ath_ortholog_gene)
            WHERE protein_id = ?
            """,
            (
                info.get("gene") or None,
                info.get("protein_name") or None,
                gpid,
                ggene or None,
                odb,
                ath_id,
                ath_gene,
                acc,
            ),
        )

    n_proteins = cconn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
    summary = {
        "slug": slug,
        "species": crop["species"],
        "proteome": crop["proteome"],
        "n_crop_proteins": n_proteins,
        "n_crop_orthodb_rows": len(crop_rows),
        "n_with_grape_ortholog": n_grape,
        "n_with_ath_via_grape": n_ath,
        "ok": True,
        "note": (
            "Fruit→grape OrthoDB join written to crop SQLite only "
            "(no new grape atlas columns). Ath via grape.ortholog_ath when present."
        ),
    }
    cconn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("grape_join", json.dumps(summary)),
    )
    cconn.commit()
    cconn.close()

    join_path = FRUIT / f"{slug}_grape_ortholog_join.json"
    join_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    ap.add_argument("--grape-db", default=str(DEFAULT_GRAPE))
    ap.add_argument("--out-dir", default=str(FRUIT))
    ap.add_argument("--cache-dir", default=str(CACHE))
    ap.add_argument("--species", default="")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--force-download", action="store_true")
    args = ap.parse_args()

    manifest = load_manifest(Path(args.manifest))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    score_species = [s for s in (manifest.get("species") or []) if s.get("status") == "score"]
    if args.all:
        selected = score_species
    elif args.species:
        want = {x.strip() for x in args.species.split(",") if x.strip()}
        selected = [s for s in score_species if s["slug"] in want]
    else:
        raise SystemExit("Specify --all or --species apple,peach")

    # Only enrich species with deposited SQLite
    selected = [
        s
        for s in selected
        if (out_dir / f"{s['slug']}_condensatome.sqlite").exists()
    ]
    if not selected:
        raise SystemExit(f"No fruit crop SQLite found under {out_dir}")

    grape_path = stream_orthodb(
        GRAPE_PROTEOME, cache / "grape_orthodb.tsv.gz", force=args.force_download
    )
    grape = parse_stream(grape_path)
    grape_inv = invert_odb(grape)
    print(f"Grape OrthoDB rows={len(grape)} groups={len(grape_inv)}", flush=True)

    gconn = sqlite3.connect(str(args.grape_db))
    grape_ath = load_grape_ath(gconn)
    gconn.close()
    print(f"Grape Ath ortholog pointers: {sum(1 for v in grape_ath.values() if v[0])}", flush=True)

    results = []
    for crop in selected:
        crop_db = out_dir / f"{crop['slug']}_condensatome.sqlite"
        try:
            results.append(
                enrich_one(
                    crop,
                    grape=grape,
                    grape_inv=grape_inv,
                    grape_ath=grape_ath,
                    crop_db_path=crop_db,
                    cache=cache,
                    force=args.force_download,
                )
            )
        except Exception as e:
            results.append({"slug": crop["slug"], "ok": False, "error": str(e)})
            print(f"[{crop['slug']}] FAILED: {e}", flush=True)

    report = {
        "n_ok": sum(1 for r in results if r.get("ok")),
        "n_failed": sum(1 for r in results if not r.get("ok")),
        "results": results,
    }
    (out_dir / "ortholog_join_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"n_ok": report["n_ok"], "n_failed": report["n_failed"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
