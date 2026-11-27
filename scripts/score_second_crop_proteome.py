#!/usr/bin/env python3
"""Score a second-crop UniProt reference proteome with PlantLLPS (composition proxy).

Default crop: soybean Glycine max UP000008827. Writes TSV.gz + slim SQLite under
data/second_crop/ (keeps grape_condensatome.sqlite under GitHub's ~100 MB limit).

Example::

    python3 scripts/score_second_crop_proteome.py
    python3 scripts/score_second_crop_proteome.py --crop tomato
    python3 scripts/score_second_crop_proteome.py --limit 2000   # smoke / subset
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sqlite3
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from _http_util import DEFAULT_UA  # noqa: E402
from plantllps.features import FEATURE_NAMES  # noqa: E402
from plantllps.model import PlantLLPSModel, assign_tier, percentile_cutoffs  # noqa: E402

DATA = ROOT / "data"
OUT_DIR = DATA / "second_crop"
DEFAULT_MODEL = DATA / "model_formal.json"

CROPS = {
    "soybean": {
        "proteome": "UP000008827",
        "taxon": 3847,
        "species": "Glycine max",
        "slug": "soybean",
        "label": "soybean",
        "odb_key": "gma",
    },
    "tomato": {
        "proteome": "UP000004994",
        "taxon": 4081,
        "species": "Solanum lycopersicum",
        "slug": "tomato",
        "label": "tomato",
        "odb_key": "sly",
    },
}

TSV_COLUMNS = [
    "protein_id",
    "length",
    "score_plantllps",
    "score_companion",
    "score_logistic",
    "tier",
    "percentile",
    *FEATURE_NAMES,
]


def download_fasta(proteome: str, dest: Path, force: bool = False) -> Path:
    if dest.exists() and not force and dest.stat().st_size > 1000:
        print(f"Using cache {dest} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    q = urllib.parse.urlencode(
        {
            "compressed": "true",
            "format": "fasta",
            "query": f"(proteome:{proteome})",
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/stream?{q}"
    print(f"Downloading FASTA {proteome} → {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": DEFAULT_UA})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=900) as resp, open(dest, "wb") as out:
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
    raise SystemExit(f"FASTA download failed for {proteome}")


def open_fasta(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_uniprot_fasta(path: Path, limit: int = 0) -> list[tuple[str, str]]:
    """Return (accession, sequence) list; UniProt headers use sp|ACC|NAME."""
    records: list[tuple[str, str]] = []
    pid = None
    chunks: list[str] = []
    with open_fasta(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if pid and chunks:
                    records.append((pid, "".join(chunks)))
                    if limit and len(records) >= limit:
                        return records
                header = line[1:].split()[0]
                parts = header.split("|")
                pid = parts[1] if len(parts) > 1 else header
                chunks = []
            else:
                chunks.append(line.strip())
        if pid and chunks and (not limit or len(records) < limit):
            records.append((pid, "".join(chunks)))
    return records


def write_slim_sqlite(
    db_path: Path,
    rows: list[dict],
    *,
    crop: dict,
    model_path: Path,
    honesty: str,
) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
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
        "crop": crop["slug"],
        "species": crop["species"],
        "proteome": crop["proteome"],
        "taxon": crop["taxon"],
        "n_scored": len(rows),
        "model": str(model_path.name),
        "disorder": "composition_proxy",
        "honesty": honesty,
    }
    for k, v in meta.items():
        conn.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?)",
            (k, json.dumps(v) if not isinstance(v, str) else v),
        )
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?)",
        ("summary_json", json.dumps(meta)),
    )
    conn.executemany(
        """
        INSERT INTO proteins(
            protein_id, length, score, score_companion, score_logistic,
            tier, percentile
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                r["protein_id"],
                int(r["length"]),
                float(r["score_plantllps"]),
                float(r["score_companion"]),
                float(r["score_logistic"]),
                r["tier"],
                float(r["percentile"]) if r.get("percentile") != "" else None,
            )
            for r in rows
        ],
    )
    conn.commit()
    conn.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--crop", choices=sorted(CROPS), default="soybean")
    ap.add_argument("--model", default=str(DEFAULT_MODEL))
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    ap.add_argument("--fasta", default="", help="Optional local FASTA(.gz); else download")
    ap.add_argument("--force-download", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="Score only first N (smoke / subset)")
    ap.add_argument("--min-length", type=int, default=80)
    ap.add_argument("--tier-a-pct", type=float, default=95.0)
    ap.add_argument("--tier-c-pct", type=float, default=50.0)
    args = ap.parse_args()

    crop = CROPS[args.crop]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    slug = crop["slug"]
    fasta = Path(args.fasta) if args.fasta else out_dir / f"{crop['proteome']}.fasta.gz"
    if not args.fasta:
        download_fasta(crop["proteome"], fasta, force=args.force_download)

    print(f"Parsing {fasta} …", flush=True)
    records = parse_uniprot_fasta(fasta, limit=args.limit or 0)
    min_len = int(args.min_length)
    clean: list[tuple[str, str]] = []
    for pid, seq in records:
        if len(seq) < min_len:
            continue
        aa = sum(1 for c in seq.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
        if aa / max(1, len(seq)) >= 0.95:
            clean.append((pid, seq))
    records = clean
    print(f"Scoring {len(records)} {crop['species']} proteins …", flush=True)

    model = PlantLLPSModel.load(args.model)
    t0 = time.time()
    rows: list[dict] = []
    for i, (pid, seq) in enumerate(records, 1):
        scored = model.score(seq, use_metapredict=False)
        rows.append({"protein_id": pid, "percentile": "", **scored})
        if i % 2000 == 0:
            elapsed = time.time() - t0
            rate = i / max(elapsed, 1e-6)
            eta = (len(records) - i) / max(rate, 1e-6)
            print(f"  {i}/{len(records)} ({rate:.0f}/s, ETA {eta/60:.1f} min)", flush=True)

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
        avg = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = avg
        i = j
    n = len(rows)
    tier_counts = {"A": 0, "B": 0, "C": 0, "D": 0}
    for idx, r in enumerate(rows):
        feats = {k: r[k] for k in FEATURE_NAMES if k in r}
        r["tier"] = assign_tier(r["score_plantllps"], feats, cutoffs)
        r["percentile"] = round(100.0 * ranks[idx] / n, 3)
        tier_counts[r["tier"]] = tier_counts.get(r["tier"], 0) + 1

    honesty = (
        f"Second-crop condensatome for {crop['species']} ({crop['proteome']}). "
        "Scores use the grape-trained PlantLLPS v13 hybrid with composition-proxy disorder "
        "(not metapredict). Tier cutoffs recalibrated on this proteome's score distribution. "
        "Not a species-specific retrained model."
    )
    if args.limit:
        honesty += f" SUBSET deposit: --limit {args.limit} (not full proteome)."

    tsv_path = out_dir / f"{slug}_condensatome.proteome.tsv.gz"
    with gzip.open(tsv_path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in TSV_COLUMNS})

    db_path = out_dir / f"{slug}_condensatome.sqlite"
    write_slim_sqlite(db_path, rows, crop=crop, model_path=Path(args.model), honesty=honesty)

    elapsed = time.time() - t0
    summary = {
        "crop": slug,
        "species": crop["species"],
        "proteome": crop["proteome"],
        "taxon": crop["taxon"],
        "n_fasta_parsed": len(records),
        "n_scored": len(rows),
        "limit": args.limit or None,
        "min_length": min_len,
        "tiers": tier_counts,
        "tier_cutoffs": cutoffs,
        "wall_seconds": round(elapsed, 1),
        "model": Path(args.model).name,
        "disorder": "composition_proxy",
        "tsv": str(tsv_path.name),
        "sqlite": str(db_path.name),
        "honesty": honesty,
        "mean_score": round(sum(scores) / max(1, len(scores)), 6),
    }
    summary_path = out_dir / f"{slug}_condensatome.summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    print(f"TSV → {tsv_path} ({tsv_path.stat().st_size} bytes)", flush=True)
    print(f"SQLite → {db_path} ({db_path.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    main()
