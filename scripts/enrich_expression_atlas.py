#!/usr/bin/env python3
"""Join EMBL-EBI Expression Atlas TPMs for Vitis vinifera into GrapeCondensatome SQLite.

Default experiment: **E-MTAB-8758** (Pinot Noir berry RNA-seq baseline; GeneIDs are
``VIT_*`` loci matching UniProt ``gene_name`` for UP000009183).

Downloads (unless cached)::

    https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-8758/E-MTAB-8758-tpms.tsv
    https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-8758/E-MTAB-8758-configuration.xml

Adds columns: expr_experiment, expr_mean_tpm, expr_max_tpm, expr_tpm_json.
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import sys
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
DEFAULT_DIR = ROOT / "data" / "expression_atlas"
EXPERIMENT = "E-MTAB-8758"
FTP_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/"
    + EXPERIMENT
)


def download(url: str, dest: Path, force: bool = False) -> Path:
    if dest.exists() and not force and dest.stat().st_size > 1000:
        print(f"Using cache → {dest}", flush=True)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url}", flush=True)
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "GrapeCondensatome/0.3 (research; Expression Atlas join)"},
    )
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=300) as resp, open(dest, "wb") as out:
                while True:
                    chunk = resp.read(1024 * 256)
                    if not chunk:
                        break
                    out.write(chunk)
            return dest
        except (urllib.error.URLError, TimeoutError) as e:
            wait = 4 * (2**attempt)
            print(f"  retry in {wait}s ({e})", flush=True)
            time.sleep(wait)
    raise SystemExit(f"download failed: {url}")


def parse_assay_labels(config_xml: Path) -> dict[str, str]:
    labels: dict[str, str] = {}
    if not config_xml.exists():
        return labels
    root = ET.parse(config_xml).getroot()
    for ag in root.findall(".//assay_group"):
        gid = ag.attrib.get("id")
        lab = ag.attrib.get("label") or gid
        if gid:
            labels[gid] = lab
    return labels


def _parse_tpm_cell(cell: str) -> float:
    """Atlas TPMs are comma-separated replicate values; use arithmetic mean."""
    cell = (cell or "").strip()
    if not cell:
        return 0.0
    vals = []
    for part in cell.split(","):
        part = part.strip()
        if not part:
            continue
        try:
            vals.append(float(part))
        except ValueError:
            continue
    return float(mean(vals)) if vals else 0.0


def load_tpm_by_gene(tpm_path: Path, assay_labels: dict[str, str]) -> dict[str, dict]:
    """GeneID → {mean, max, by_group: {label: tpm}}."""
    out: dict[str, dict] = {}
    with tpm_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        group_cols = [c for c in (reader.fieldnames or []) if c not in {"GeneID", "Gene Name"}]
        for row in reader:
            gid = (row.get("GeneID") or "").strip()
            if not gid or gid.startswith("ENSRNA"):
                continue
            by_group = {}
            vals = []
            for g in group_cols:
                v = _parse_tpm_cell(row.get(g, ""))
                label = assay_labels.get(g, g)
                by_group[label] = round(v, 4)
                vals.append(v)
            if not vals:
                continue
            out[gid] = {
                "mean": round(mean(vals), 4),
                "max": round(max(vals), 4),
                "by_group": by_group,
            }
    return out


def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    for col, decl in (
        ("expr_experiment", "TEXT"),
        ("expr_mean_tpm", "REAL"),
        ("expr_max_tpm", "REAL"),
        ("expr_tpm_json", "TEXT"),
    ):
        if col not in cols:
            conn.execute(f"ALTER TABLE proteins ADD COLUMN {col} {decl}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--data-dir", default=str(DEFAULT_DIR))
    ap.add_argument("--force-download", action="store_true")
    ap.add_argument("--experiment", default=EXPERIMENT)
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    tpm_path = data_dir / f"{args.experiment}-tpms.tsv"
    cfg_path = data_dir / f"{args.experiment}-configuration.xml"
    if args.experiment == EXPERIMENT:
        download(f"{FTP_BASE}/{EXPERIMENT}-tpms.tsv", tpm_path, force=args.force_download)
        download(
            f"{FTP_BASE}/{EXPERIMENT}-configuration.xml",
            cfg_path,
            force=args.force_download,
        )
    elif not tpm_path.exists():
        raise SystemExit(f"TPM file not found: {tpm_path}")

    labels = parse_assay_labels(cfg_path)
    expr = load_tpm_by_gene(tpm_path, labels)
    print(f"Loaded expression for {len(expr)} GeneIDs from {tpm_path.name}", flush=True)

    conn = sqlite3.connect(str(args.db))
    ensure_columns(conn)
    rows = list(conn.execute("SELECT protein_id, gene_name FROM proteins"))
    matched = 0
    for pid, gene in rows:
        gene = (gene or "").strip()
        hit = expr.get(gene)
        if not hit:
            continue
        matched += 1
        conn.execute(
            """
            UPDATE proteins SET
                expr_experiment = ?,
                expr_mean_tpm = ?,
                expr_max_tpm = ?,
                expr_tpm_json = ?
            WHERE protein_id = ?
            """,
            (
                args.experiment,
                hit["mean"],
                hit["max"],
                json.dumps(hit["by_group"], sort_keys=True),
                pid,
            ),
        )
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (
            "expression_atlas",
            json.dumps(
                {
                    "experiment": args.experiment,
                    "n_genes_in_atlas": len(expr),
                    "n_proteins_matched": matched,
                    "source": f"{FTP_BASE}/{args.experiment}-tpms.tsv",
                }
            ),
        ),
    )
    conn.commit()
    print(
        f"DB {args.db}: matched={matched}/{len(rows)} "
        f"({100.0 * matched / max(1, len(rows)):.1f}%)",
        flush=True,
    )
    for row in conn.execute(
        "SELECT protein_id, gene_name, expr_mean_tpm, expr_max_tpm FROM proteins "
        "WHERE protein_id='F6HET4' OR expr_mean_tpm IS NOT NULL "
        "ORDER BY expr_mean_tpm DESC LIMIT 5"
    ):
        print(" ", row)
    conn.close()


if __name__ == "__main__":
    main()
