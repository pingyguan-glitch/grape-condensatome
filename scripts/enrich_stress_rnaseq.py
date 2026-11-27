#!/usr/bin/env python3
"""Join Expression Atlas *stress* RNA-seq differentials into GrapeCondensatome SQLite.

Default accessions (public EMBL-EBI Expression Atlas; clear reuse under Atlas terms)::

  * **E-MTAB-4324** — Pinot Noir leaves; Trichoderma harzianum T39 / Plasmopara viticola
    (downy mildew) induced-resistance RNA-seq (Perazzolli et al. 2012).
  * **E-GEOD-67932** — Semillon berries; Botrytis cinerea noble-rot stages vs control
    RNA-seq (Blanco-Ulate et al. 2015).

Both deposit Atlas ``*-analytics.tsv`` with ``VIT_*`` Gene IDs matching UniProt
``gene_name`` for UP000009183 (same join key as baseline E-MTAB-8758).

Writes gzip sidecar ``data/expression_atlas/stress_by_protein.json.gz`` (keeps the
main SQLite under GitHub's 100 MiB limit) and meta ``stress_rnaseq`` with match
rates. Optionally also fills ``expr_stress_json`` when ``--also-sqlite-column``.
Baseline berry TPM (E-MTAB-8758) is unchanged.

Example::

    python3 scripts/enrich_stress_rnaseq.py
    python3 scripts/enrich_stress_rnaseq.py --experiments E-MTAB-4324,E-GEOD-67932
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

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
DEFAULT_DIR = ROOT / "data" / "expression_atlas"
FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"

# Accession → short human label for UI / docs
DEFAULT_EXPERIMENTS = {
    "E-MTAB-4324": {
        "title": "Downy mildew / Trichoderma induced resistance (leaves)",
        "stress": "biotic: Plasmopara viticola ± Trichoderma harzianum T39",
        "technology": "RNA-seq differential",
        "citation": "Perazzolli et al. 2012 (BMC Plant Biol / Atlas E-MTAB-4324)",
        "license_note": "EMBL-EBI Expression Atlas public experiment files",
    },
    "E-GEOD-67932": {
        "title": "Noble rot Botrytis stages vs control (berries)",
        "stress": "biotic: Botrytis cinerea noble rot (stages 1–3)",
        "technology": "RNA-seq differential",
        "citation": "Blanco-Ulate et al. 2015 (Plant Physiol / Atlas E-GEOD-67932)",
        "license_note": "EMBL-EBI Expression Atlas public experiment files",
    },
}

# Atlas analytics often omit adj.p; treat nominal p < 0.05 and |log2FC| >= 1 as call
P_THR = 0.05
LFC_THR = 1.0


def download(url: str, dest: Path, force: bool = False) -> Path:
    if dest.exists() and not force and dest.stat().st_size > 1000:
        print(f"Using cache → {dest}", flush=True)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url}", flush=True)
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "Mozilla/5.0 (compatible; GrapeCondensatome/0.7; "
                "+https://github.com/pingyguan-glitch/grape-condensatome)"
            )
        },
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


def parse_contrast_labels(config_xml: Path) -> dict[str, str]:
    """contrast_id → 'reference_label vs test_label'."""
    labels: dict[str, str] = {}
    assay: dict[str, str] = {}
    if not config_xml.exists():
        return labels
    root = ET.parse(config_xml).getroot()
    for ag in root.findall(".//assay_group"):
        gid = ag.attrib.get("id")
        lab = (ag.attrib.get("label") or gid or "").strip()
        if gid:
            assay[gid] = lab
    for ct in root.findall(".//contrast"):
        cid = ct.attrib.get("id")
        if not cid:
            continue
        ref = (ct.findtext("reference_assay_group") or "").strip()
        test = (ct.findtext("test_assay_group") or "").strip()
        # Atlas contrast id is often reference_test (g4_g1 = none vs pathogen)
        ref_lab = assay.get(ref, ref)
        test_lab = assay.get(test, test)
        labels[cid] = f"{ref_lab} vs {test_lab}"
    return labels


def load_analytics(path: Path, contrast_labels: dict[str, str]) -> dict[str, dict]:
    """GeneID → {comparisons: {label: {log2fc, pvalue, significant}}}."""
    out: dict[str, dict] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        contrasts: list[tuple[str, str, str]] = []
        for f in fields:
            if f.endswith(".log2foldchange"):
                cid = f[: -len(".log2foldchange")]
                pcol = f"{cid}.p-value"
                if pcol not in fields:
                    pcol = f"{cid}.pvalue"
                contrasts.append((cid, f, pcol if pcol in fields else ""))
        for row in reader:
            gid = (row.get("Gene ID") or row.get("GeneID") or "").strip()
            if not gid or gid.startswith("ENSRNA"):
                continue
            comps: dict[str, dict] = {}
            for cid, lfc_col, p_col in contrasts:
                raw_lfc = (row.get(lfc_col) or "").strip()
                if not raw_lfc:
                    continue
                try:
                    lfc = float(raw_lfc)
                except ValueError:
                    continue
                pval = None
                if p_col:
                    raw_p = (row.get(p_col) or "").strip()
                    if raw_p:
                        try:
                            pval = float(raw_p)
                        except ValueError:
                            pval = None
                lab = contrast_labels.get(cid, cid)
                sig = bool(
                    pval is not None and pval < P_THR and abs(lfc) >= LFC_THR
                )
                comps[lab] = {
                    "log2fc": round(lfc, 4),
                    "pvalue": None if pval is None else round(pval, 6),
                    "contrast_id": cid,
                    "significant": sig,
                }
            if comps:
                out[gid] = {"comparisons": comps}
    return out


def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    if "expr_stress_json" not in cols:
        conn.execute("ALTER TABLE proteins ADD COLUMN expr_stress_json TEXT")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--data-dir", default=str(DEFAULT_DIR))
    ap.add_argument(
        "--experiments",
        default=",".join(DEFAULT_EXPERIMENTS.keys()),
        help="Comma-separated Atlas accessions",
    )
    ap.add_argument("--force-download", action="store_true")
    ap.add_argument(
        "--also-sqlite-column",
        action="store_true",
        help="Also write expr_stress_json into SQLite (large; exceeds GH 100MiB with full proteome)",
    )
    args = ap.parse_args()

    accessions = [a.strip() for a in args.experiments.split(",") if a.strip()]
    data_dir = Path(args.data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    per_acc: dict[str, dict] = {}
    for acc in accessions:
        meta = dict(DEFAULT_EXPERIMENTS.get(acc) or {})
        meta.setdefault("title", acc)
        meta.setdefault("technology", "Expression Atlas differential")
        base = f"{FTP_BASE}/{acc}"
        analytics = data_dir / f"{acc}-analytics.tsv"
        config = data_dir / f"{acc}-configuration.xml"
        download(f"{base}/{acc}-analytics.tsv", analytics, force=args.force_download)
        download(f"{base}/{acc}-configuration.xml", config, force=args.force_download)
        labels = parse_contrast_labels(config)
        genes = load_analytics(analytics, labels)
        per_acc[acc] = {
            "meta": meta,
            "genes": genes,
            "n_genes": len(genes),
            "contrast_labels": labels,
            "source": f"{base}/{acc}-analytics.tsv",
        }
        print(
            f"{acc}: {len(genes)} genes, contrasts={list(labels.values()) or list(labels)}",
            flush=True,
        )

    conn = sqlite3.connect(str(args.db), timeout=120.0)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=60000")
    ensure_columns(conn)
    rows = list(conn.execute("SELECT protein_id, gene_name, expr_stress_json FROM proteins"))

    # Merge per protein across accessions
    matched_by_acc = {acc: 0 for acc in accessions}
    n_any = 0
    sidecar: dict[str, dict] = {}
    for pid, gene, prev_raw in rows:
        gene = (gene or "").strip()
        blob: dict = {}
        hit_any = False
        for acc, bundle in per_acc.items():
            hit = bundle["genes"].get(gene) if gene else None
            if not hit:
                continue
            matched_by_acc[acc] += 1
            hit_any = True
            entry = {
                "title": bundle["meta"].get("title"),
                "stress": bundle["meta"].get("stress"),
                "technology": bundle["meta"].get("technology"),
                "comparisons": hit["comparisons"],
                "n_significant": sum(
                    1 for c in hit["comparisons"].values() if c.get("significant")
                ),
            }
            blob[acc] = entry
        if hit_any:
            n_any += 1
            sidecar[pid] = blob
            if args.also_sqlite_column:
                conn.execute(
                    "UPDATE proteins SET expr_stress_json = ? WHERE protein_id = ?",
                    (json.dumps(blob, sort_keys=True), pid),
                )
        elif args.also_sqlite_column and prev_raw:
            conn.execute(
                "UPDATE proteins SET expr_stress_json = NULL WHERE protein_id = ?",
                (pid,),
            )

    import gzip

    sidecar_path = data_dir / "stress_by_protein.json.gz"
    with gzip.open(sidecar_path, "wt", encoding="utf-8") as fh:
        json.dump(sidecar, fh, separators=(",", ":"))
    print(f"Wrote sidecar {sidecar_path} ({sidecar_path.stat().st_size} bytes)", flush=True)

    summary = {
        "experiments": {
            acc: {
                **bundle["meta"],
                "n_genes_in_atlas": bundle["n_genes"],
                "n_proteins_matched": matched_by_acc[acc],
                "match_rate": round(
                    100.0 * matched_by_acc[acc] / max(1, len(rows)), 2
                ),
                "contrasts": bundle["contrast_labels"],
                "source": bundle["source"],
                "call_rule": f"significant if p<{P_THR} and |log2FC|>={LFC_THR}",
            }
            for acc, bundle in per_acc.items()
        },
        "n_proteins": len(rows),
        "n_proteins_with_any_stress": n_any,
        "accessions": accessions,
        "sidecar": str(sidecar_path.relative_to(ROOT)),
        "sqlite_column": bool(args.also_sqlite_column),
        "choice_rationale": (
            "Selected two public Atlas RNA-seq differentials with VIT_* GeneIDs "
            "and clear biotic-stress contrasts (downy mildew / Trichoderma; "
            "Botrytis noble rot). Abiotic drought/salt Atlas entries for grape "
            "are microarray (E-GEOD-31594 / E-GEOD-31677) — not joined as RNA-seq. "
            "Full DE blobs ship as gzip sidecar to keep grape_condensatome.sqlite "
            "under GitHub's 100 MiB file limit."
        ),
    }
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("stress_rnaseq", json.dumps(summary)),
    )
    conn.commit()
    report_path = data_dir / "stress_rnaseq_join_report.json"
    report_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(
        f"DB matched any-stress={n_any}/{len(rows)}; per-acc={matched_by_acc}",
        flush=True,
    )
    for row in conn.execute(
        "SELECT protein_id, gene_name, substr(expr_stress_json,1,120) FROM proteins "
        "WHERE protein_id='F6HET4' OR expr_stress_json IS NOT NULL LIMIT 3"
    ):
        print(" ", row)
    conn.close()


if __name__ == "__main__":
    main()
