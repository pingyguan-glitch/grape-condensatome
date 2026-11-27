#!/usr/bin/env python3
"""Build a Zenodo deposit bundle for PlantLLPS + GrapeCondensatome.

Copies model weights, SQLite atlases, optional Parquet/TSV exports, CITATION.cff,
LICENSE, checksums, and a deposit README into ``dist/zenodo_<version>/``.

Does **not** invent or mint a DOI. After uploading to Zenodo, paste the real DOI
into ``CITATION.cff`` (see ``docs/ZENODO.md``).

Stdlib-only by default. Parquet export is optional:

    pip install pyarrow
    python3 scripts/export_zenodo_bundle.py --parquet
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sqlite3
import sys
import tarfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
DEFAULT_MODEL = ROOT / "data" / "model_formal.json"
DEFAULT_VERSION = "1.1.0"

# Slim columnar export (skip huge JSON blobs that are already in SQLite).
PROTEINS_PARQUET_COLS = [
    "protein_id",
    "gene_name",
    "protein_name",
    "length",
    "score",
    "score_companion",
    "score_logistic",
    "tier",
    "go_terms",
    "keywords",
    "expr_experiment",
    "expr_mean_tpm",
    "expr_max_tpm",
    "score_psap",
    "score_llphyscore_percentile",
    "plaac_nllr_approx",
    "plaac_qn_frac",
    "af_mean_plddt",
    "af_frac_very_low",
    "af_frac_low",
    "af_frac_confident",
    "af_frac_very_high",
    "af_model_id",
    "af_n_disorder_residues",
    "ortholog_ath",
    "ortholog_ath_gene",
    "ortholog_osa",
    "ortholog_osa_gene",
    "orthodb_id",
    "cross_crop_family",
    "frac_disorder_metapredict",
    "domain_pfam_ids",
    "domain_interpro_ids",
    "idr_n_residues",
    "idr_sources",
    "llps_gold",
    "string_id",
    "string_n_neighbors",
    "string_enrichment_fold",
    "string_hypergeom_pvalue",
    "ortholog_gma",
    "ortholog_gma_gene",
    "ortholog_gma_score",
    "ortholog_gma_tier",
]

CORE_FILES = [
    ("CITATION.cff", "CITATION.cff"),
    ("LICENSE", "LICENSE"),
    ("README.md", "REPO_README.md"),
    ("data/model_formal.json", "model_formal.json"),
    ("data/model_formal.metrics.v13.json", "model_formal.metrics.v13.json"),
    ("data/train_formal_v13.csv", "train_formal_v13.csv"),
    ("data/proteome_summary.json", "proteome_summary.json"),
    ("data/cross_crop_hooks.json", "cross_crop_hooks.json"),
]

OPTIONAL_FILES = [
    ("data/grape_condensatome.proteome.tsv", "grape_condensatome.proteome.tsv"),
    ("data/map_screen/map_screen_spearman.json", "map_screen/map_screen_spearman.json"),
    ("data/map_screen/README.md", "map_screen/README.md"),
    ("data/curated/llps_experimental_gold.tsv", "curated/llps_experimental_gold.tsv"),
    ("data/benchmarks/family_holdout_metrics.json", "benchmarks/family_holdout_metrics.json"),
    ("data/second_crop/soybean_condensatome.sqlite", "second_crop/soybean_condensatome.sqlite"),
    ("data/second_crop/soybean_condensatome.summary.json", "second_crop/soybean_condensatome.summary.json"),
    ("data/second_crop/soybean_grape_ortholog_join.json", "second_crop/soybean_grape_ortholog_join.json"),
    ("data/string/enrichment_summary.json", "string/enrichment_summary.json"),
]


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def copy_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def write_deposit_readme(out_dir: Path, version: str, n_proteins: int) -> None:
    text = f"""# PlantLLPS + GrapeCondensatome — Zenodo deposit (v{version})

This folder is the **deposit package** prepared by
`scripts/export_zenodo_bundle.py`. It is **not** a DOI.

## Contents

| File | Description |
|------|-------------|
| `model_formal.json` | PlantLLPS **v13** hybrid logistic + companion weights |
| `model_formal.metrics.v13.json` | Training metrics snapshot |
| `train_formal_v13.csv` | Formal training set (n=162) |
| `grape_condensatome.sqlite` | *Vitis vinifera* UP000009183 atlas (**{n_proteins}** proteins) |
| `grape_condensatome.proteins.tsv.gz` | Slim proteins-table export (TSV) |
| `grape_condensatome.proteins.parquet` | Same slim table (optional; needs pyarrow at export) |
| `second_crop/` | Soybean condensatome SQLite + join JSON (when present) |
| `CITATION.cff` | Citation metadata (**DOI field empty until you publish**) |
| `LICENSE` | MIT |
| `SHA256SUMS.txt` | Checksums of every file in this bundle |
| `bundle_manifest.json` | Machine-readable inventory |

## Honesty

- PlantLLPS gold-panel AUCs are **partially in-sample** — prefer family hold-out / transfer metrics.
- Default disorder in model weights is a **composition proxy** (metapredict is an annotation track).
- MAP-screen Spearman ρ in this deposit is from an **annotation reconstruction**, not the unpublished companion Supplemental Data S1 table.
- No fabricated public hosted URL ships with this deposit; see `docs/HOSTING.md` in the GitHub repo.

## How to get a real DOI

Follow **`docs/ZENODO.md`** in the repository:

1. Upload this folder (or the `.zip` / `.tar.gz` archive) to Zenodo.
2. Publish the record — Zenodo mints `10.5281/zenodo.<id>`.
3. Paste that DOI into `CITATION.cff` → `identifiers` and re-tag the GitHub release if needed.

**Do not invent a DOI.** Until Zenodo publishes the record, cite the GitHub tag `v{version}`
and this software release via `CITATION.cff` without a Zenodo identifier.

## Source

- GitHub: https://github.com/pingyguan-glitch/grape-condensatome
- Tag / version: `v{version}`
- Generated (UTC): {datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
"""
    (out_dir / "README.md").write_text(text, encoding="utf-8")


def export_proteins_tsv_gz(db: Path, out_path: Path, cols: list[str]) -> int:
    import gzip

    con = sqlite3.connect(str(db))
    con.row_factory = sqlite3.Row
    # Keep only columns that exist
    existing = {r[1] for r in con.execute("PRAGMA table_info(proteins)")}
    use = [c for c in cols if c in existing]
    if not use:
        con.close()
        raise SystemExit("No matching proteins columns for export")
    sql = f"SELECT {', '.join(use)} FROM proteins ORDER BY protein_id"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(out_path, "wt", encoding="utf-8", newline="") as fh:
        import csv

        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(use)
        for row in con.execute(sql):
            w.writerow([row[c] for c in use])
            n += 1
    con.close()
    return n


def export_proteins_parquet(db: Path, out_path: Path, cols: list[str]) -> int:
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise SystemExit(
            "Parquet export requires pyarrow. Install with: pip install pyarrow\n"
            "Or re-run without --parquet (TSV.gz is always written)."
        ) from exc

    con = sqlite3.connect(str(db))
    existing = {r[1] for r in con.execute("PRAGMA table_info(proteins)")}
    use = [c for c in cols if c in existing]
    rows = list(con.execute(f"SELECT {', '.join(use)} FROM proteins ORDER BY protein_id"))
    con.close()
    arrays = []
    for i, _name in enumerate(use):
        col = [r[i] for r in rows]
        if all(v is None or isinstance(v, bool) for v in col):
            arrays.append(pa.array(col, type=pa.bool_()))
        elif all(v is None or (isinstance(v, int) and not isinstance(v, bool)) for v in col):
            arrays.append(pa.array(col, type=pa.int64()))
        elif all(v is None or isinstance(v, (int, float)) and not isinstance(v, bool) for v in col):
            arrays.append(pa.array([None if v is None else float(v) for v in col], type=pa.float64()))
        else:
            arrays.append(pa.array([None if v is None else str(v) for v in col], type=pa.string()))
    table = pa.Table.from_arrays(arrays, names=use)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out_path, compression="zstd")
    return table.num_rows


def write_checksums(out_dir: Path) -> Path:
    lines = []
    for path in sorted(p for p in out_dir.rglob("*") if p.is_file()):
        if path.name in {"SHA256SUMS.txt", "bundle_manifest.json"}:
            continue
        rel = path.relative_to(out_dir).as_posix()
        lines.append(f"{sha256_file(path)}  {rel}")
    out = out_dir / "SHA256SUMS.txt"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def make_archives(out_dir: Path, stem: str) -> list[Path]:
    archives: list[Path] = []
    zip_path = out_dir.parent / f"{stem}.zip"
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for path in sorted(p for p in out_dir.rglob("*") if p.is_file()):
            zf.write(path, arcname=f"{stem}/{path.relative_to(out_dir).as_posix()}")
    archives.append(zip_path)

    tar_path = out_dir.parent / f"{stem}.tar.gz"
    if tar_path.exists():
        tar_path.unlink()
    with tarfile.open(tar_path, "w:gz") as tf:
        tf.add(out_dir, arcname=stem)
    archives.append(tar_path)
    return archives


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--version", default=DEFAULT_VERSION, help="Bundle version label (default: 1.1.0)")
    ap.add_argument("--db", type=Path, default=DEFAULT_DB, help="GrapeCondensatome SQLite path")
    ap.add_argument("--model", type=Path, default=DEFAULT_MODEL, help="PlantLLPS model JSON")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: dist/zenodo_<version>/)",
    )
    ap.add_argument("--parquet", action="store_true", help="Also write proteins Parquet (needs pyarrow)")
    ap.add_argument("--no-archive", action="store_true", help="Skip .zip / .tar.gz archives")
    ap.add_argument("--dry-run", action="store_true", help="Print planned copies; write nothing")
    args = ap.parse_args(argv)

    version = args.version.lstrip("v")
    out_dir = args.out_dir or (ROOT / "dist" / f"zenodo_{version}")
    stem = f"grape-condensatome-zenodo-v{version}"

    if not args.db.is_file():
        raise SystemExit(f"Missing SQLite: {args.db}")
    if not args.model.is_file():
        raise SystemExit(f"Missing model: {args.model}")

    planned: list[tuple[Path, Path]] = []
    for src_rel, dst_rel in CORE_FILES + OPTIONAL_FILES:
        src = ROOT / src_rel
        if src.is_file():
            planned.append((src, out_dir / dst_rel))
        elif (src_rel, dst_rel) in CORE_FILES:
            raise SystemExit(f"Missing required file: {src}")

    # Always include grape SQLite under a stable name
    planned.append((args.db.resolve(), out_dir / "grape_condensatome.sqlite"))
    # Allow --model override when it differs from the default core copy
    default_model = (ROOT / "data" / "model_formal.json").resolve()
    if args.model.resolve() != default_model:
        planned = [(s, d) for s, d in planned if d.name != "model_formal.json"]
        planned.append((args.model.resolve(), out_dir / "model_formal.json"))

    # De-dupe by destination (last source wins)
    by_dst: dict[Path, Path] = {dst: src for src, dst in planned}
    planned = [(src, dst) for dst, src in by_dst.items()]

    print(f"Bundle → {out_dir}")
    for src, dst in planned:
        try:
            src_show = src.relative_to(ROOT)
        except ValueError:
            src_show = src
        print(f"  {src_show} → {dst.relative_to(out_dir)}")

    if args.dry_run:
        print("Dry run only.")
        return 0

    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for src, dst in planned:
        copy_file(src, dst)

    con = sqlite3.connect(str(args.db))
    n_proteins = int(con.execute("SELECT COUNT(*) FROM proteins").fetchone()[0])
    con.close()

    tsv_path = out_dir / "grape_condensatome.proteins.tsv.gz"
    n_tsv = export_proteins_tsv_gz(args.db, tsv_path, PROTEINS_PARQUET_COLS)
    print(f"Wrote {tsv_path.name} ({n_tsv} rows)")

    parquet_path = None
    if args.parquet:
        parquet_path = out_dir / "grape_condensatome.proteins.parquet"
        n_pq = export_proteins_parquet(args.db, parquet_path, PROTEINS_PARQUET_COLS)
        print(f"Wrote {parquet_path.name} ({n_pq} rows)")

    write_deposit_readme(out_dir, version, n_proteins)

    files = sorted(
        p.relative_to(out_dir).as_posix()
        for p in out_dir.rglob("*")
        if p.is_file() and p.name not in {"SHA256SUMS.txt", "bundle_manifest.json"}
    )
    write_checksums(out_dir)

    manifest = {
        "title": "PlantLLPS and GrapeCondensatome",
        "version": version,
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "repository": "https://github.com/pingyguan-glitch/grape-condensatome",
        "n_proteins": n_proteins,
        "doi": None,
        "doi_note": "Do not invent a DOI. Paste the Zenodo DOI after publish (docs/ZENODO.md).",
        "files": files,
        "exports": {
            "sqlite": "grape_condensatome.sqlite",
            "proteins_tsv_gz": "grape_condensatome.proteins.tsv.gz",
            "proteins_parquet": "grape_condensatome.proteins.parquet" if parquet_path else None,
            "model": "model_formal.json",
        },
    }
    (out_dir / "bundle_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    archives: list[Path] = []
    if not args.no_archive:
        archives = make_archives(out_dir, stem)
        for a in archives:
            print(f"Archive: {a} ({a.stat().st_size / (1 << 20):.1f} MiB)")

    print("Done.")
    print("Next: upload the folder or archive to Zenodo — see docs/ZENODO.md")
    print("DOI status: NOT SET (paste after Zenodo publish)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
