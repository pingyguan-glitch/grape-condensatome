#!/usr/bin/env python3
"""Fruit-crop cross-species high-score comparison layer.

Reads existing ``data/fruit_crops/*_condensatome.sqlite`` score tables (no
re-scoring) and deposits conservative vs species-specific high-score tables
under ``data/fruit_crops/cross_species/``.

Grouping keys (priority, documented in the deposit README):

1. OrthoDB group id (``orthodb_id``) when present
2. else grape UniProt accession from OrthoDB fruit→grape join
3. else normalized UniProt gene symbol proxy (``gene:SYMBOL``)
4. else singleton accession (``acc:ID@slug``)

Example::

    python3 scripts/analyze_fruit_cross_species.py
    python3 scripts/analyze_fruit_cross_species.py --min-conserved 6 --max-specific 1
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sqlite3
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FRUIT = ROOT / "data" / "fruit_crops"
OUT_DIR = FRUIT / "cross_species"
HIGH = frozenset({"A", "B"})

# Defaults match deposit README / METHODS.
DEFAULT_MIN_CONSERVED = 5
DEFAULT_MAX_SPECIFIC = 2
DEFAULT_MIN_PRESENT_FOR_PEAK = 3  # high in few, present (any tier) in ≥ this many
BROWSE_CONSERVED_CAP = 500
BROWSE_SPECIFIC_PER_SPECIES = 40


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_index(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def scored_fruit_slugs(index: dict, fruit_dir: Path) -> list[str]:
    """Prefer index.json status=score; fall back to on-disk SQLite deposits."""
    slugs: list[str] = []
    for sp in index.get("species") or []:
        if sp.get("status") == "score" and sp.get("slug"):
            slugs.append(sp["slug"])
    if slugs:
        return sorted(slugs)
    return sorted(
        p.name[: -len("_condensatome.sqlite")]
        for p in fruit_dir.glob("*_condensatome.sqlite")
    )


def species_meta(index: dict, slug: str, fruit_dir: Path) -> dict:
    for sp in index.get("species") or []:
        if sp.get("slug") == slug:
            return {
                "slug": slug,
                "label": sp.get("label") or slug,
                "species": sp.get("species") or "",
                "proteome": sp.get("proteome") or "",
                "taxon": sp.get("taxon"),
            }
    db = fruit_dir / f"{slug}_condensatome.sqlite"
    meta = {"slug": slug, "label": slug, "species": "", "proteome": "", "taxon": None}
    if not db.exists():
        return meta
    con = sqlite3.connect(str(db))
    rows = dict(con.execute("SELECT key, value FROM meta").fetchall())
    con.close()
    meta["label"] = rows.get("label") or rows.get("crop") or slug
    meta["species"] = rows.get("species") or ""
    meta["proteome"] = rows.get("proteome") or ""
    return meta


def normalize_gene(raw: str | None) -> str:
    if not raw:
        return ""
    g = raw.split()[0].strip().upper()
    if g in {"", "-", "NA", "N/A", "NONE", "NULL", "UNKNOWN"}:
        return ""
    return g


def cluster_key(orthodb_id: str | None, grape_ortholog: str | None, gene_name: str | None, protein_id: str, slug: str) -> tuple[str, str]:
    """Return (key, key_type)."""
    odb = (orthodb_id or "").split(";")[0].strip()
    if odb:
        return f"odb:{odb}", "orthodb"
    grape = (grape_ortholog or "").split(";")[0].strip()
    if grape:
        return f"grape:{grape}", "grape_ortholog"
    gene = normalize_gene(gene_name)
    if gene:
        return f"gene:{gene}", "gene_proxy"
    return f"acc:{protein_id}@{slug}", "singleton"


def load_species_proteins(db_path: Path, slug: str) -> list[dict]:
    if not db_path.exists():
        return []
    con = sqlite3.connect(str(db_path))
    con.row_factory = sqlite3.Row
    rows = []
    for r in con.execute(
        """
        SELECT protein_id, gene_name, protein_name, length, score, tier, percentile,
               grape_ortholog, grape_ortholog_gene, orthodb_id, ath_ortholog, ath_ortholog_gene
        FROM proteins
        """
    ):
        rows.append(
            {
                "slug": slug,
                "protein_id": r["protein_id"],
                "gene_name": r["gene_name"] or "",
                "protein_name": (r["protein_name"] or "")[:160],
                "length": r["length"],
                "score": float(r["score"]) if r["score"] is not None else None,
                "tier": r["tier"] or "",
                "percentile": float(r["percentile"]) if r["percentile"] is not None else None,
                "grape_ortholog": r["grape_ortholog"] or "",
                "grape_ortholog_gene": r["grape_ortholog_gene"] or "",
                "orthodb_id": r["orthodb_id"] or "",
                "ath_ortholog": r["ath_ortholog"] or "",
                "ath_ortholog_gene": r["ath_ortholog_gene"] or "",
            }
        )
    con.close()
    return rows


def best_member(members: list[dict]) -> dict:
    return max(
        members,
        key=lambda m: (
            1 if m["tier"] in HIGH else 0,
            m["score"] if m["score"] is not None else -1.0,
            m["percentile"] if m["percentile"] is not None else -1.0,
        ),
    )


def write_tsv_gz(path: Path, cols: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, "") for c in cols})


def methods_readme(
    *,
    n_species: int,
    min_conserved: int,
    max_specific: int,
    min_present: int,
    n_clusters: int,
    n_conserved: int,
    n_specific: int,
    n_with_odb_species: int,
) -> str:
    return f"""# Fruit-crop cross-species high-score comparison

Deposit produced by `scripts/analyze_fruit_cross_species.py` from existing
`data/fruit_crops/*_condensatome.sqlite` scores (**no re-scoring**).

## Panel

- **{n_species}** scored fruit-crop proteomes under `data/fruit_crops/` (status=`score`
  in `index.json`). Linked grape / soybean / tomato atlases are **not** folded into
  these tables (they remain on `/hooks` separately).
- Unavailable taxa stay in the fruit-crop manifest `unavailable` block — not invented.

## High-score definition

- **High** = PlantLLPS Tier **A∪B** on that species' recalibrated cutoffs
  (grape-trained v13 + composition-proxy disorder; not species-retrained).

## Orthology / grouping (approximation)

Cluster keys, in priority order:

1. **`orthodb:…`** — UniProt OrthoDB group id written by
   `scripts/enrich_fruit_crop_orthologs.py` (eukaryote-level; can be broad).
2. **`grape:…`** — fruit→grape OrthoDB join accession when OrthoDB id is missing
   but a grape pointer exists.
3. **`gene:…`** — normalized UniProt primary gene symbol proxy when OrthoDB xrefs
   are absent (common for several fruit proteomes). Gene symbols are **not** true
   orthology; cross-species collisions are rare in this panel.
4. **`acc:…@slug`** — singleton accession when no joinable key exists.

Species with dense OrthoDB xrefs dominate multi-species clusters
(~{n_with_odb_species} / {n_species} species contribute ≥1 OrthoDB-linked Tier A/B
row in this deposit). Proteomes lacking OrthoDB xrefs mostly enter as gene-proxy
or singleton clusters and therefore enrich the **species-specific** table — that is
expected, not a bug.

## Conservative LLPS candidates

- Cluster has high-score members in **≥ {min_conserved}** panel species.
- Ranked by `n_species_high` ↓, then `median_high_score` ↓, then `mean_high_score` ↓.
- Per-species representative = highest-score member in that species (Tier A preferred
  over B only via score/tier sort already used for “best member”).

## Species-specific peaks

A cluster is a **species-specific peak** when:

- high-score in **≤ {max_specific}** panel species, **and**
- either:
  - **panel-contrast peak**: present (any tier) in **≥ {min_present}** species, with
    the non-peak species not Tier A/B; or
  - **orphan high**: high in 1–{max_specific} species and present in fewer than
    {min_present} species (typical for gene-proxy / singleton / OrthoDB-sparse taxa).

Peak species = the species (or species list) carrying the high-score member(s).

## Files

| File | Role |
|------|------|
| `conservative_candidates.tsv.gz` | Conservative cluster table |
| `species_specific_peaks.tsv.gz` | Species-specific peak table |
| `cluster_members.tsv.gz` | Per-member long table for clusters in either class |
| `summary.json` | Counts, thresholds, honesty, species list |
| `browse.json` | UI-sized slice (top conserved + per-species peaks) |
| `README.md` | This methods note |

## Honesty

- Ranks are **hypotheses** from cross-species transfer of grape-trained weights.
- OrthoDB eukaryote groups and gene-name proxies are **near-ortholog approximations**,
  not curated 1:1 orthologs.
- Do **not** treat conservation as experimental LLPS evidence.

## Reproduce

```bash
python3 scripts/analyze_fruit_cross_species.py
# optional threshold knobs:
python3 scripts/analyze_fruit_cross_species.py --min-conserved {min_conserved} --max-specific {max_specific}
```

This deposit: **{n_clusters}** scored clusters · conservative **{n_conserved}** ·
species-specific **{n_specific}**.
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fruit-dir", default=str(FRUIT))
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    ap.add_argument("--min-conserved", type=int, default=DEFAULT_MIN_CONSERVED)
    ap.add_argument("--max-specific", type=int, default=DEFAULT_MAX_SPECIFIC)
    ap.add_argument("--min-present-for-peak", type=int, default=DEFAULT_MIN_PRESENT_FOR_PEAK)
    args = ap.parse_args()

    fruit_dir = Path(args.fruit_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    index = load_index(fruit_dir / "index.json")
    slugs = scored_fruit_slugs(index, fruit_dir)
    if not slugs:
        raise SystemExit(f"No scored fruit-crop SQLite under {fruit_dir}")

    metas = {s: species_meta(index, s, fruit_dir) for s in slugs}
    clusters: dict[str, dict] = {}
    n_proteins = 0
    n_high = 0
    species_with_odb_ab: set[str] = set()

    for slug in slugs:
        rows = load_species_proteins(fruit_dir / f"{slug}_condensatome.sqlite", slug)
        n_proteins += len(rows)
        for r in rows:
            if r["tier"] in HIGH:
                n_high += 1
            key, key_type = cluster_key(
                r["orthodb_id"], r["grape_ortholog"], r["gene_name"], r["protein_id"], slug
            )
            if key_type == "orthodb" and r["tier"] in HIGH:
                species_with_odb_ab.add(slug)
            cl = clusters.get(key)
            if cl is None:
                cl = {
                    "cluster_id": key,
                    "key_type": key_type,
                    "orthodb_id": r["orthodb_id"].split(";")[0].strip() if r["orthodb_id"] else "",
                    "grape_ortholog": r["grape_ortholog"].split(";")[0].strip() if r["grape_ortholog"] else "",
                    "gene_symbol": normalize_gene(r["gene_name"]),
                    "by_species": defaultdict(list),
                }
                clusters[key] = cl
            else:
                # Prefer richer metadata when later rows carry OrthoDB/grape pointers
                if not cl["orthodb_id"] and r["orthodb_id"]:
                    cl["orthodb_id"] = r["orthodb_id"].split(";")[0].strip()
                if not cl["grape_ortholog"] and r["grape_ortholog"]:
                    cl["grape_ortholog"] = r["grape_ortholog"].split(";")[0].strip()
                if not cl["gene_symbol"]:
                    cl["gene_symbol"] = normalize_gene(r["gene_name"])
            cl["by_species"][slug].append(r)

    conserved_rows: list[dict] = []
    specific_rows: list[dict] = []
    member_rows: list[dict] = []

    for key, cl in clusters.items():
        by_sp = cl["by_species"]
        present_slugs = sorted(by_sp.keys())
        high_slugs = sorted(s for s, ms in by_sp.items() if any(m["tier"] in HIGH for m in ms))
        n_present = len(present_slugs)
        n_high_sp = len(high_slugs)
        if n_high_sp == 0:
            continue

        high_scores: list[float] = []
        best_by_species: dict[str, dict] = {}
        for s, ms in by_sp.items():
            bm = best_member(ms)
            best_by_species[s] = bm
            if bm["tier"] in HIGH and bm["score"] is not None:
                high_scores.append(bm["score"])

        median_high = statistics.median(high_scores) if high_scores else None
        mean_high = statistics.mean(high_scores) if high_scores else None
        max_high = max(high_scores) if high_scores else None

        # Representative label: prefer a gene symbol seen on a high member
        label_gene = ""
        label_name = ""
        for s in high_slugs:
            bm = best_by_species[s]
            if normalize_gene(bm["gene_name"]):
                label_gene = normalize_gene(bm["gene_name"])
                label_name = bm["protein_name"]
                break
        if not label_gene:
            label_gene = cl["gene_symbol"]
        if not label_name and high_slugs:
            label_name = best_by_species[high_slugs[0]]["protein_name"]

        species_scores = []
        species_tiers = []
        species_ids = []
        for s in slugs:
            if s not in best_by_species:
                species_scores.append("")
                species_tiers.append("")
                species_ids.append("")
                continue
            bm = best_by_species[s]
            species_scores.append("" if bm["score"] is None else f"{bm['score']:.6f}")
            species_tiers.append(bm["tier"])
            species_ids.append(bm["protein_id"])

        base = {
            "cluster_id": key,
            "key_type": cl["key_type"],
            "orthodb_id": cl["orthodb_id"],
            "grape_ortholog": cl["grape_ortholog"],
            "gene_symbol": label_gene,
            "protein_name_example": label_name,
            "n_species_present": n_present,
            "n_species_high": n_high_sp,
            "species_high": ";".join(high_slugs),
            "species_present": ";".join(present_slugs),
            "median_high_score": "" if median_high is None else round(median_high, 6),
            "mean_high_score": "" if mean_high is None else round(mean_high, 6),
            "max_high_score": "" if max_high is None else round(max_high, 6),
            "conservation_rank": "",  # filled after sort
            "peak_species": ";".join(high_slugs),
            "class": "",
        }
        # Wide matrix helpers for TSV consumers
        for i, s in enumerate(slugs):
            base[f"score__{s}"] = species_scores[i]
            base[f"tier__{s}"] = species_tiers[i]
            base[f"id__{s}"] = species_ids[i]

        is_conserved = n_high_sp >= args.min_conserved
        is_contrast_peak = (
            n_high_sp <= args.max_specific and n_present >= args.min_present_for_peak
        )
        is_orphan_high = n_high_sp <= args.max_specific and n_present < args.min_present_for_peak
        is_specific = is_contrast_peak or is_orphan_high

        if is_conserved:
            row = dict(base)
            row["class"] = "conservative"
            conserved_rows.append(row)
        if is_specific:
            row = dict(base)
            row["class"] = "panel_contrast" if is_contrast_peak else "orphan_high"
            # Prefer a single peak label when one species dominates by score
            if len(high_slugs) == 1:
                row["peak_species"] = high_slugs[0]
            else:
                # order by best score desc
                ordered = sorted(
                    high_slugs,
                    key=lambda s: best_by_species[s]["score"] or -1.0,
                    reverse=True,
                )
                row["peak_species"] = ";".join(ordered)
            specific_rows.append(row)

        if is_conserved or is_specific:
            for s, ms in by_sp.items():
                bm = best_member(ms)
                member_rows.append(
                    {
                        "cluster_id": key,
                        "class": "conservative"
                        if is_conserved
                        else ("panel_contrast" if is_contrast_peak else "orphan_high"),
                        "slug": s,
                        "protein_id": bm["protein_id"],
                        "gene_name": bm["gene_name"],
                        "protein_name": bm["protein_name"],
                        "tier": bm["tier"],
                        "score": "" if bm["score"] is None else round(bm["score"], 6),
                        "percentile": "" if bm["percentile"] is None else bm["percentile"],
                        "length": bm["length"] or "",
                        "is_high": "1" if bm["tier"] in HIGH else "0",
                        "grape_ortholog": bm["grape_ortholog"],
                        "orthodb_id": bm["orthodb_id"],
                        "key_type": cl["key_type"],
                    }
                )

    conserved_rows.sort(
        key=lambda r: (
            -int(r["n_species_high"]),
            -(float(r["median_high_score"]) if r["median_high_score"] != "" else -1.0),
            -(float(r["mean_high_score"]) if r["mean_high_score"] != "" else -1.0),
            r["cluster_id"],
        )
    )
    for i, r in enumerate(conserved_rows, 1):
        r["conservation_rank"] = i

    specific_rows.sort(
        key=lambda r: (
            r["peak_species"].split(";")[0],
            -int(r["n_species_high"] == 1),
            -(float(r["max_high_score"]) if r["max_high_score"] != "" else -1.0),
            r["cluster_id"],
        )
    )

    matrix_cols = []
    for s in slugs:
        matrix_cols.extend([f"score__{s}", f"tier__{s}", f"id__{s}"])

    conserved_cols = [
        "conservation_rank",
        "cluster_id",
        "key_type",
        "class",
        "orthodb_id",
        "grape_ortholog",
        "gene_symbol",
        "protein_name_example",
        "n_species_high",
        "n_species_present",
        "species_high",
        "species_present",
        "median_high_score",
        "mean_high_score",
        "max_high_score",
        "peak_species",
    ] + matrix_cols

    specific_cols = [
        "cluster_id",
        "key_type",
        "class",
        "peak_species",
        "orthodb_id",
        "grape_ortholog",
        "gene_symbol",
        "protein_name_example",
        "n_species_high",
        "n_species_present",
        "species_high",
        "species_present",
        "median_high_score",
        "mean_high_score",
        "max_high_score",
    ] + matrix_cols

    member_cols = [
        "cluster_id",
        "class",
        "key_type",
        "slug",
        "protein_id",
        "gene_name",
        "protein_name",
        "tier",
        "score",
        "percentile",
        "length",
        "is_high",
        "grape_ortholog",
        "orthodb_id",
    ]

    write_tsv_gz(out_dir / "conservative_candidates.tsv.gz", conserved_cols, conserved_rows)
    write_tsv_gz(out_dir / "species_specific_peaks.tsv.gz", specific_cols, specific_rows)
    write_tsv_gz(out_dir / "cluster_members.tsv.gz", member_cols, member_rows)

    # UI browse slice
    browse_conserved = []
    for r in conserved_rows[:BROWSE_CONSERVED_CAP]:
        high_set = set(r["species_high"].split(";")) if r["species_high"] else set()
        browse_conserved.append(
            {
                "conservation_rank": r["conservation_rank"],
                "cluster_id": r["cluster_id"],
                "key_type": r["key_type"],
                "gene_symbol": r["gene_symbol"],
                "protein_name_example": r["protein_name_example"],
                "orthodb_id": r["orthodb_id"],
                "grape_ortholog": r["grape_ortholog"],
                "n_species_high": r["n_species_high"],
                "n_species_present": r["n_species_present"],
                "species_high": r["species_high"].split(";") if r["species_high"] else [],
                "median_high_score": r["median_high_score"],
                "mean_high_score": r["mean_high_score"],
                "max_high_score": r["max_high_score"],
                "members": {
                    s: {
                        "protein_id": r.get(f"id__{s}") or None,
                        "tier": r.get(f"tier__{s}") or None,
                        "score": float(r[f"score__{s}"]) if r.get(f"score__{s}") else None,
                    }
                    for s in slugs
                    if s in high_set and r.get(f"id__{s}")
                },
            }
        )

    def specific_browse_priority(r: dict) -> tuple:
        """Prefer panel-contrast OrthoDB peaks, then Tier-A orphans, then score."""
        peak = (r["peak_species"] or "").split(";")[0]
        peak_tier = r.get(f"tier__{peak}") or ""
        class_rank = 0 if r["class"] == "panel_contrast" else 1
        key_rank = 0 if r["key_type"] == "orthodb" else (1 if r["key_type"] == "grape_ortholog" else 2)
        tier_rank = 0 if peak_tier == "A" else 1
        score = float(r["max_high_score"]) if r["max_high_score"] != "" else -1.0
        return (class_rank, key_rank, tier_rank, -score, r["cluster_id"])

    specific_sorted = sorted(specific_rows, key=specific_browse_priority)
    by_peak: dict[str, list] = defaultdict(list)
    for r in specific_sorted:
        peak = (r["peak_species"] or "").split(";")[0]
        if not peak:
            continue
        if len(by_peak[peak]) >= BROWSE_SPECIFIC_PER_SPECIES:
            continue
        # Compact members: peak species + other high species only (UI filter still has species list).
        member_slugs = set(r["species_high"].split(";")) if r["species_high"] else set()
        member_slugs.add(peak)
        by_peak[peak].append(
            {
                "cluster_id": r["cluster_id"],
                "key_type": r["key_type"],
                "class": r["class"],
                "gene_symbol": r["gene_symbol"],
                "protein_name_example": r["protein_name_example"],
                "orthodb_id": r["orthodb_id"],
                "grape_ortholog": r["grape_ortholog"],
                "peak_species": r["peak_species"].split(";") if r["peak_species"] else [],
                "n_species_high": r["n_species_high"],
                "n_species_present": r["n_species_present"],
                "species_high": r["species_high"].split(";") if r["species_high"] else [],
                "max_high_score": r["max_high_score"],
                "median_high_score": r["median_high_score"],
                "members": {
                    s: {
                        "protein_id": r.get(f"id__{s}") or None,
                        "tier": r.get(f"tier__{s}") or None,
                        "score": float(r[f"score__{s}"]) if r.get(f"score__{s}") else None,
                    }
                    for s in sorted(member_slugs)
                    if r.get(f"id__{s}")
                },
            }
        )

    browse = {
        "version": 1,
        "generated_at_utc": utc_now(),
        "thresholds": {
            "high_tiers": sorted(HIGH),
            "min_species_conserved": args.min_conserved,
            "max_species_specific": args.max_specific,
            "min_present_for_panel_contrast": args.min_present_for_peak,
        },
        "species": [metas[s] for s in slugs],
        "conservative": browse_conserved,
        "species_specific_by_peak": {s: by_peak.get(s, []) for s in slugs},
        "n_conservative_total": len(conserved_rows),
        "n_species_specific_total": len(specific_rows),
        "honesty": (
            "Near-ortholog approximation via OrthoDB / grape join / gene-name proxy. "
            "Grape-trained PlantLLPS tiers transferred across fruit proteomes — not experimental LLPS."
        ),
    }
    (out_dir / "browse.json").write_text(json.dumps(browse, indent=2) + "\n", encoding="utf-8")

    n_clusters_scored = sum(
        1
        for cl in clusters.values()
        if any(any(m["tier"] in HIGH for m in ms) for ms in cl["by_species"].values())
    )
    key_type_counts = defaultdict(int)
    for r in conserved_rows:
        key_type_counts[f"conserved:{r['key_type']}"] += 1
    for r in specific_rows:
        key_type_counts[f"specific:{r['key_type']}"] += 1
        key_type_counts[f"specific_class:{r['class']}"] += 1

    summary = {
        "version": 1,
        "generated_at_utc": utc_now(),
        "tool": "scripts/analyze_fruit_cross_species.py",
        "n_panel_species": len(slugs),
        "species": [metas[s] for s in slugs],
        "n_proteins_scanned": n_proteins,
        "n_tier_ab_proteins": n_high,
        "n_clusters_with_high": n_clusters_scored,
        "n_conservative": len(conserved_rows),
        "n_species_specific": len(specific_rows),
        "n_species_specific_panel_contrast": sum(
            1 for r in specific_rows if r["class"] == "panel_contrast"
        ),
        "n_species_specific_orphan_high": sum(
            1 for r in specific_rows if r["class"] == "orphan_high"
        ),
        "n_species_with_orthodb_ab": len(species_with_odb_ab),
        "species_with_orthodb_ab": sorted(species_with_odb_ab),
        "thresholds": {
            "high_tiers": sorted(HIGH),
            "min_species_conserved": args.min_conserved,
            "max_species_specific": args.max_specific,
            "min_present_for_panel_contrast": args.min_present_for_peak,
        },
        "key_type_counts": dict(sorted(key_type_counts.items())),
        "files": {
            "conservative_candidates": "conservative_candidates.tsv.gz",
            "species_specific_peaks": "species_specific_peaks.tsv.gz",
            "cluster_members": "cluster_members.tsv.gz",
            "browse": "browse.json",
            "readme": "README.md",
        },
        "honesty": (
            "Cross-species transfer of grape-trained PlantLLPS; OrthoDB eukaryote groups "
            "and gene-name proxies are near-ortholog approximations. Do not treat ranks as "
            "experimental LLPS calls. Species lacking OrthoDB xrefs enrich orphan-high peaks."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (out_dir / "README.md").write_text(
        methods_readme(
            n_species=len(slugs),
            min_conserved=args.min_conserved,
            max_specific=args.max_specific,
            min_present=args.min_present_for_peak,
            n_clusters=n_clusters_scored,
            n_conserved=len(conserved_rows),
            n_specific=len(specific_rows),
            n_with_odb_species=len(species_with_odb_ab),
        ),
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "ok": True,
                "out_dir": str(out_dir),
                "n_panel_species": len(slugs),
                "n_conservative": len(conserved_rows),
                "n_species_specific": len(specific_rows),
                "n_clusters_with_high": n_clusters_scored,
                "species_with_orthodb_ab": sorted(species_with_odb_ab),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
