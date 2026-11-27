#!/usr/bin/env python3
"""Cross-crop conserved condensatome: grape ↔ soybean ↔ Arabidopsis Tier A/B.

Joins grape atlas Tier A∪B proteins to soybean (OrthoDB + deposited scores) and
Arabidopsis (OrthoDB accession + optional deposited Ath PlantLLPS tiers).

Deposits under data/cross_crop/:
  - conserved_tier_ab_overlap.tsv.gz   per-grape-protein overlap table
  - conserved_orthodb_clusters.tsv.gz  OrthoDB group-level summary
  - conserved_summary.json             counts + honesty

Example::

    python3 scripts/analyze_cross_crop_conserved.py
    python3 scripts/score_second_crop_proteome.py --crop arabidopsis --out-dir data/cross_crop
    python3 scripts/analyze_cross_crop_conserved.py  # refresh with Ath tiers
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sqlite3
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DEFAULT_GRAPE = DATA / "grape_condensatome.sqlite"
DEFAULT_SOY = DATA / "second_crop" / "soybean_condensatome.sqlite"
DEFAULT_ATH = DATA / "cross_crop" / "arabidopsis_condensatome.sqlite"
OUT_DIR = DATA / "cross_crop"
HOOKS = DATA / "cross_crop_hooks.json"

HIGH = frozenset({"A", "B"})


def load_scores(db_path: Path) -> dict[str, dict]:
    if not db_path.exists():
        return {}
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    out = {
        r["protein_id"]: {
            "score": r["score"],
            "tier": r["tier"],
            "length": r["length"] if "length" in r.keys() else None,
            "percentile": r["percentile"] if "percentile" in r.keys() else None,
        }
        for r in conn.execute(
            "SELECT protein_id, score, tier, length, percentile FROM proteins"
        )
    }
    conn.close()
    return out


def load_hooks(path: Path) -> dict:
    if not path.exists():
        return {"hooks": []}
    return json.loads(path.read_text())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grape-db", default=str(DEFAULT_GRAPE))
    ap.add_argument("--soy-db", default=str(DEFAULT_SOY))
    ap.add_argument("--ath-db", default=str(DEFAULT_ATH))
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    ap.add_argument("--hooks", default=str(HOOKS))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    soy_scores = load_scores(Path(args.soy_db))
    ath_scores = load_scores(Path(args.ath_db))
    hooks = load_hooks(Path(args.hooks))

    hook_by_grape: dict[str, str] = {}
    for h in hooks.get("hooks") or []:
        fam = h.get("family") or ""
        for gid in h.get("grape") or []:
            hook_by_grape[gid] = fam

    conn = sqlite3.connect(str(args.grape_db))
    conn.row_factory = sqlite3.Row
    grape_rows = list(
        conn.execute(
            """
            SELECT protein_id, gene_name, protein_name, tier, score, length,
                   orthodb_id, ortholog_ath, ortholog_ath_gene,
                   ortholog_gma, ortholog_gma_gene, ortholog_gma_score, ortholog_gma_tier,
                   cross_crop_family, ortholog_json
            FROM proteins
            """
        )
    )
    n_grape = len(grape_rows)
    n_grape_ab = sum(1 for r in grape_rows if r["tier"] in HIGH)
    conn.close()

    overlap_rows: list[dict] = []
    cluster_members: dict[str, dict] = defaultdict(
        lambda: {
            "grape_ab": [],
            "grape_all": [],
            "soy_ids": set(),
            "soy_ab": set(),
            "ath_ids": set(),
            "ath_ab": set(),
        }
    )

    n_grape_ab_soy_ab = 0
    n_grape_ab_ath_link = 0
    n_grape_ab_ath_ab = 0
    n_triple_ab = 0  # grape AB ∩ soy AB ∩ ath AB
    n_triple_odb = 0  # grape AB ∩ soy AB ∩ ath OrthoDB present
    n_grape_specific_ab = 0  # grape AB with neither soy AB nor ath AB (or ath link)

    for r in grape_rows:
        pid = r["protein_id"]
        tier = r["tier"]
        odb = r["orthodb_id"] or ""
        gma = r["ortholog_gma"]
        gma_tier = r["ortholog_gma_tier"]
        ath = r["ortholog_ath"]
        ath_gene = r["ortholog_ath_gene"] or ""

        # Prefer live soy DB tier when available
        soy_tier = gma_tier
        soy_score = r["ortholog_gma_score"]
        if gma and gma in soy_scores:
            soy_tier = soy_scores[gma]["tier"]
            soy_score = soy_scores[gma]["score"]

        ath_tier = None
        ath_score = None
        if ath and ath in ath_scores:
            ath_tier = ath_scores[ath]["tier"]
            ath_score = ath_scores[ath]["score"]

        soy_ab = bool(gma and soy_tier in HIGH)
        ath_link = bool(ath)
        ath_ab = bool(ath and ath_tier in HIGH)

        if odb:
            cm = cluster_members[odb]
            cm["grape_all"].append(pid)
            if tier in HIGH:
                cm["grape_ab"].append(pid)
            if gma:
                cm["soy_ids"].add(gma)
                if soy_ab:
                    cm["soy_ab"].add(gma)
            if ath:
                cm["ath_ids"].add(ath)
                if ath_ab:
                    cm["ath_ab"].add(ath)

        if tier not in HIGH:
            continue

        class_label = "grape_specific"
        if soy_ab and ath_ab:
            class_label = "triple_tier_ab"
            n_triple_ab += 1
        elif soy_ab and ath_link:
            class_label = "grape_soy_ab_ath_ortholog"
            n_triple_odb += 1
        elif soy_ab:
            class_label = "grape_soy_ab"
        elif ath_ab:
            class_label = "grape_ath_ab"
        elif ath_link:
            class_label = "grape_ab_ath_ortholog_only"
        else:
            n_grape_specific_ab += 1

        if soy_ab:
            n_grape_ab_soy_ab += 1
        if ath_link:
            n_grape_ab_ath_link += 1
        if ath_ab:
            n_grape_ab_ath_ab += 1

        fam = r["cross_crop_family"] or hook_by_grape.get(pid) or ""
        overlap_rows.append(
            {
                "grape_id": pid,
                "grape_gene": r["gene_name"] or "",
                "grape_tier": tier,
                "grape_score": round(float(r["score"]), 6) if r["score"] is not None else "",
                "orthodb_id": odb,
                "soybean_id": gma or "",
                "soybean_gene": r["ortholog_gma_gene"] or "",
                "soybean_tier": soy_tier or "",
                "soybean_score": round(float(soy_score), 6) if soy_score is not None else "",
                "ath_id": ath or "",
                "ath_gene": ath_gene,
                "ath_tier": ath_tier or "",
                "ath_score": round(float(ath_score), 6) if ath_score is not None else "",
                "overlap_class": class_label,
                "cross_crop_family": fam,
            }
        )

    # Cluster-level: OrthoDB groups with ≥1 grape Tier A/B
    cluster_rows: list[dict] = []
    n_clusters_grape_ab = 0
    n_clusters_soy_ab = 0
    n_clusters_ath_ab = 0
    n_clusters_triple_ab = 0
    n_clusters_triple_odb = 0

    for oid, cm in sorted(cluster_members.items()):
        if not cm["grape_ab"]:
            continue
        n_clusters_grape_ab += 1
        has_soy_ab = bool(cm["soy_ab"])
        has_ath_ab = bool(cm["ath_ab"])
        has_ath = bool(cm["ath_ids"])
        if has_soy_ab:
            n_clusters_soy_ab += 1
        if has_ath_ab:
            n_clusters_ath_ab += 1
        if has_soy_ab and has_ath_ab:
            n_clusters_triple_ab += 1
        if has_soy_ab and has_ath:
            n_clusters_triple_odb += 1
        if has_soy_ab and has_ath_ab:
            cclass = "triple_tier_ab"
        elif has_soy_ab and has_ath:
            cclass = "grape_soy_ab_ath_ortholog"
        elif has_soy_ab:
            cclass = "grape_soy_ab"
        elif has_ath_ab:
            cclass = "grape_ath_ab"
        elif has_ath:
            cclass = "grape_ab_ath_ortholog_only"
        else:
            cclass = "grape_specific"
        cluster_rows.append(
            {
                "orthodb_id": oid,
                "n_grape_ab": len(cm["grape_ab"]),
                "n_grape_all": len(cm["grape_all"]),
                "n_soy": len(cm["soy_ids"]),
                "n_soy_ab": len(cm["soy_ab"]),
                "n_ath": len(cm["ath_ids"]),
                "n_ath_ab": len(cm["ath_ab"]),
                "grape_ab_ids": ";".join(cm["grape_ab"][:20]),
                "soy_ab_ids": ";".join(sorted(cm["soy_ab"])[:20]),
                "ath_ab_ids": ";".join(sorted(cm["ath_ab"])[:20]),
                "cluster_class": cclass,
            }
        )

    # Fix grape-specific count: proteins labeled grape_specific among Tier A/B
    n_grape_specific_ab = sum(1 for r in overlap_rows if r["overlap_class"] == "grape_specific")
    # n_triple_odb counted only grape_soy_ab_ath_ortholog; also include triple_tier_ab in "shared odb"
    n_shared_soy_ath_odb = sum(
        1
        for r in overlap_rows
        if r["overlap_class"] in ("triple_tier_ab", "grape_soy_ab_ath_ortholog")
    )

    overlap_path = out_dir / "conserved_tier_ab_overlap.tsv.gz"
    overlap_cols = [
        "grape_id",
        "grape_gene",
        "grape_tier",
        "grape_score",
        "orthodb_id",
        "soybean_id",
        "soybean_gene",
        "soybean_tier",
        "soybean_score",
        "ath_id",
        "ath_gene",
        "ath_tier",
        "ath_score",
        "overlap_class",
        "cross_crop_family",
    ]
    with gzip.open(overlap_path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=overlap_cols, delimiter="\t")
        w.writeheader()
        for row in sorted(overlap_rows, key=lambda x: (-float(x["grape_score"] or 0), x["grape_id"])):
            w.writerow(row)

    cluster_path = out_dir / "conserved_orthodb_clusters.tsv.gz"
    cluster_cols = [
        "orthodb_id",
        "n_grape_ab",
        "n_grape_all",
        "n_soy",
        "n_soy_ab",
        "n_ath",
        "n_ath_ab",
        "grape_ab_ids",
        "soy_ab_ids",
        "ath_ab_ids",
        "cluster_class",
    ]
    with gzip.open(cluster_path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cluster_cols, delimiter="\t")
        w.writeheader()
        for row in sorted(cluster_rows, key=lambda x: (-x["n_grape_ab"], x["orthodb_id"])):
            w.writerow(row)

    class_counts = defaultdict(int)
    for r in overlap_rows:
        class_counts[r["overlap_class"]] += 1

    honesty = (
        "Conserved condensatome overlap uses grape PlantLLPS tiers (atlas SQLite), "
        "soybean OrthoDB joins + deposited soybean tiers (data/second_crop/), and "
        "Arabidopsis OrthoDB accessions on grape rows. When arabidopsis_condensatome.sqlite "
        "is present, Ath Tier A/B is assigned with grape-trained PlantLLPS v13 "
        "(composition-proxy disorder; tiers recalibrated on Ath). OrthoDB eukaryote-level "
        "groups can be broad; curated cross_crop_hooks.json overlays are noted via "
        "cross_crop_family when present. Not a multi-species retrained model."
    )

    summary = {
        "version": 1,
        "grape_db": Path(args.grape_db).name,
        "soy_db": Path(args.soy_db).name if Path(args.soy_db).exists() else None,
        "ath_db": Path(args.ath_db).name if Path(args.ath_db).exists() else None,
        "n_grape_proteins": n_grape,
        "n_grape_tier_ab": n_grape_ab,
        "n_soy_scored": len(soy_scores),
        "n_ath_scored": len(ath_scores),
        "protein_overlap": {
            "grape_ab_with_soy_ab": n_grape_ab_soy_ab,
            "grape_ab_with_ath_ortholog": n_grape_ab_ath_link,
            "grape_ab_with_ath_ab": n_grape_ab_ath_ab,
            "triple_tier_ab": n_triple_ab,
            "grape_soy_ab_and_ath_ortholog": n_shared_soy_ath_odb,
            "grape_specific_ab": n_grape_specific_ab,
            "by_class": dict(sorted(class_counts.items())),
        },
        "orthodb_clusters": {
            "n_with_grape_ab": n_clusters_grape_ab,
            "n_with_soy_ab": n_clusters_soy_ab,
            "n_with_ath_ab": n_clusters_ath_ab,
            "n_triple_tier_ab": n_clusters_triple_ab,
            "n_grape_soy_ab_ath_ortholog": n_clusters_triple_odb,
        },
        "deposits": {
            "overlap_tsv": overlap_path.name,
            "clusters_tsv": cluster_path.name,
        },
        "honesty": honesty,
    }
    summary_path = out_dir / "conserved_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    # Also stamp grape meta (thin pointer; no large tables in SQLite)
    gconn = sqlite3.connect(str(args.grape_db))
    gconn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("cross_crop_conserved", json.dumps(summary)),
    )
    gconn.commit()
    gconn.close()

    print(json.dumps(summary, indent=2), flush=True)
    print(f"Overlap → {overlap_path} ({overlap_path.stat().st_size} bytes)", flush=True)
    print(f"Clusters → {cluster_path} ({cluster_path.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    main()
