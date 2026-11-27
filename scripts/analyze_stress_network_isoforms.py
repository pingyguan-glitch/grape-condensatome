#!/usr/bin/env python3
"""Atlas analyses 2 / 6 / 7: stress deep-dive, STRING/stress hubs, isoform deltas.

Builds on existing stress_modules, STRING neighborhood, and UniProt GO/keyword joins
already deposited in the grape condensatome (no new network downloads required).

1. **Stress module deep dive (analysis 2)**
   - Primary focus: FDR-significant Botrytis down contrast
     (E-GEOD-67932 ``g1_g4`` down ∩ Tier A∪B).
   - UniProt GO-ID + keyword hypergeometric enrichment (KEGG-style pathway
     annotation is not in the DB; keywords stand in for functional classes).
   - Heatmap-ready protein × contrast log2FC matrix for the primary set.
   - Other contrast×direction sets are run and labeled **exploratory**.

2. **STRING / stress seed hubs (analysis 6)**
   - Rank Tier A∪B STRING 1-hop neighborhood enrichment (min neighbor count).
   - Rank stress log2FC coexpression seed hubs (Tier A seeds from modules pass).
   - Deposit a combined top-hubs TSV.

3. **Isoform / multi-protein genes (analysis 7)**
   - Group scored UniProt accessions by shared ``VIT_*`` locus ID.
   - Deposit high score-delta multi-accession genes (honest: UP000009183 is
     nearly 1:1 VIT↔accession).

Example::

    python3 scripts/analyze_stress_network_isoforms.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
STRESS_SIDECAR = ROOT / "data" / "expression_atlas" / "stress_by_protein.json.gz"
STRESS_ENRICH = ROOT / "data" / "stress_modules" / "de_enrichment.tsv"
STRESS_MEMBERS = ROOT / "data" / "stress_modules" / "module_members.tsv.gz"
STRESS_SEEDS = ROOT / "data" / "stress_modules" / "seed_neighborhoods.tsv.gz"
STRING_STATS = ROOT / "data" / "string" / "neighborhood_stats.tsv.gz"

DEEPDIVE_DIR = ROOT / "data" / "stress_modules" / "deepdive"
HUBS_TSV = ROOT / "data" / "string" / "stress_seed_hubs.tsv"
ISOFORM_DIR = ROOT / "data" / "isoforms"

P_THR = 0.05
LFC_THR = 1.0
VIT_RE = re.compile(r"VIT_[0-9]{2}s[0-9]{4}g[0-9]{5}", re.IGNORECASE)

# Skip ubiquitous / non-informative UniProt keywords for enrichment.
KW_SKIP = {
    "Reference proteome",
    "3D-structure",
    "Direct protein sequencing",
    "Complete proteome",
}


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(N, K, n)."""
    if n == 0 or K == 0 or N == 0:
        return 1.0
    k = max(0, min(k, n, K))
    try:
        log_denom = math.lgamma(N + 1) - math.lgamma(n + 1) - math.lgamma(N - n + 1)
    except ValueError:
        return 1.0
    total = 0.0
    lo = max(0, n - (N - K))
    for i in range(max(k, lo), min(n, K) + 1):
        log_p = (
            math.lgamma(K + 1)
            - math.lgamma(i + 1)
            - math.lgamma(K - i + 1)
            + math.lgamma(N - K + 1)
            - math.lgamma(n - i + 1)
            - math.lgamma(N - K - (n - i) + 1)
            - log_denom
        )
        total += math.exp(log_p)
    if total <= 0:
        return 1.0
    return min(1.0, total)


def bh_fdr(pvals: list[float]) -> list[float]:
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        idx = order[rank - 1]
        val = min(prev, pvals[idx] * m / rank)
        adj[idx] = min(1.0, val)
        prev = adj[idx]
    return adj


def load_json_gz(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def contrast_axis(stress: dict[str, dict]) -> list[tuple[str, str, str]]:
    keys: dict[tuple[str, str], str] = {}
    for blob in stress.values():
        for acc, info in blob.items():
            for lab, comp in (info.get("comparisons") or {}).items():
                cid = str(comp.get("contrast_id") or lab)
                keys.setdefault((acc, lab), cid)
    out = [(acc, lab, cid) for (acc, lab), cid in keys.items()]
    out.sort(key=lambda t: (t[0], t[2], t[1]))
    return out


def parse_terms(raw: str | None, sep: str = ";") -> list[str]:
    if not raw:
        return []
    return [t.strip() for t in raw.split(sep) if t.strip()]


def term_enrichment(
    foreground: set[str],
    background: set[str],
    term_map: dict[str, list[str]],
    *,
    kind: str,
    min_term_count: int,
    min_overlap: int,
) -> list[dict]:
    """Hypergeometric enrichment of annotation terms in foreground vs background."""
    bg_term_proteins: dict[str, set[str]] = defaultdict(set)
    for pid in background:
        for term in term_map.get(pid, []):
            bg_term_proteins[term].add(pid)
    N = len(background)
    n = len(foreground)
    rows: list[dict] = []
    for term, bg_set in bg_term_proteins.items():
        K = len(bg_set)
        if K < min_term_count:
            continue
        overlap = foreground & bg_set
        k = len(overlap)
        if k < min_overlap:
            continue
        bg_rate = K / N if N else 0.0
        fold = (k / n / bg_rate) if n and bg_rate else None
        pval = hypergeom_sf(k, N, K, n) if n else 1.0
        rows.append(
            {
                "term_kind": kind,
                "term": term,
                "n_background_with_term": K,
                "n_foreground": n,
                "n_overlap": k,
                "background_rate": round(bg_rate, 6),
                "foreground_rate": round(k / n, 6) if n else None,
                "enrichment_fold": round(fold, 6) if fold is not None else None,
                "hypergeom_pvalue": pval,
                "overlap_proteins": ",".join(sorted(overlap)[:25]),
            }
        )
    pvals = [r["hypergeom_pvalue"] for r in rows]
    for r, fdr in zip(rows, bh_fdr(pvals)):
        r["bh_fdr"] = round(fdr, 6)
        r["hypergeom_pvalue"] = round(r["hypergeom_pvalue"], 6)
    rows.sort(key=lambda r: (r["bh_fdr"], r["hypergeom_pvalue"], -(r["enrichment_fold"] or 0)))
    return rows


def load_set_members(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["set_id"]].append(row["protein_id"])
    return out


def analysis_stress_deepdive(
    con: sqlite3.Connection,
    stress: dict[str, dict],
    *,
    out_dir: Path,
    min_term_count: int,
    min_overlap: int,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = list(
        con.execute(
            "SELECT protein_id, gene_name, tier, score, go_terms, keywords FROM proteins"
        )
    )
    tiers = {r["protein_id"]: (r["tier"] or "") for r in rows}
    scores = {r["protein_id"]: float(r["score"] or 0.0) for r in rows}
    genes = {r["protein_id"]: (r["gene_name"] or "") for r in rows}
    go_map = {r["protein_id"]: parse_terms(r["go_terms"]) for r in rows}
    kw_map = {
        r["protein_id"]: [k for k in parse_terms(r["keywords"]) if k not in KW_SKIP]
        for r in rows
    }

    atlas_ids = set(tiers)
    bg = set(stress.keys()) & atlas_ids
    tier_ab = {pid for pid, t in tiers.items() if t in ("A", "B")}
    bg_ab = bg & tier_ab

    enrich_tsv = list(csv.DictReader(STRESS_ENRICH.open(encoding="utf-8"), delimiter="\t"))
    primary = [
        r
        for r in enrich_tsv
        if r["accession"] == "E-GEOD-67932"
        and r["contrast_id"] == "g1_g4"
        and r["direction"] == "down"
        and float(r["bh_fdr"]) < 0.05
    ]
    if not primary:
        raise SystemExit(
            "Expected FDR-significant E-GEOD-67932 g1_g4 down enrichment row missing"
        )
    primary_row = primary[0]
    primary_set_id = primary_row["gene_set_id"]

    members = load_set_members(STRESS_MEMBERS)
    if primary_set_id not in members:
        raise SystemExit(f"Gene set not in module_members: {primary_set_id}")

    # Build DE∩TierAB sets for all contrasts (reuse deposited membership when present)
    axis = contrast_axis(stress)
    set_specs: list[dict] = []
    for acc, lab, cid in axis:
        for direction in ("up", "down", "any"):
            set_id = f"DE_{acc}_{cid}_{direction}_TierAB"
            is_primary = set_id == primary_set_id
            if set_id in members:
                fg = set(members[set_id])
            else:
                fg = set()
                for pid in bg_ab:
                    blob = stress.get(pid) or {}
                    info = blob.get(acc) or {}
                    comp = (info.get("comparisons") or {}).get(lab)
                    if not comp or not comp.get("significant"):
                        continue
                    lfc = float(comp.get("log2fc") or 0.0)
                    if direction == "any":
                        fg.add(pid)
                    elif direction == "up" and lfc >= LFC_THR:
                        fg.add(pid)
                    elif direction == "down" and lfc <= -LFC_THR:
                        fg.add(pid)
            de_fdr = next(
                (
                    float(r["bh_fdr"])
                    for r in enrich_tsv
                    if r["gene_set_id"] == set_id
                ),
                None,
            )
            set_specs.append(
                {
                    "set_id": set_id,
                    "accession": acc,
                    "contrast_id": cid,
                    "contrast_label": lab,
                    "direction": direction,
                    "scope": "primary" if is_primary else "exploratory",
                    "de_tier_ab_fdr": de_fdr,
                    "members": sorted(fg),
                }
            )

    enrich_rows: list[dict] = []
    for spec in set_specs:
        fg = set(spec["members"])
        if len(fg) < 5:
            continue
        for kind, tmap in (("go", go_map), ("keyword", kw_map)):
            for row in term_enrichment(
                fg,
                bg,
                tmap,
                kind=kind,
                min_term_count=min_term_count,
                min_overlap=min_overlap,
            ):
                enrich_rows.append(
                    {
                        "scope": spec["scope"],
                        "set_id": spec["set_id"],
                        "accession": spec["accession"],
                        "contrast_id": spec["contrast_id"],
                        "contrast_label": spec["contrast_label"],
                        "direction": spec["direction"],
                        "de_tier_ab_fdr": spec["de_tier_ab_fdr"],
                        **row,
                    }
                )

    enrich_path = out_dir / "go_keyword_enrichment.tsv"
    enrich_cols = [
        "scope",
        "set_id",
        "accession",
        "contrast_id",
        "contrast_label",
        "direction",
        "de_tier_ab_fdr",
        "term_kind",
        "term",
        "n_background_with_term",
        "n_foreground",
        "n_overlap",
        "background_rate",
        "foreground_rate",
        "enrichment_fold",
        "hypergeom_pvalue",
        "bh_fdr",
        "overlap_proteins",
    ]
    with enrich_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=enrich_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in enrich_rows:
            w.writerow(r)

    # Heatmap-ready matrix for primary set: proteins × contrasts (log2FC)
    primary_members = members[primary_set_id]
    col_keys = [(acc, lab, cid) for acc, lab, cid in axis]
    heat_path = out_dir / "heatmap_matrix.tsv.gz"
    with gzip.open(heat_path, "wt", encoding="utf-8", newline="") as fh:
        header = [
            "protein_id",
            "gene_name",
            "tier",
            "score",
            *[f"{acc}|{cid}|log2fc" for acc, _lab, cid in col_keys],
            *[f"{acc}|{cid}|significant" for acc, _lab, cid in col_keys],
        ]
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for pid in sorted(primary_members, key=lambda p: (-scores.get(p, 0.0), p)):
            blob = stress.get(pid) or {}
            lfc_vals: list[str] = []
            sig_vals: list[str] = []
            for acc, lab, _cid in col_keys:
                info = blob.get(acc) or {}
                comp = (info.get("comparisons") or {}).get(lab)
                if not comp or comp.get("log2fc") is None:
                    lfc_vals.append("")
                    sig_vals.append("")
                else:
                    lfc_vals.append(f"{float(comp['log2fc']):.4f}")
                    sig_vals.append("1" if comp.get("significant") else "0")
            w.writerow(
                [
                    pid,
                    genes.get(pid) or "",
                    tiers.get(pid) or "",
                    f"{scores.get(pid, 0.0):.6f}",
                    *lfc_vals,
                    *sig_vals,
                ]
            )

    primary_enrich = [r for r in enrich_rows if r["scope"] == "primary"]
    primary_sig = [r for r in primary_enrich if r["bh_fdr"] < 0.05]
    top_primary = primary_enrich[:15]

    summary = {
        "layer": "stress_module_deepdive",
        "version": "1.5.0",
        "primary_set_id": primary_set_id,
        "primary_contrast": {
            "accession": primary_row["accession"],
            "contrast_id": primary_row["contrast_id"],
            "contrast_label": primary_row["contrast_label"],
            "direction": primary_row["direction"],
            "n_de_tier_ab": int(primary_row["n_de_tier_ab"]),
            "enrichment_fold": float(primary_row["enrichment_fold"]),
            "bh_fdr": float(primary_row["bh_fdr"]),
        },
        "n_primary_members": len(primary_members),
        "n_enrichment_rows": len(enrich_rows),
        "n_primary_term_tests": len(primary_enrich),
        "n_primary_term_fdr_lt_0.05": len(primary_sig),
        "n_exploratory_term_tests": sum(1 for r in enrich_rows if r["scope"] == "exploratory"),
        "annotation_note": (
            "GO IDs and UniProt keywords from proteins.go_terms / keywords. "
            "No KEGG pathway table is deposited in this atlas; keyword enrichment "
            "is the KEGG-style functional-class proxy."
        ),
        "top_primary_terms": top_primary,
        "paths": {
            "go_keyword_enrichment": str(enrich_path.relative_to(ROOT)),
            "heatmap_matrix": str(heat_path.relative_to(ROOT)),
        },
        "honesty": (
            "Primary scope = FDR-significant Botrytis stage-3 down ∩ Tier A∪B only. "
            "All other contrast×direction GO/keyword tests are exploratory and must "
            "not be cited as confirmed pathway hits. Enrichment is descriptive "
            "hypergeometric overlap on UniProt annotations — not g:Profiler/clusterProfiler "
            "with ontology propagation, and not KEGG pathway analysis."
        ),
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    (out_dir / "README.md").write_text(
        f"""# Stress module deep dive (analysis 2)

GO-ID + UniProt **keyword** enrichment for condensatome × stress DE gene sets, plus a
heatmap-ready log2FC matrix for the FDR-significant Botrytis down set.

## Primary set

- **{primary_set_id}** — E-GEOD-67932 stage-3 Botrytis down ∩ Tier A∪B
- Members: **{len(primary_members):,}**
- DE∩TierAB FDR: **{primary_row['bh_fdr']}** (fold ≈ {primary_row['enrichment_fold']})

## Deposits

| File | Contents |
|------|----------|
| `go_keyword_enrichment.tsv` | GO + keyword hypergeometric tests; `scope=primary\\|exploratory` |
| `heatmap_matrix.tsv.gz` | Primary members × contrast log2FC (+ significance flags) |
| `summary.json` | Run summary + honesty |

## Reproduce

```bash
python3 scripts/analyze_stress_network_isoforms.py
```

## Honesty

{summary["honesty"]}

{summary["annotation_note"]}
""",
        encoding="utf-8",
    )
    print(
        f"Deep dive: primary members={len(primary_members):,} "
        f"term rows={len(enrich_rows):,} primary FDR<0.05 terms={len(primary_sig):,}",
        flush=True,
    )
    return summary


def analysis_hubs(
    con: sqlite3.Connection,
    *,
    out_path: Path,
    top_n: int,
    min_string_neighbors: int,
) -> dict:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(con.execute("SELECT protein_id, gene_name, tier, score FROM proteins"))
    genes = {r["protein_id"]: (r["gene_name"] or "") for r in rows}
    tiers = {r["protein_id"]: (r["tier"] or "") for r in rows}
    scores = {r["protein_id"]: float(r["score"] or 0.0) for r in rows}
    tier_ab = {pid for pid, t in tiers.items() if t in ("A", "B")}

    string_rows: list[dict] = []
    with gzip.open(STRING_STATS, "rt", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["protein_id"] not in tier_ab:
                continue
            n = int(r["n_neighbors"] or 0)
            if n < min_string_neighbors:
                continue
            fold = float(r["enrichment_fold"]) if r["enrichment_fold"] else None
            pval = float(r["hypergeom_pvalue"]) if r["hypergeom_pvalue"] else None
            string_rows.append(
                {
                    "protein_id": r["protein_id"],
                    "string_n_neighbors": n,
                    "string_n_tier_ab": int(r["n_tier_ab"] or 0),
                    "string_frac_tier_ab": float(r["frac_tier_ab"] or 0),
                    "string_enrichment_fold": fold,
                    "string_hypergeom_pvalue": pval,
                    "string_mean_neighbor_score": (
                        float(r["mean_neighbor_score"])
                        if r["mean_neighbor_score"]
                        else None
                    ),
                }
            )
    string_rows.sort(
        key=lambda r: (
            -(r["string_enrichment_fold"] or 0),
            r["string_hypergeom_pvalue"] if r["string_hypergeom_pvalue"] is not None else 1,
            -r["string_n_neighbors"],
            -scores.get(r["protein_id"], 0.0),
        )
    )
    for i, r in enumerate(string_rows, start=1):
        r["string_rank"] = i

    stress_rows: list[dict] = []
    with gzip.open(STRESS_SEEDS, "rt", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            stress_rows.append(
                {
                    "protein_id": r["protein_id"],
                    "stress_n_partners_r_ge_0.7": int(r["n_partners_r_ge_0.7"] or 0),
                    "stress_n_partners_r_ge_thr": int(r["n_partners_r_ge_thr"] or 0),
                    "stress_top_partner_id": r["top_partner_id"],
                    "stress_top_partner_r": float(r["top_partner_r"] or 0),
                }
            )
    stress_rows.sort(
        key=lambda r: (
            -r["stress_n_partners_r_ge_0.7"],
            -r["stress_top_partner_r"],
            -scores.get(r["protein_id"], 0.0),
        )
    )
    for i, r in enumerate(stress_rows, start=1):
        r["stress_rank"] = i

    by_id: dict[str, dict] = {}
    for r in string_rows:
        by_id[r["protein_id"]] = {**r}
    for r in stress_rows:
        by_id.setdefault(r["protein_id"], {"protein_id": r["protein_id"]}).update(r)

    combined: list[dict] = []
    for pid, r in by_id.items():
        s_rank = r.get("string_rank")
        t_rank = r.get("stress_rank")
        sources = []
        if s_rank is not None:
            sources.append("string_neighborhood")
        if t_rank is not None:
            sources.append("stress_coexpression")
        # Rank score: harmonic-ish blend of available ranks (lower better)
        parts = [x for x in (s_rank, t_rank) if x is not None]
        if not parts:
            continue
        rank_score = len(parts) / sum(1.0 / x for x in parts)
        combined.append(
            {
                "protein_id": pid,
                "gene_name": genes.get(pid) or "",
                "tier": tiers.get(pid) or "",
                "score": round(scores.get(pid, 0.0), 6),
                "hub_sources": "+".join(sources),
                "string_rank": s_rank if s_rank is not None else "",
                "stress_rank": t_rank if t_rank is not None else "",
                "combined_rank_score": round(rank_score, 4),
                "string_n_neighbors": r.get("string_n_neighbors", ""),
                "string_n_tier_ab": r.get("string_n_tier_ab", ""),
                "string_enrichment_fold": (
                    round(r["string_enrichment_fold"], 6)
                    if r.get("string_enrichment_fold") is not None
                    else ""
                ),
                "string_hypergeom_pvalue": (
                    f"{r['string_hypergeom_pvalue']:.6g}"
                    if r.get("string_hypergeom_pvalue") is not None
                    else ""
                ),
                "stress_n_partners_r_ge_0.7": r.get("stress_n_partners_r_ge_0.7", ""),
                "stress_top_partner_id": r.get("stress_top_partner_id", ""),
                "stress_top_partner_r": r.get("stress_top_partner_r", ""),
            }
        )
    combined.sort(
        key=lambda r: (
            r["combined_rank_score"],
            r["string_rank"] if r["string_rank"] != "" else 10**9,
            r["stress_rank"] if r["stress_rank"] != "" else 10**9,
            -r["score"],
        )
    )
    for i, r in enumerate(combined, start=1):
        r["combined_rank"] = i

    top = combined[:top_n]
    cols = [
        "combined_rank",
        "protein_id",
        "gene_name",
        "tier",
        "score",
        "hub_sources",
        "combined_rank_score",
        "string_rank",
        "string_n_neighbors",
        "string_n_tier_ab",
        "string_enrichment_fold",
        "string_hypergeom_pvalue",
        "stress_rank",
        "stress_n_partners_r_ge_0.7",
        "stress_top_partner_id",
        "stress_top_partner_r",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in top:
            w.writerow(r)

    # Also write full ranked tables as gzip sidecars next to hubs TSV
    full_path = out_path.with_name("stress_seed_hubs_full.tsv.gz")
    with gzip.open(full_path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in combined:
            w.writerow(r)

    summary = {
        "layer": "string_stress_seed_hubs",
        "version": "1.5.0",
        "min_string_neighbors": min_string_neighbors,
        "n_string_tier_ab_ranked": len(string_rows),
        "n_stress_seeds_ranked": len(stress_rows),
        "n_combined": len(combined),
        "top_n_deposited": len(top),
        "paths": {
            "top_hubs": str(out_path.relative_to(ROOT)),
            "full_hubs": str(full_path.relative_to(ROOT)),
        },
        "top_hubs": top[:10],
        "honesty": (
            "STRING combined_score neighborhoods are descriptive adjacency "
            "(experiments/databases/coexpression/textmining mixed) — not LLPS evidence. "
            "Stress hubs use Pearson r on short Atlas log2FC vectors from the modules "
            "pass (Tier A seeds). Combined rank blends available STRING and stress ranks."
        ),
    }
    summary_path = out_path.with_name("stress_seed_hubs_summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        f"Hubs: STRING Tier A∪B ranked={len(string_rows):,} "
        f"(min n={min_string_neighbors}); stress seeds={len(stress_rows):,}; "
        f"top deposited={len(top)} → {out_path}",
        flush=True,
    )
    return summary


def extract_vit_ids(gene_name: str | None, annotations_json: str | None) -> set[str]:
    texts = [gene_name or ""]
    if annotations_json:
        try:
            a = json.loads(annotations_json)
            texts.append(a.get("gene_names") or "")
            texts.append(a.get("gene_name") or "")
        except json.JSONDecodeError:
            pass
    hits = set()
    for t in texts:
        for m in VIT_RE.findall(t):
            hits.add(m)
    return hits


def analysis_isoforms(
    con: sqlite3.Connection,
    *,
    out_dir: Path,
    min_delta: float,
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    by_vit: dict[str, dict[str, tuple]] = defaultdict(dict)
    n_with_vit = 0
    for r in con.execute(
        "SELECT protein_id, gene_name, tier, score, annotations, protein_name "
        "FROM proteins"
    ):
        vits = extract_vit_ids(r["gene_name"], r["annotations"])
        if not vits:
            continue
        n_with_vit += 1
        for vit in vits:
            by_vit[vit][r["protein_id"]] = (
                float(r["score"] or 0.0),
                r["tier"] or "",
                r["gene_name"] or "",
                r["protein_name"] or "",
            )

    multi = {vit: accs for vit, accs in by_vit.items() if len(accs) >= 2}
    delta_rows: list[dict] = []
    for vit, accs in multi.items():
        items = sorted(
            (
                {
                    "protein_id": pid,
                    "score": sc,
                    "tier": tier,
                    "gene_name": gn,
                    "protein_name": pn,
                }
                for pid, (sc, tier, gn, pn) in accs.items()
            ),
            key=lambda x: (-x["score"], x["protein_id"]),
        )
        scores = [x["score"] for x in items]
        delta = max(scores) - min(scores)
        tiers = sorted({x["tier"] for x in items if x["tier"]})
        delta_rows.append(
            {
                "vit_gene": vit,
                "n_accessions": len(items),
                "score_max": round(max(scores), 6),
                "score_min": round(min(scores), 6),
                "score_delta": round(delta, 6),
                "score_mean": round(sum(scores) / len(scores), 6),
                "tiers": ",".join(tiers),
                "max_protein_id": items[0]["protein_id"],
                "max_gene_name": items[0]["gene_name"],
                "min_protein_id": items[-1]["protein_id"],
                "min_gene_name": items[-1]["gene_name"],
                "accessions_json": json.dumps(items, separators=(",", ":")),
            }
        )
    delta_rows.sort(key=lambda r: (-r["score_delta"], -r["n_accessions"], r["vit_gene"]))

    high = [r for r in delta_rows if r["score_delta"] >= min_delta]
    # If the proteome is nearly 1:1, still deposit all multi-accession VIT genes.
    deposit = high if high else delta_rows

    tsv_path = out_dir / "multi_accession_vit_score_deltas.tsv"
    cols = [
        "vit_gene",
        "n_accessions",
        "score_max",
        "score_min",
        "score_delta",
        "score_mean",
        "tiers",
        "max_protein_id",
        "max_gene_name",
        "min_protein_id",
        "min_gene_name",
        "accessions_json",
    ]
    with tsv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in deposit:
            w.writerow(r)

    summary = {
        "layer": "isoform_multi_protein_genes",
        "version": "1.5.0",
        "n_proteins_with_vit": n_with_vit,
        "n_vit_loci": len(by_vit),
        "n_multi_accession_vit": len(multi),
        "min_delta": min_delta,
        "n_high_delta_deposited": len(deposit),
        "max_observed_delta": delta_rows[0]["score_delta"] if delta_rows else 0.0,
        "paths": {"deltas": str(tsv_path.relative_to(ROOT))},
        "top_deltas": deposit[:20],
        "honesty": (
            "Grouping key = VIT_* locus IDs parsed from UniProt gene_name / gene_names. "
            "UP000009183 is nearly one UniProt accession per VIT locus, so multi-accession "
            "isoform sets are rare. Score deltas are PlantLLPS hybrid score spreads across "
            "accessions of the same VIT — not experimentally confirmed splice-isoform LLPS "
            "differences. Shared gene symbols that map to distinct VIT IDs (paralogs) are "
            "intentionally excluded from this deposit."
        ),
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (out_dir / "README.md").write_text(
        f"""# Isoform / multi-protein gene score deltas (analysis 7)

PlantLLPS score spreads for UniProt accessions that share the same ``VIT_*`` locus.

## Headline

- Proteins with a parseable VIT ID: **{n_with_vit:,}**
- Distinct VIT loci: **{len(by_vit):,}**
- Multi-accession VIT loci: **{len(multi):,}**
- Deposited rows (delta ≥ {min_delta} or all multi if none): **{len(deposit):,}**
- Max observed delta: **{summary['max_observed_delta']}**

## Deposit

| File | Contents |
|------|----------|
| `multi_accession_vit_score_deltas.tsv` | High-delta (or all multi-accession) VIT genes |
| `summary.json` | Counts + honesty |

## Reproduce

```bash
python3 scripts/analyze_stress_network_isoforms.py
```

## Honesty

{summary["honesty"]}
""",
        encoding="utf-8",
    )
    print(
        f"Isoforms: VIT loci={len(by_vit):,} multi={len(multi):,} "
        f"deposited={len(deposit):,} max_delta={summary['max_observed_delta']}",
        flush=True,
    )
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--stress", type=Path, default=STRESS_SIDECAR)
    ap.add_argument("--deepdive-dir", type=Path, default=DEEPDIVE_DIR)
    ap.add_argument("--hubs-tsv", type=Path, default=HUBS_TSV)
    ap.add_argument("--isoform-dir", type=Path, default=ISOFORM_DIR)
    ap.add_argument("--top-hubs", type=int, default=100)
    ap.add_argument("--min-string-neighbors", type=int, default=5)
    ap.add_argument("--min-term-count", type=int, default=10)
    ap.add_argument("--min-overlap", type=int, default=3)
    ap.add_argument("--min-isoform-delta", type=float, default=0.05)
    args = ap.parse_args()

    for req in (
        args.db,
        args.stress,
        STRESS_ENRICH,
        STRESS_MEMBERS,
        STRESS_SEEDS,
        STRING_STATS,
    ):
        if not Path(req).exists():
            raise SystemExit(f"Required input missing: {req}")

    stress = load_json_gz(args.stress)
    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row

    print("=== Analysis 2: stress module deep dive ===", flush=True)
    deep = analysis_stress_deepdive(
        con,
        stress,
        out_dir=args.deepdive_dir,
        min_term_count=args.min_term_count,
        min_overlap=args.min_overlap,
    )

    print("=== Analysis 6: STRING / stress seed hubs ===", flush=True)
    hubs = analysis_hubs(
        con,
        out_path=args.hubs_tsv,
        top_n=args.top_hubs,
        min_string_neighbors=args.min_string_neighbors,
    )

    print("=== Analysis 7: isoform / multi-protein genes ===", flush=True)
    iso = analysis_isoforms(
        con,
        out_dir=args.isoform_dir,
        min_delta=args.min_isoform_delta,
    )

    combined = {
        "version": "1.5.0",
        "layer": "atlas_analyses_stress_network_isoforms",
        "analyses": {
            "2_stress_module_deepdive": {
                k: deep[k]
                for k in deep
                if k not in ("top_primary_terms",)
            },
            "6_string_stress_seed_hubs": {
                k: hubs[k] for k in hubs if k not in ("top_hubs",)
            },
            "7_isoform_multi_protein_genes": {
                k: iso[k] for k in iso if k not in ("top_deltas",)
            },
        },
        "honesty": {
            "deepdive": deep["honesty"],
            "hubs": hubs["honesty"],
            "isoforms": iso["honesty"],
        },
    }
    meta_path = ROOT / "data" / "atlas_analyses_summary.json"
    meta_path.write_text(json.dumps(combined, indent=2) + "\n", encoding="utf-8")

    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("atlas_analyses_stress_network_isoforms", json.dumps(combined)),
    )
    con.commit()
    con.close()

    print(json.dumps(combined, indent=2)[:2500], flush=True)
    print(f"Wrote {meta_path}", flush=True)


if __name__ == "__main__":
    main()
