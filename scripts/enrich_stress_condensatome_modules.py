#!/usr/bin/env python3
"""Condensatome × stress DE enrichment + coexpression modules (NAR+ omics layer).

Builds on the stress RNA-seq joins (E-MTAB-4324, E-GEOD-67932) already deposited
by ``enrich_stress_rnaseq.py``:

1. **DE × Tier A∪B enrichment** — for each Atlas contrast and direction (up / down /
   any), test whether significant DE genes are enriched for PlantLLPS Tier A∪B
   proteins (hypergeometric vs the stress-matched atlas background).
2. **Coexpression / response modules for high-tier proteins** —
   (a) **Pattern modules**: Tier A∪B proteins sharing the same signed significant
   DE signature across contrasts (−1 / 0 / +1); (b) **Seed neighborhoods**:
   Pearson r of stress log2FC profiles for Tier A seeds (top-N partners, r ≥ 0.7).
3. **Gene-set deposit** — GMT + membership TSV for DE∩Tier sets and modules.

Writes under ``data/stress_modules/`` (sidecars; keeps ``grape_condensatome.sqlite``
under GitHub's 100 MiB limit). Records a compact summary in SQLite ``meta``.

Honesty
-------
- Call rule matches stress join: p < 0.05 and |log2FC| ≥ 1.
- Enrichment is descriptive overlap, not causal condensate membership.
- Pattern modules group shared stress DE call signatures (few contrasts), not a
  dedicated multi-condition coexpression atlas. Berry TPM (E-MTAB-8758) has only
  6 assay groups and is **not** used as the primary module signal.
- Seed neighborhoods use Pearson r on the same short log2FC vectors.
- Modules are hypothesis-generating; not WGCNA / soft-threshold networks.

Example::

    python3 scripts/enrich_stress_condensatome_modules.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
STRESS_SIDECAR = ROOT / "data" / "expression_atlas" / "stress_by_protein.json.gz"
OUT_DIR = ROOT / "data" / "stress_modules"

P_THR = 0.05
LFC_THR = 1.0
DEFAULT_R_THR = 0.85  # seed-neighborhood reporting threshold (edges at r≥0.7 listed)
DEFAULT_MIN_MODULE = 5
DEFAULT_TOP_N = 15
DEFAULT_MIN_CONTRASTS = 4


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
    """Benjamini–Hochberg FDR (independent / positive regression)."""
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


def pearson(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 3:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    num = 0.0
    dx = 0.0
    dy = 0.0
    for a, b in zip(xs, ys):
        xa = a - mx
        yb = b - my
        num += xa * yb
        dx += xa * xa
        dy += yb * yb
    if dx <= 0 or dy <= 0:
        return None
    return num / math.sqrt(dx * dy)


def load_stress(path: Path) -> dict[str, dict]:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def contrast_keys(stress: dict[str, dict]) -> list[tuple[str, str, str]]:
    """Return sorted (accession, contrast_label, contrast_id) keys present in data."""
    keys: dict[tuple[str, str], str] = {}
    for blob in stress.values():
        for acc, info in blob.items():
            for lab, comp in (info.get("comparisons") or {}).items():
                cid = str(comp.get("contrast_id") or lab)
                keys.setdefault((acc, lab), cid)
    out = [(acc, lab, cid) for (acc, lab), cid in keys.items()]
    out.sort(key=lambda t: (t[0], t[2], t[1]))
    return out


def pattern_label(sig: tuple[int, ...], axis: list[tuple[str, str]]) -> str:
    parts = []
    for bit, (acc, lab) in zip(sig, axis):
        if bit == 0:
            continue
        short = acc.replace("E-", "")
        parts.append(f"{short}:{'up' if bit > 0 else 'down'}")
    return "; ".join(parts) if parts else "no_DE"


def slug(s: str) -> str:
    out = []
    for ch in s.lower():
        if ch.isalnum():
            out.append(ch)
        elif ch in (" ", "-", "/", ";", ":", ",", "."):
            out.append("_")
    text = "".join(out)
    while "__" in text:
        text = text.replace("__", "_")
    return text.strip("_")[:80] or "set"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--stress", type=Path, default=STRESS_SIDECAR)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--r-threshold", type=float, default=DEFAULT_R_THR)
    ap.add_argument("--min-module-size", type=int, default=DEFAULT_MIN_MODULE)
    ap.add_argument("--top-n", type=int, default=DEFAULT_TOP_N)
    ap.add_argument("--min-contrasts", type=int, default=DEFAULT_MIN_CONTRASTS)
    ap.add_argument(
        "--seeds",
        choices=("A", "AB"),
        default="A",
        help="High-tier seed set for coexpression neighborhoods (default: Tier A)",
    )
    args = ap.parse_args()

    if not args.db.exists():
        raise SystemExit(f"DB not found: {args.db}")
    if not args.stress.exists():
        raise SystemExit(
            f"Stress sidecar missing: {args.stress}\n"
            "Run: python3 scripts/enrich_stress_rnaseq.py"
        )

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    stress = load_stress(args.stress)
    print(f"Loaded stress sidecar: {len(stress):,} proteins", flush=True)

    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row
    rows = list(
        con.execute(
            "SELECT protein_id, gene_name, tier, score FROM proteins"
        )
    )
    tiers = {r["protein_id"]: (r["tier"] or "") for r in rows}
    scores = {r["protein_id"]: float(r["score"] or 0.0) for r in rows}
    genes = {r["protein_id"]: (r["gene_name"] or "") for r in rows}
    atlas_ids = set(tiers)
    tier_ab = {pid for pid, t in tiers.items() if t in ("A", "B")}
    seed_tiers = {"A"} if args.seeds == "A" else {"A", "B"}
    seeds = {pid for pid, t in tiers.items() if t in seed_tiers}
    # Background = proteins with any stress DE join (same universe as DE calls)
    bg = set(stress.keys()) & atlas_ids
    bg_ab = bg & tier_ab
    N = len(bg)
    K = len(bg_ab)
    bg_rate = K / N if N else 0.0
    print(
        f"Background stress-matched={N:,}; Tier A∪B in bg={K:,} ({bg_rate:.4f}); "
        f"seeds={args.seeds} n={len(seeds):,}",
        flush=True,
    )

    contrasts = contrast_keys(stress)
    print(f"Contrasts: {len(contrasts)}", flush=True)

    # ---- DE enrichment -------------------------------------------------
    enrich_rows: list[dict] = []
    gene_sets: dict[str, dict] = {}  # set_id → {name, desc, members:set}

    for acc, lab, cid in contrasts:
        buckets = {
            "up": set(),
            "down": set(),
            "any": set(),
        }
        for pid in bg:
            blob = stress.get(pid) or {}
            info = blob.get(acc) or {}
            comp = (info.get("comparisons") or {}).get(lab)
            if not comp or not comp.get("significant"):
                continue
            lfc = float(comp.get("log2fc") or 0.0)
            buckets["any"].add(pid)
            if lfc >= LFC_THR:
                buckets["up"].add(pid)
            elif lfc <= -LFC_THR:
                buckets["down"].add(pid)

        for direction, de_set in buckets.items():
            n = len(de_set)
            overlap = de_set & bg_ab
            k = len(overlap)
            fold = (k / n / bg_rate) if n and bg_rate else None
            pval = hypergeom_sf(k, N, K, n) if n else 1.0
            set_id = f"DE_{acc}_{cid}_{direction}_TierAB"
            set_name = f"{acc} {cid} {direction} ∩ Tier A∪B"
            gene_sets[set_id] = {
                "name": set_name,
                "description": (
                    f"Significant DE ({direction}) in {acc} contrast '{lab}' "
                    f"intersect PlantLLPS Tier A∪B; call p<{P_THR} |log2FC|>={LFC_THR}"
                ),
                "members": sorted(overlap),
                "kind": "de_tier_overlap",
                "accession": acc,
                "contrast_id": cid,
                "contrast_label": lab,
                "direction": direction,
            }
            # Also deposit full DE set (not just TierAB) for reuse
            full_id = f"DE_{acc}_{cid}_{direction}"
            gene_sets[full_id] = {
                "name": f"{acc} {cid} {direction}",
                "description": (
                    f"Significant DE ({direction}) in {acc} contrast '{lab}'; "
                    f"call p<{P_THR} |log2FC|>={LFC_THR}"
                ),
                "members": sorted(de_set),
                "kind": "de_full",
                "accession": acc,
                "contrast_id": cid,
                "contrast_label": lab,
                "direction": direction,
            }
            enrich_rows.append(
                {
                    "accession": acc,
                    "contrast_id": cid,
                    "contrast_label": lab,
                    "direction": direction,
                    "n_background": N,
                    "n_tier_ab_background": K,
                    "background_tier_ab_rate": round(bg_rate, 6),
                    "n_de": n,
                    "n_de_tier_ab": k,
                    "frac_de_tier_ab": round(k / n, 6) if n else None,
                    "enrichment_fold": round(fold, 6) if fold is not None else None,
                    "hypergeom_pvalue": pval,
                    "gene_set_id": set_id,
                    "full_de_set_id": full_id,
                }
            )

    pvals = [r["hypergeom_pvalue"] for r in enrich_rows]
    fdrs = bh_fdr(pvals)
    for r, fdr in zip(enrich_rows, fdrs):
        r["bh_fdr"] = round(fdr, 6)
        r["hypergeom_pvalue"] = round(r["hypergeom_pvalue"], 6)

    enrich_rows.sort(
        key=lambda r: (r["bh_fdr"], r["hypergeom_pvalue"], -(r["enrichment_fold"] or 0))
    )

    # ---- Pattern modules + seed coexpression neighborhoods -------------
    axis = [(acc, lab) for acc, lab, _cid in contrasts]
    axis_cids = [cid for _acc, _lab, cid in contrasts]

    # Signed significant-DE patterns among Tier A∪B (−1 / 0 / +1 per contrast)
    pattern_members: dict[tuple[int, ...], list[str]] = defaultdict(list)
    profiles: dict[str, list[float | None]] = {}
    for pid in bg & tier_ab:
        blob = stress.get(pid) or {}
        sig: list[int] = []
        vec: list[float | None] = []
        for acc, lab in axis:
            info = blob.get(acc) or {}
            comp = (info.get("comparisons") or {}).get(lab)
            if not comp or comp.get("log2fc") is None:
                sig.append(0)
                vec.append(None)
                continue
            lfc = float(comp["log2fc"])
            vec.append(lfc)
            if comp.get("significant"):
                if lfc >= LFC_THR:
                    sig.append(1)
                elif lfc <= -LFC_THR:
                    sig.append(-1)
                else:
                    sig.append(0)
            else:
                sig.append(0)
        key = tuple(sig)
        if any(bit != 0 for bit in key):
            pattern_members[key].append(pid)
        n_obs = sum(1 for v in vec if v is not None)
        if n_obs >= args.min_contrasts:
            profiles[pid] = vec

    print(
        f"Non-empty DE patterns among Tier A∪B: {len(pattern_members):,}; "
        f"profiled for correlation: {len(profiles):,}",
        flush=True,
    )

    module_rows: list[dict] = []
    pattern_items = sorted(
        pattern_members.items(),
        key=lambda kv: (
            -len(kv[1]),
            -sum(scores.get(p, 0.0) for p in kv[1]),
            kv[0],
        ),
    )
    mi = 0
    for sig, members in pattern_items:
        if len(members) < args.min_module_size:
            continue
        mi += 1
        mid = f"CX_stress_mod_{mi:03d}"
        members_sorted = sorted(members, key=lambda p: (-scores.get(p, 0.0), p))
        n_a = sum(1 for p in members if tiers.get(p) == "A")
        n_b = sum(1 for p in members if tiers.get(p) == "B")
        hub = members_sorted[0]
        plab = pattern_label(sig, axis)
        gene_sets[mid] = {
            "name": f"Stress response pattern module {mi}",
            "description": (
                f"Tier A∪B proteins sharing signed significant DE pattern "
                f"[{plab}] across {', '.join(axis_cids)}; hub={hub}"
            ),
            "members": members_sorted,
            "kind": "response_pattern_module",
            "hub": hub,
            "pattern": list(sig),
            "pattern_label": plab,
            "n_tier_a": n_a,
            "n_tier_b": n_b,
        }
        module_rows.append(
            {
                "module_id": mid,
                "n_members": len(members_sorted),
                "n_tier_a": n_a,
                "n_tier_b": n_b,
                "pattern": "".join("+" if b > 0 else ("-" if b < 0 else "0") for b in sig),
                "pattern_label": plab,
                "hub_protein_id": hub,
                "hub_gene_name": genes.get(hub) or "",
                "hub_tier": tiers.get(hub) or "",
                "hub_score": round(scores.get(hub, 0.0), 6),
                "mean_score": round(
                    sum(scores.get(p, 0.0) for p in members_sorted) / len(members_sorted),
                    6,
                ),
                "top_members": ",".join(members_sorted[:10]),
            }
        )
    print(
        f"Pattern modules (size≥{args.min_module_size}): {len(module_rows)}",
        flush=True,
    )

    # Seed neighborhoods via Pearson r on log2FC profiles
    r_thr = args.r_threshold
    seed_nb: dict[str, list[tuple[str, float]]] = {
        s: [] for s in seeds if s in profiles
    }
    seed_pids = sorted(seed_nb.keys(), key=lambda p: (-scores.get(p, 0.0), p))
    partner_pool = sorted(profiles.keys(), key=lambda p: (-scores.get(p, 0.0), p))
    print(
        f"Computing seed neighborhoods for {len(seed_pids):,} Tier {args.seeds} "
        f"seeds vs {len(partner_pool):,} profiled Tier A∪B …",
        flush=True,
    )
    for i, a in enumerate(seed_pids):
        if i % 200 == 0:
            print(f"  seed neighborhoods {i}/{len(seed_pids)} …", flush=True)
        va = profiles[a]
        for b in partner_pool:
            if a == b:
                continue
            vb = profiles[b]
            xs: list[float] = []
            ys: list[float] = []
            for xa, yb in zip(va, vb):
                if xa is None or yb is None:
                    continue
                xs.append(xa)
                ys.append(yb)
            if len(xs) < args.min_contrasts:
                continue
            r = pearson(xs, ys)
            if r is None or r < 0.7:
                continue
            seed_nb[a].append((b, r))

    seed_rows: list[dict] = []
    edge_ge_thr = 0
    for sid in seed_pids:
        partners = sorted(seed_nb[sid], key=lambda t: (-t[1], -scores.get(t[0], 0.0)))
        n_ge_thr = sum(1 for _p, r in partners if r >= r_thr)
        edge_ge_thr += n_ge_thr
        top = partners[: args.top_n]
        if not top:
            continue
        set_id = f"CX_seed_{sid}"
        members = [sid] + [p for p, _r in top]
        gene_sets[set_id] = {
            "name": f"Coexpression neighborhood of {sid}",
            "description": (
                f"Tier {args.seeds} seed {sid} + top-{args.top_n} partners with "
                f"Pearson r≥0.7 on stress log2FC profiles"
            ),
            "members": members,
            "kind": "seed_neighborhood",
            "seed": sid,
        }
        seed_rows.append(
            {
                "protein_id": sid,
                "gene_name": genes.get(sid) or "",
                "tier": tiers.get(sid) or "",
                "score": round(scores.get(sid, 0.0), 6),
                "n_partners_r_ge_0.7": len(partners),
                "n_partners_r_ge_thr": n_ge_thr,
                "top_partner_id": top[0][0],
                "top_partner_r": round(top[0][1], 4),
                "top_partners_json": json.dumps(
                    [
                        {
                            "protein_id": p,
                            "r": round(r, 4),
                            "tier": tiers.get(p),
                            "score": round(scores.get(p, 0.0), 6),
                            "gene_name": genes.get(p) or "",
                        }
                        for p, r in top
                    ],
                    separators=(",", ":"),
                ),
                "gene_set_id": set_id,
            }
        )
    print(
        f"Seed neighborhoods deposited: {len(seed_rows):,} "
        f"(partner links r≥{r_thr}: {edge_ge_thr:,})",
        flush=True,
    )
    edge_count = edge_ge_thr

    # ---- Write deposits ------------------------------------------------
    enrich_tsv = out_dir / "de_enrichment.tsv"
    with enrich_tsv.open("w", newline="") as fh:
        cols = [
            "accession",
            "contrast_id",
            "contrast_label",
            "direction",
            "n_background",
            "n_tier_ab_background",
            "background_tier_ab_rate",
            "n_de",
            "n_de_tier_ab",
            "frac_de_tier_ab",
            "enrichment_fold",
            "hypergeom_pvalue",
            "bh_fdr",
            "gene_set_id",
            "full_de_set_id",
        ]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in enrich_rows:
            w.writerow(r)

    modules_tsv = out_dir / "coexpression_modules.tsv"
    with modules_tsv.open("w", newline="") as fh:
        cols = [
            "module_id",
            "n_members",
            "n_tier_a",
            "n_tier_b",
            "pattern",
            "pattern_label",
            "hub_protein_id",
            "hub_gene_name",
            "hub_tier",
            "hub_score",
            "mean_score",
            "top_members",
        ]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in module_rows:
            w.writerow(r)

    seeds_tsv = out_dir / "seed_neighborhoods.tsv.gz"
    with gzip.open(seeds_tsv, "wt", encoding="utf-8", newline="") as fh:
        cols = [
            "protein_id",
            "gene_name",
            "tier",
            "score",
            "n_partners_r_ge_0.7",
            "n_partners_r_ge_thr",
            "top_partner_id",
            "top_partner_r",
            "top_partners_json",
            "gene_set_id",
        ]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in seed_rows:
            w.writerow(r)

    members_path = out_dir / "module_members.tsv.gz"
    with gzip.open(members_path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["set_id", "kind", "protein_id", "gene_name", "tier", "score", "rank_in_set"]
        )
        for set_id, info in sorted(gene_sets.items()):
            for rank, pid in enumerate(info["members"], start=1):
                w.writerow(
                    [
                        set_id,
                        info.get("kind") or "",
                        pid,
                        genes.get(pid) or "",
                        tiers.get(pid) or "",
                        round(scores.get(pid, 0.0), 6),
                        rank,
                    ]
                )

    gmt_path = out_dir / "gene_sets.gmt.gz"
    with gzip.open(gmt_path, "wt", encoding="utf-8") as fh:
        for set_id, info in sorted(gene_sets.items()):
            # GMT: name, description, members… (use gene_name when present else protein_id)
            mem_labels = []
            for pid in info["members"]:
                g = (genes.get(pid) or "").strip()
                mem_labels.append(g if g else pid)
            desc = (info.get("description") or "").replace("\t", " ")
            fh.write("\t".join([set_id, desc, *mem_labels]) + "\n")

    # Compact UI JSON (enrichment + modules + seed table head)
    ui_payload = {
        "enrichment": enrich_rows,
        "modules": module_rows,
        "seed_neighborhoods_n": len(seed_rows),
        "seed_neighborhoods_top": seed_rows[:50],
        "parameters": {
            "p_threshold": P_THR,
            "lfc_threshold": LFC_THR,
            "r_threshold": r_thr,
            "min_module_size": args.min_module_size,
            "min_contrasts": args.min_contrasts,
            "top_n": args.top_n,
            "seeds": args.seeds,
            "n_background": N,
            "n_tier_ab_background": K,
            "background_tier_ab_rate": round(bg_rate, 6),
            "n_profiled_tier_ab": len(profiles),
            "n_seed_partner_links_r_ge_thr": edge_count,
            "n_modules": len(module_rows),
            "n_gene_sets": len(gene_sets),
            "module_definition": (
                "Tier A∪B proteins sharing the same signed significant DE pattern "
                f"across contrasts (min size {args.min_module_size})"
            ),
            "contrasts": [
                {"accession": a, "contrast_id": c, "contrast_label": l}
                for a, l, c in contrasts
            ],
            "accessions": sorted({a for a, _l, _c in contrasts}),
        },
    }
    ui_json = out_dir / "ui_tables.json"
    ui_json.write_text(json.dumps(ui_payload, indent=2) + "\n", encoding="utf-8")

    summary = {
        "version": "1.2.0",
        "layer": "stress_condensatome_modules",
        "accessions": sorted({a for a, _l, _c in contrasts}),
        "n_enrichment_tests": len(enrich_rows),
        "n_enrichment_fdr_lt_0.05": sum(1 for r in enrich_rows if r["bh_fdr"] < 0.05),
        "n_modules": len(module_rows),
        "n_seed_neighborhoods": len(seed_rows),
        "n_gene_sets": len(gene_sets),
        "n_profiled_tier_ab": len(profiles),
        "n_seed_partner_links_r_ge_thr": edge_count,
        "parameters": ui_payload["parameters"],
        "paths": {
            "de_enrichment": str(enrich_tsv.relative_to(ROOT)),
            "coexpression_modules": str(modules_tsv.relative_to(ROOT)),
            "seed_neighborhoods": str(seeds_tsv.relative_to(ROOT)),
            "module_members": str(members_path.relative_to(ROOT)),
            "gene_sets_gmt": str(gmt_path.relative_to(ROOT)),
            "ui_tables": str(ui_json.relative_to(ROOT)),
        },
        "top_enrichment": enrich_rows[:8],
        "honesty": (
            "DE call rule matches stress join (p<0.05, |log2FC|≥1). Hypergeometric "
            "enrichment of Tier A∪B among DE genes is descriptive. Response-pattern "
            "modules group Tier A∪B proteins by shared signed significant DE calls "
            "across Atlas contrasts; seed neighborhoods use Pearson r on the same short "
            "log2FC vectors — not a dedicated multi-condition coexpression atlas. Berry "
            "TPM E-MTAB-8758 has only 6 assay groups and is not the primary module "
            "signal. Not WGCNA."
        ),
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    readme = out_dir / "README.md"
    readme.write_text(
        f"""# Stress condensatome modules (NAR+ omics layer)

Condensatome × stress DE enrichment for **E-MTAB-4324** and **E-GEOD-67932**, plus
coexpression modules among high-tier (PlantLLPS Tier A∪B) proteins.

## Reproduce

```bash
python3 scripts/enrich_stress_rnaseq.py          # prerequisite sidecar
python3 scripts/enrich_stress_condensatome_modules.py
```

## Deposits

| File | Contents |
|------|----------|
| `de_enrichment.tsv` | Per-contrast × direction hypergeometric enrichment vs Tier A∪B |
| `coexpression_modules.tsv` | Response-pattern modules (shared signed DE calls) |
| `seed_neighborhoods.tsv.gz` | Tier {args.seeds} seed top-{args.top_n} partners (r≥0.7) |
| `module_members.tsv.gz` | Membership for every gene set |
| `gene_sets.gmt.gz` | GMT deposit (gene names; UniProt fallback) |
| `ui_tables.json` | Compact payload for Node `/modules` |
| `summary.json` | Run summary + honesty note |

## Parameters

- Background: stress-matched atlas proteins (**{N:,}**); Tier A∪B in background **{K:,}** (rate {bg_rate:.4f})
- DE call: p < {P_THR} and |log2FC| ≥ {LFC_THR}
- Pattern modules: shared −1/0/+1 signature across contrasts; size ≥ {args.min_module_size}
- Seed neighborhoods: Pearson r on stress log2FC (≥{args.min_contrasts} shared); top-{args.top_n} at r≥0.7
- Profiled Tier A∪B: **{len(profiles):,}**; seed partner links r≥{r_thr}: **{edge_count:,}**; modules: **{len(module_rows)}**; gene sets: **{len(gene_sets)}**

## Honesty

{summary["honesty"]}
""",
        encoding="utf-8",
    )

    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("stress_condensatome_modules", json.dumps(summary)),
    )
    con.commit()
    con.close()

    print(json.dumps({k: summary[k] for k in summary if k != "top_enrichment"}, indent=2))
    print("Top enrichment rows:", flush=True)
    for r in enrich_rows[:6]:
        print(
            f"  {r['accession']} {r['contrast_id']} {r['direction']}: "
            f"n_de={r['n_de']} ∩AB={r['n_de_tier_ab']} fold={r['enrichment_fold']} "
            f"p={r['hypergeom_pvalue']} FDR={r['bh_fdr']}",
            flush=True,
        )
    print(f"Wrote {out_dir}", flush=True)


if __name__ == "__main__":
    main()
