#!/usr/bin/env python3
"""Structure–score joint: high PlantLLPS ∩ low-pLDDT segments.

Overlaps PlantLLPS Tier A∪B with AlphaFold low-pLDDT (pLDDT<50) disorder
fraction / residue counts deposited in the atlas. Writes summary stats,
optional hypergeometric enrichment of high low-pLDDT proteins among Tier A∪B,
and a per-protein flag table under ``data/structure_score/``.

Honesty
-------
- Uses deposited ``af_frac_very_low`` / ``af_n_disorder_residues`` / IDR segments
  with ``src=af_plddt_lt50`` (full-residue AF JSON is gitignored for size).
- High low-pLDDT among high PlantLLPS is expected for IDR-rich condensates but
  is **not** proof of LLPS; treat as descriptive structure–score co-occurrence.
- FDR applies only to the discrete threshold tests listed in the enrichment
  table; the per-protein flag list is exploratory.

Example::

    python3 scripts/enrich_structure_score_joint.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
OUT_DIR = ROOT / "data" / "structure_score"

# Thresholds tested for enrichment (BH across this list)
AF_FRAC_THRESHOLDS = (0.15, 0.25, 0.35, 0.50)


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
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


def mean(xs: list[float]) -> float | None:
    return sum(xs) / len(xs) if xs else None


def median(xs: list[float]) -> float | None:
    if not xs:
        return None
    s = sorted(xs)
    m = len(s) // 2
    if len(s) % 2:
        return s[m]
    return 0.5 * (s[m - 1] + s[m])


def af_idr_residues(idr_json: str | None) -> int:
    """Count residues covered by af_plddt_lt50 IDR segments (union approx via sum)."""
    if not idr_json:
        return 0
    try:
        segs = json.loads(idr_json)
    except json.JSONDecodeError:
        return 0
    if not isinstance(segs, list):
        return 0
    n = 0
    for s in segs:
        if not isinstance(s, dict):
            continue
        src = str(s.get("src") or s.get("source") or "")
        if "af_plddt" not in src:
            continue
        start = s.get("s", s.get("start"))
        end = s.get("e", s.get("end"))
        if start is None or end is None:
            continue
        n += max(0, int(end) - int(start) + 1)
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument(
        "--flag-threshold",
        type=float,
        default=0.25,
        help="af_frac_very_low threshold for per-protein high_low_plddt flag",
    )
    args = ap.parse_args()

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(str(args.db))
    rows = list(
        con.execute(
            """
            SELECT protein_id, gene_name, protein_name, tier, score, length,
                   af_mean_plddt, af_frac_very_low, af_n_disorder_residues,
                   idr_segments_json, idr_n_residues
            FROM proteins
            """
        )
    )

    N = len(rows)
    fg = {r[0] for r in rows if r[3] in ("A", "B")}
    K = len(fg)

    af_covered = [r for r in rows if r[7] is not None]
    N_af = len(af_covered)
    fg_af = [r for r in af_covered if r[0] in fg]
    bg_only_af = [r for r in af_covered if r[0] not in fg]

    enrich_rows = []
    for thr in AF_FRAC_THRESHOLDS:
        high = {r[0] for r in af_covered if float(r[7]) >= thr}
        n = len(high)
        k = len(high & fg)
        # Background for this test: AF-covered proteome
        fold = (k / K) / (n / N_af) if K and n and N_af else 0.0
        # Hypergeometric: drawing n high-low-pLDDT proteins from N_af, K Tier A∪B in AF set
        K_af = len(fg_af)
        pval = hypergeom_sf(k, N_af, K_af, n)
        enrich_rows.append(
            {
                "test": f"af_frac_very_low>={thr}",
                "threshold": thr,
                "n_background_af": N_af,
                "n_tier_ab_af": K_af,
                "n_high_low_plddt": n,
                "n_high_intersect_tier_ab": k,
                "enrichment_fold": round(fold, 6),
                "hypergeom_pvalue": pval,
            }
        )
    fdrs = bh_fdr([r["hypergeom_pvalue"] for r in enrich_rows])
    for r, fdr in zip(enrich_rows, fdrs):
        r["bh_fdr"] = round(fdr, 6)
        r["hypergeom_pvalue"] = round(r["hypergeom_pvalue"], 6)

    flag_thr = float(args.flag_threshold)
    flag_rows = []
    for (
        pid,
        gene,
        pname,
        tier,
        score,
        length,
        mean_p,
        frac,
        n_dis,
        idr_json,
        idr_n,
    ) in rows:
        if tier not in ("A", "B"):
            continue
        if frac is None:
            continue
        af_idr_n = af_idr_residues(idr_json)
        high_flag = float(frac) >= flag_thr
        has_seg = (n_dis or 0) > 0 or af_idr_n > 0
        flag_rows.append(
            {
                "protein_id": pid,
                "gene_name": gene or "",
                "protein_name": pname or "",
                "tier": tier,
                "score_plantllps": round(float(score), 8),
                "length": length or "",
                "af_mean_plddt": round(float(mean_p), 4) if mean_p is not None else "",
                "af_frac_very_low": round(float(frac), 6),
                "af_n_disorder_residues": int(n_dis or 0),
                "af_idr_segment_residues": af_idr_n,
                "idr_n_residues": int(idr_n or 0),
                "flag_high_low_plddt": int(high_flag),
                "flag_has_low_plddt_segment": int(has_seg),
                "joint_high_plantllps_low_plddt": int(high_flag),
            }
        )

    flag_rows.sort(
        key=lambda r: (
            -r["joint_high_plantllps_low_plddt"],
            -r["af_frac_very_low"],
            -r["score_plantllps"],
        )
    )
    for i, r in enumerate(flag_rows, start=1):
        r["rank"] = i

    enrich_tsv = out_dir / "enrichment_thresholds.tsv"
    with enrich_tsv.open("w", encoding="utf-8", newline="") as fh:
        fields = [
            "test",
            "threshold",
            "n_background_af",
            "n_tier_ab_af",
            "n_high_low_plddt",
            "n_high_intersect_tier_ab",
            "enrichment_fold",
            "hypergeom_pvalue",
            "bh_fdr",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in enrich_rows:
            w.writerow(r)

    flags_tsv = out_dir / "tier_ab_low_plddt_flags.tsv"
    with flags_tsv.open("w", encoding="utf-8", newline="") as fh:
        fields = [
            "rank",
            "protein_id",
            "gene_name",
            "protein_name",
            "tier",
            "score_plantllps",
            "length",
            "af_mean_plddt",
            "af_frac_very_low",
            "af_n_disorder_residues",
            "af_idr_segment_residues",
            "idr_n_residues",
            "flag_high_low_plddt",
            "flag_has_low_plddt_segment",
            "joint_high_plantllps_low_plddt",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in flag_rows:
            w.writerow(r)

    def tier_stats(subset: list) -> dict:
        fracs = [float(r[7]) for r in subset if r[7] is not None]
        means = [float(r[6]) for r in subset if r[6] is not None]
        return {
            "n": len(subset),
            "n_with_af": len(fracs),
            "mean_af_frac_very_low": round(mean(fracs), 6) if fracs else None,
            "median_af_frac_very_low": round(median(fracs), 6) if fracs else None,
            "mean_af_mean_plddt": round(mean(means), 4) if means else None,
            "frac_af_frac_ge_flag_thr": round(
                sum(1 for x in fracs if x >= flag_thr) / len(fracs), 6
            )
            if fracs
            else None,
        }

    by_tier = {}
    for t in ("A", "B", "C", "D"):
        by_tier[t] = tier_stats([r for r in rows if r[3] == t])

    n_joint = sum(1 for r in flag_rows if r["joint_high_plantllps_low_plddt"])
    honesty = (
        "Joint overlap of PlantLLPS Tier A∪B with AlphaFold low-pLDDT "
        f"(af_frac_very_low ≥ {flag_thr} for per-protein flags). Hypergeometric tests "
        "across discrete fraction thresholds are BH-controlled; the ranked flag table "
        "is exploratory co-occurrence, not proof of LLPS. AF coverage incomplete "
        f"({N_af}/{N} proteins)."
    )

    summary = {
        "version": "1.5.0",
        "layer": "structure_score_joint",
        "n_proteins": N,
        "n_af_covered": N_af,
        "n_tier_ab": K,
        "n_tier_ab_af": len(fg_af),
        "flag_threshold": flag_thr,
        "n_flag_rows": len(flag_rows),
        "n_joint_high_plantllps_low_plddt": n_joint,
        "by_tier": by_tier,
        "enrichment": enrich_rows,
        "n_enrichment_fdr_lt_0.05": sum(1 for r in enrich_rows if r["bh_fdr"] < 0.05),
        "paths": {
            "enrichment": str(enrich_tsv.relative_to(ROOT)),
            "flags": str(flags_tsv.relative_to(ROOT)),
            "ui_tables": "data/structure_score/ui_tables.json",
        },
        "honesty": honesty,
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    ui = {
        "parameters": {
            "n_proteins": N,
            "n_af_covered": N_af,
            "n_tier_ab": K,
            "n_tier_ab_af": len(fg_af),
            "flag_threshold": flag_thr,
            "n_joint": n_joint,
            "n_enrichment_fdr_lt_0.05": summary["n_enrichment_fdr_lt_0.05"],
        },
        "honesty": honesty,
        "by_tier": by_tier,
        "enrichment": enrich_rows,
        "flags": flag_rows[:200],
    }
    (out_dir / "ui_tables.json").write_text(
        json.dumps(ui, indent=2) + "\n", encoding="utf-8"
    )

    (out_dir / "README.md").write_text(
        f"""# Structure–score joint (PlantLLPS × low-pLDDT)

Overlap of high PlantLLPS (Tier A∪B) with AlphaFold low-pLDDT segments /
``af_frac_very_low``.

## Reproduce

```bash
python3 scripts/enrich_structure_score_joint.py
```

## Deposits

| File | Contents |
|------|----------|
| `enrichment_thresholds.tsv` | Hypergeometric tests at af_frac thresholds + BH FDR |
| `tier_ab_low_plddt_flags.tsv` | Per-protein flags for Tier A∪B with AF coverage |
| `ui_tables.json` | Compact payload for Node `/analyses` |
| `summary.json` | Run summary + per-tier stats |

## Headline

- AF-covered **{N_af:,}** / {N:,}; Tier A∪B with AF **{len(fg_af):,}**
- Joint flags (af_frac ≥ {flag_thr}): **{n_joint:,}**
- Enrichment tests FDR&lt;0.05: **{summary["n_enrichment_fdr_lt_0.05"]}** / {len(enrich_rows)}

## Honesty

{honesty}
""",
        encoding="utf-8",
    )

    meta = {k: summary[k] for k in summary if k != "enrichment"}
    meta["enrichment"] = enrich_rows
    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("structure_score_joint", json.dumps(meta)),
    )
    con.commit()
    con.close()

    print(json.dumps({k: summary[k] for k in summary if k not in ("by_tier", "enrichment")}, indent=2))
    print("Enrichment:", flush=True)
    for r in enrich_rows:
        print(
            f"  {r['test']}: ∩={r['n_high_intersect_tier_ab']}/{r['n_high_low_plddt']} "
            f"fold={r['enrichment_fold']} FDR={r['bh_fdr']}",
            flush=True,
        )
    print(f"Wrote {out_dir} ({len(flag_rows)} flag rows)", flush=True)


if __name__ == "__main__":
    main()
