#!/usr/bin/env python3
"""Tier A∪B Pfam/InterPro family enrichment (hypergeometric + BH FDR).

Foreground = PlantLLPS Tier A∪B proteins. Background = full scored proteome.
Each Pfam / InterPro accession present in ≥ ``--min-foreground`` foreground proteins
is tested for over-representation (hypergeometric survival). FDR is controlled
within each database (Pfam and InterPro separately).

Deposits under ``data/family_enrichment/`` (sidecars; SQLite meta pointer only).

Honesty
-------
- Enrichment is descriptive domain bias among high-score proteins, not causal
  condensate membership.
- Families with FDR ≥ 0.05 are exploratory; cite FDR-controlled hits as primary.
- Domain IDs come from deposited UniProt/InterPro xrefs (``domain_pfam_ids`` /
  ``domain_interpro_ids``); coordinate labels are best-effort from ``domains_json``.

Example::

    python3 scripts/enrich_family_enrichment.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
import xml.sax.saxutils as xml
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
OUT_DIR = ROOT / "data" / "family_enrichment"


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


def split_ids(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [x.strip() for x in str(raw).split(";") if x.strip()]


def label_map(con: sqlite3.Connection) -> dict[str, str]:
    """Best-effort accession → human label from domains_json."""
    out: dict[str, str] = {}
    for (blob,) in con.execute(
        "SELECT domains_json FROM proteins "
        "WHERE domains_json IS NOT NULL AND domains_json != ''"
    ):
        try:
            entries = json.loads(blob)
        except json.JSONDecodeError:
            continue
        if not isinstance(entries, list):
            continue
        for e in entries:
            if not isinstance(e, dict):
                continue
            acc = e.get("accession") or e.get("id")
            lab = e.get("label") or e.get("name") or e.get("description")
            if acc and lab and acc not in out:
                out[str(acc)] = str(lab)
    return out


def enrich_db(
    *,
    db_name: str,
    id_col: str,
    rows: list[tuple],
    fg: set[str],
    N: int,
    K: int,
    min_foreground: int,
    labels: dict[str, str],
) -> list[dict]:
    bg_counts: Counter[str] = Counter()
    fg_counts: Counter[str] = Counter()
    members: dict[str, list[str]] = defaultdict(list)

    for pid, _tier, raw in rows:
        ids = set(split_ids(raw))
        for fam in ids:
            bg_counts[fam] += 1
            if pid in fg:
                fg_counts[fam] += 1
                members[fam].append(pid)

    tested: list[dict] = []
    for fam, n in bg_counts.items():
        k = fg_counts.get(fam, 0)
        if k < min_foreground:
            continue
        fold = (k / K) / (n / N) if K and n and N else 0.0
        pval = hypergeom_sf(k, N, K, n)
        tested.append(
            {
                "db": db_name,
                "family_id": fam,
                "label": labels.get(fam, ""),
                "n_background": n,
                "n_tier_ab": k,
                "frac_tier_ab": round(k / n, 6) if n else 0.0,
                "background_tier_ab_rate": round(K / N, 6) if N else 0.0,
                "enrichment_fold": round(fold, 6),
                "hypergeom_pvalue": pval,
                "example_proteins": ",".join(sorted(members[fam])[:8]),
            }
        )

    pvals = [r["hypergeom_pvalue"] for r in tested]
    fdrs = bh_fdr(pvals)
    for r, fdr in zip(tested, fdrs):
        r["bh_fdr"] = round(fdr, 6)
        r["hypergeom_pvalue"] = round(r["hypergeom_pvalue"], 6)
        r["neg_log10_fdr"] = (
            round(-math.log10(max(r["bh_fdr"], 1e-300)), 4) if r["bh_fdr"] > 0 else 300.0
        )

    tested.sort(
        key=lambda r: (r["bh_fdr"], r["hypergeom_pvalue"], -r["enrichment_fold"], -r["n_tier_ab"])
    )
    for i, r in enumerate(tested, start=1):
        r["rank"] = i
    return tested


def write_svg(pfam_rows: list[dict], path: Path, top_n: int = 20) -> None:
    """Figure-ready horizontal bar chart of top Pfam families by −log10(FDR)."""
    rows = [r for r in pfam_rows if r["bh_fdr"] < 0.05][:top_n]
    if not rows:
        rows = pfam_rows[:top_n]
    width, left, right, top, row_h = 920, 210, 40, 56, 22
    height = top + 28 + row_h * max(len(rows), 1) + 40
    max_val = max((r["neg_log10_fdr"] for r in rows), default=1.0) or 1.0
    bar_w = width - left - right

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="#fbfaf7"/>',
        '<text x="460" y="28" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" '
        'font-size="16" font-weight="700" fill="#1a1a1a">'
        "Tier A∪B Pfam enrichment (−log10 FDR)</text>",
        '<text x="460" y="46" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" '
        'font-size="11" fill="#555">'
        "Hypergeometric vs scored proteome · BH FDR within Pfam tests · exploratory if FDR≥0.05"
        "</text>",
    ]
    for i, r in enumerate(rows):
        y = top + i * row_h
        w = max(2.0, bar_w * (r["neg_log10_fdr"] / max_val))
        color = "#0B6E4F" if r["bh_fdr"] < 0.05 else "#8a8a8a"
        lab = r["family_id"]
        if r["label"]:
            short = r["label"] if len(r["label"]) <= 28 else r["label"][:25] + "…"
            lab = f"{r['family_id']}  {short}"
        parts.append(
            f'<text x="8" y="{y + 14}" font-family="Helvetica,Arial,sans-serif" '
            f'font-size="11" fill="#222">{xml.escape(lab)}</text>'
        )
        parts.append(
            f'<rect x="{left}" y="{y}" width="{w:.1f}" height="16" fill="{color}" opacity="0.85"/>'
        )
        parts.append(
            f'<text x="{left + w + 6:.1f}" y="{y + 13}" '
            f'font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#333">'
            f'fold={r["enrichment_fold"]:.2f}  k={r["n_tier_ab"]}/{r["n_background"]}'
            f'</text>'
        )
    parts.append("</svg>\n")
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--min-foreground", type=int, default=3)
    ap.add_argument("--top-n", type=int, default=50)
    args = ap.parse_args()

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row
    tiers = {
        r["protein_id"]: r["tier"]
        for r in con.execute("SELECT protein_id, tier FROM proteins")
    }
    fg = {pid for pid, t in tiers.items() if t in ("A", "B")}
    N = len(tiers)
    K = len(fg)
    labels = label_map(con)

    pfam_src = [
        (r["protein_id"], r["tier"], r["domain_pfam_ids"])
        for r in con.execute(
            "SELECT protein_id, tier, domain_pfam_ids FROM proteins"
        )
    ]
    ipr_src = [
        (r["protein_id"], r["tier"], r["domain_interpro_ids"])
        for r in con.execute(
            "SELECT protein_id, tier, domain_interpro_ids FROM proteins"
        )
    ]

    pfam_rows = enrich_db(
        db_name="Pfam",
        id_col="domain_pfam_ids",
        rows=pfam_src,
        fg=fg,
        N=N,
        K=K,
        min_foreground=args.min_foreground,
        labels=labels,
    )
    ipr_rows = enrich_db(
        db_name="InterPro",
        id_col="domain_interpro_ids",
        rows=ipr_src,
        fg=fg,
        N=N,
        K=K,
        min_foreground=args.min_foreground,
        labels=labels,
    )
    all_rows = pfam_rows + ipr_rows

    fieldnames = [
        "rank",
        "db",
        "family_id",
        "label",
        "n_background",
        "n_tier_ab",
        "frac_tier_ab",
        "background_tier_ab_rate",
        "enrichment_fold",
        "hypergeom_pvalue",
        "bh_fdr",
        "neg_log10_fdr",
        "example_proteins",
    ]
    tsv_path = out_dir / "top_families.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in all_rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    svg_path = out_dir / "top_families_pfam.svg"
    write_svg(pfam_rows, svg_path, top_n=20)
    public_fig = ROOT / "public" / "figures" / "top_families_pfam.svg"
    public_fig.parent.mkdir(parents=True, exist_ok=True)
    public_fig.write_text(svg_path.read_text(encoding="utf-8"), encoding="utf-8")

    honesty = (
        "Hypergeometric enrichment of Pfam/InterPro families among PlantLLPS Tier A∪B "
        f"vs the scored proteome (N={N}, K={K}). BH FDR controlled within each database. "
        "FDR<0.05 hits are the primary claim; remaining ranks are exploratory domain bias, "
        "not causal condensate membership. Labels from deposited domains_json when present."
    )

    summary = {
        "version": "1.5.0",
        "layer": "family_enrichment",
        "n_background": N,
        "n_tier_ab": K,
        "background_tier_ab_rate": round(K / N, 6) if N else 0.0,
        "min_foreground": args.min_foreground,
        "n_pfam_tested": len(pfam_rows),
        "n_interpro_tested": len(ipr_rows),
        "n_pfam_fdr_lt_0.05": sum(1 for r in pfam_rows if r["bh_fdr"] < 0.05),
        "n_interpro_fdr_lt_0.05": sum(1 for r in ipr_rows if r["bh_fdr"] < 0.05),
        "top_pfam": pfam_rows[:15],
        "top_interpro": ipr_rows[:15],
        "paths": {
            "top_families": str(tsv_path.relative_to(ROOT)),
            "figure": str(svg_path.relative_to(ROOT)),
            "ui_tables": "data/family_enrichment/ui_tables.json",
        },
        "honesty": honesty,
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    ui = {
        "parameters": {
            "n_background": N,
            "n_tier_ab": K,
            "background_tier_ab_rate": round(K / N, 6) if N else 0.0,
            "min_foreground": args.min_foreground,
        },
        "honesty": honesty,
        "pfam": pfam_rows[: args.top_n],
        "interpro": ipr_rows[: args.top_n],
        "figure": "/figures/top_families_pfam.svg",
    }
    (out_dir / "ui_tables.json").write_text(
        json.dumps(ui, indent=2) + "\n", encoding="utf-8"
    )

    (out_dir / "README.md").write_text(
        f"""# Tier A∪B Pfam / InterPro family enrichment

Hypergeometric enrichment of domain families among PlantLLPS **Tier A∪B** vs the
scored proteome background.

## Reproduce

```bash
python3 scripts/enrich_family_enrichment.py
```

## Deposits

| File | Contents |
|------|----------|
| `top_families.tsv` | Ranked Pfam + InterPro tests (fold, p, BH FDR) |
| `top_families_pfam.svg` | Figure-ready top Pfam −log10(FDR) bars |
| `ui_tables.json` | Compact payload for Node `/analyses` |
| `summary.json` | Run summary + honesty |

## Headline

- Background **{N:,}**; Tier A∪B **{K:,}** (rate {K / N:.4f})
- Pfam tested **{len(pfam_rows)}** (FDR&lt;0.05: **{summary["n_pfam_fdr_lt_0.05"]}**)
- InterPro tested **{len(ipr_rows)}** (FDR&lt;0.05: **{summary["n_interpro_fdr_lt_0.05"]}**)

## Honesty

{honesty}
""",
        encoding="utf-8",
    )

    meta = {k: summary[k] for k in summary if k not in ("top_pfam", "top_interpro")}
    meta["top_pfam"] = [
        {
            "family_id": r["family_id"],
            "label": r["label"],
            "n_tier_ab": r["n_tier_ab"],
            "n_background": r["n_background"],
            "enrichment_fold": r["enrichment_fold"],
            "bh_fdr": r["bh_fdr"],
        }
        for r in pfam_rows[:8]
    ]
    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("family_enrichment", json.dumps(meta)),
    )
    con.commit()
    con.close()

    print(json.dumps({k: summary[k] for k in summary if k not in ("top_pfam", "top_interpro")}, indent=2))
    print("Top Pfam:", flush=True)
    for r in pfam_rows[:8]:
        print(
            f"  {r['family_id']}: k={r['n_tier_ab']}/{r['n_background']} "
            f"fold={r['enrichment_fold']} FDR={r['bh_fdr']} {r['label'][:40]}",
            flush=True,
        )
    print(f"Wrote {out_dir}", flush=True)


if __name__ == "__main__":
    main()
