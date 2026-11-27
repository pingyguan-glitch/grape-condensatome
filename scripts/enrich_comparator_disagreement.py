#!/usr/bin/env python3
"""Comparator disagreement set: PlantLLPS high / PSAP and/or LLPhyScore low.

Defines the set of Tier A∪B proteins whose native proteome PSAP and/or
LLPhyScore percentile fall in the bottom quartile of the scored proteome
(Ser-rich / MAP-class focus flagged). Deposits a ranked list + method note
under ``data/comparator_disagreement/`` (SQLite meta pointer only).

Honesty
-------
- Thresholds are atlas-relative quartiles on deposited native comparator columns,
  not published cutoffs from PSAP/LLPhyScore papers.
- Ser-rich / MAP focus uses PlantLLPS features (``map65_tail_like``, ``frac_ser``,
  ``c_term_ser_frac_100``) plus optional ``cross_crop_family`` MAP tags — hypothesis
  generating, not an experimental MAP screen call.
- Disagreement is exploratory unless interpreted together with FDR-controlled
  family enrichment / gold benchmarks elsewhere in the resource.

Example::

    python3 scripts/enrich_comparator_disagreement.py
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
OUT_DIR = ROOT / "data" / "comparator_disagreement"


def percentile_rank(sorted_vals: list[float], x: float) -> float:
    """Empirical CDF percentile in [0, 100] for value x among sorted_vals."""
    if not sorted_vals:
        return float("nan")
    # fraction of values strictly < x, plus half of ties
    lo = 0
    hi = len(sorted_vals)
    while lo < hi:
        mid = (lo + hi) // 2
        if sorted_vals[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    n_lt = lo
    hi = len(sorted_vals)
    while lo < hi:
        mid = (lo + hi) // 2
        if sorted_vals[mid] <= x:
            lo = mid + 1
        else:
            hi = mid
    n_le = lo
    n_eq = n_le - n_lt
    return 100.0 * (n_lt + 0.5 * n_eq) / len(sorted_vals)


def quantile(sorted_vals: list[float], q: float) -> float:
    if not sorted_vals:
        return float("nan")
    if q <= 0:
        return sorted_vals[0]
    if q >= 1:
        return sorted_vals[-1]
    idx = (len(sorted_vals) - 1) * q
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return sorted_vals[lo]
    w = idx - lo
    return sorted_vals[lo] * (1 - w) + sorted_vals[hi] * w


def is_ser_map(
    feats: dict, cross_crop: str | None, protein_name: str | None, gold: int
) -> tuple[bool, str]:
    """Stricter Ser-rich / MAP-class focus (map65_tail_like alone is too common)."""
    reasons = []
    fs = float(feats.get("frac_ser") or 0)
    ct = float(feats.get("c_term_ser_frac_100") or 0)
    if fs >= 0.12:
        reasons.append("frac_ser>=0.12")
    if ct >= 0.15:
        reasons.append("c_term_ser>=0.15")
    cc = (cross_crop or "").upper()
    if "MAP" in cc:
        reasons.append("cross_crop_MAP")
    pn = (protein_name or "").lower()
    if "microtubule-associated" in pn:
        reasons.append("name_MAP")
    if gold:
        reasons.append("llps_gold")
    return (bool(reasons), ",".join(reasons))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument(
        "--low-quantile",
        type=float,
        default=0.25,
        help="Proteome quantile for 'low' PSAP / LLPhyScore (default 0.25)",
    )
    args = ap.parse_args()

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(str(args.db))
    rows = list(
        con.execute(
            """
            SELECT protein_id, gene_name, protein_name, tier, score,
                   score_psap, score_llphyscore_percentile, features_json,
                   cross_crop_family, llps_gold
            FROM proteins
            """
        )
    )

    psap_all = sorted(float(r[5]) for r in rows if r[5] is not None)
    llphy_all = sorted(float(r[6]) for r in rows if r[6] is not None)
    psap_lo = quantile(psap_all, args.low_quantile)
    llphy_lo = quantile(llphy_all, args.low_quantile)

    records = []
    for (
        pid,
        gene,
        pname,
        tier,
        score,
        psap,
        llphy,
        feats_raw,
        cross_crop,
        gold,
    ) in rows:
        if tier not in ("A", "B"):
            continue
        if psap is None and llphy is None:
            continue
        feats = {}
        if feats_raw:
            try:
                feats = json.loads(feats_raw)
            except json.JSONDecodeError:
                feats = {}
        ser_map, ser_reasons = is_ser_map(feats, cross_crop, pname, int(gold or 0))
        low_psap = psap is not None and float(psap) <= psap_lo
        low_llphy = llphy is not None and float(llphy) <= llphy_lo
        if not (low_psap or low_llphy):
            continue
        psap_pct = percentile_rank(psap_all, float(psap)) if psap is not None else float("nan")
        # LLPhyScore column is already a proteome percentile (0–100)
        llphy_pct = float(llphy) if llphy is not None else float("nan")
        plant_pct = float(feats.get("percentile") or float("nan"))
        # disagreement: high PlantLLPS percentile vs low native comparator percentiles
        comp_floor = []
        if low_psap and math.isfinite(psap_pct):
            comp_floor.append(psap_pct)
        if low_llphy and math.isfinite(llphy_pct):
            comp_floor.append(llphy_pct)
        disagree_gap = (
            plant_pct - min(comp_floor) if math.isfinite(plant_pct) and comp_floor else float("nan")
        )
        records.append(
            {
                "protein_id": pid,
                "gene_name": gene or "",
                "protein_name": pname or "",
                "tier": tier,
                "score_plantllps": round(float(score), 8),
                "plantllps_percentile": round(plant_pct, 4) if math.isfinite(plant_pct) else "",
                "score_psap": round(float(psap), 8) if psap is not None else "",
                "psap_percentile": round(psap_pct, 4) if math.isfinite(psap_pct) else "",
                "score_llphyscore_percentile": round(float(llphy), 8)
                if llphy is not None
                else "",
                "low_psap": int(low_psap),
                "low_llphyscore": int(low_llphy),
                "low_both": int(low_psap and low_llphy),
                "ser_map_focus": int(ser_map),
                "ser_map_reasons": ser_reasons,
                "frac_ser": round(float(feats.get("frac_ser") or 0), 6),
                "c_term_ser_frac_100": round(
                    float(feats.get("c_term_ser_frac_100") or 0), 6
                ),
                "map65_tail_like": round(float(feats.get("map65_tail_like") or 0), 6),
                "cross_crop_family": cross_crop or "",
                "llps_gold": int(gold or 0),
                "disagreement_gap": round(disagree_gap, 4)
                if math.isfinite(disagree_gap)
                else "",
            }
        )

    # Rank: both-low first, then ser/MAP focus, then gap / PlantLLPS score
    records.sort(
        key=lambda r: (
            -r["low_both"],
            -r["ser_map_focus"],
            -(r["disagreement_gap"] if r["disagreement_gap"] != "" else -1),
            -r["score_plantllps"],
        )
    )
    for i, r in enumerate(records, start=1):
        r["rank"] = i

    fieldnames = [
        "rank",
        "protein_id",
        "gene_name",
        "protein_name",
        "tier",
        "score_plantllps",
        "plantllps_percentile",
        "score_psap",
        "psap_percentile",
        "score_llphyscore_percentile",
        "low_psap",
        "low_llphyscore",
        "low_both",
        "ser_map_focus",
        "ser_map_reasons",
        "frac_ser",
        "c_term_ser_frac_100",
        "map65_tail_like",
        "cross_crop_family",
        "llps_gold",
        "disagreement_gap",
    ]
    tsv_path = out_dir / "disagreement_ranked.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in records:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    n_ab = sum(1 for r in rows if r[3] in ("A", "B"))
    n_ser = sum(1 for r in records if r["ser_map_focus"])
    n_both = sum(1 for r in records if r["low_both"])
    n_psap = sum(1 for r in records if r["low_psap"])
    n_llphy = sum(1 for r in records if r["low_llphyscore"])
    n_gold = sum(1 for r in records if r["llps_gold"])

    method = f"""# Comparator disagreement method note

## Definition

**High PlantLLPS:** Tier A∪B on the deposited grape condensatome
(`n` = {n_ab}).

**Low comparator:** native proteome columns from
`scripts/score_psap_proteome.py` / `scripts/score_llphyscore_proteome.py`:

- PSAP ≤ proteome quantile **{args.low_quantile:.2f}** (threshold = **{psap_lo:.8f}**)
- LLPhyScore percentile ≤ proteome quantile **{args.low_quantile:.2f}**
  (threshold = **{llphy_lo:.6f}**)

A protein enters the disagreement set if it is Tier A∪B and low on **either**
comparator (`low_psap` or `low_llphyscore`). `low_both` marks the stricter
intersection.

## Ser-rich / MAP-class focus

Flag `ser_map_focus` when any of:

- `frac_ser` ≥ 0.12
- `c_term_ser_frac_100` ≥ 0.15
- `cross_crop_family` contains MAP
- protein name contains “microtubule-associated”
- curated `llps_gold`

`map65_tail_like` is reported as a column but **not** used alone for the focus
flag (it is common among Tier A∪B). Thresholds are atlas heuristics for the
companion MAP65 / Ser-rich grammar story (e.g. F6HET4), not experimental
MAP-screen calls.

## Ranking

1. `low_both` first
2. then `ser_map_focus`
3. then `disagreement_gap` = PlantLLPS feature percentile − min(low comparator percentiles)
4. then PlantLLPS score

## Honesty

Quartile cuts are **atlas-relative** on deposited native scores — not published
PSAP/LLPhyScore decision thresholds. The set is **exploratory** (hypothesis
generating for plant Ser-rich / MAP rescue). Prefer FDR-controlled family
enrichment and leakage-controlled gold benchmarks for primary claims.
"""
    (out_dir / "METHOD.md").write_text(method, encoding="utf-8")

    honesty = (
        "Tier A∪B proteins with native PSAP and/or LLPhyScore percentile in the "
        f"bottom proteome quartile (q={args.low_quantile}). Ser-rich/MAP focus = "
        "frac_ser≥0.12 or C-term Ser≥0.15 or MAP annotation/name or gold "
        "(map65_tail_like alone is not sufficient). Exploratory disagreement set — "
        "not an FDR-controlled claim."
    )

    summary = {
        "version": "1.5.0",
        "layer": "comparator_disagreement",
        "n_tier_ab": n_ab,
        "low_quantile": args.low_quantile,
        "psap_low_threshold": round(psap_lo, 8),
        "llphyscore_low_threshold": round(llphy_lo, 6),
        "n_disagreement": len(records),
        "n_low_psap": n_psap,
        "n_low_llphyscore": n_llphy,
        "n_low_both": n_both,
        "n_ser_map_focus": n_ser,
        "n_llps_gold": n_gold,
        "anchor_F6HET4": next((r for r in records if r["protein_id"] == "F6HET4"), None),
        "paths": {
            "ranked": str(tsv_path.relative_to(ROOT)),
            "method": "data/comparator_disagreement/METHOD.md",
            "ui_tables": "data/comparator_disagreement/ui_tables.json",
        },
        "honesty": honesty,
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    ui = {
        "parameters": {
            "n_tier_ab": n_ab,
            "low_quantile": args.low_quantile,
            "psap_low_threshold": round(psap_lo, 8),
            "llphyscore_low_threshold": round(llphy_lo, 6),
            "n_disagreement": len(records),
            "n_low_both": n_both,
            "n_ser_map_focus": n_ser,
        },
        "honesty": honesty,
        "rows": records[:200],
    }
    (out_dir / "ui_tables.json").write_text(
        json.dumps(ui, indent=2) + "\n", encoding="utf-8"
    )

    (out_dir / "README.md").write_text(
        f"""# Comparator disagreement set

PlantLLPS **Tier A∪B** proteins with **low** native PSAP and/or LLPhyScore
(proteome bottom quartile), with Ser-rich / MAP-class focus flags.

## Reproduce

```bash
python3 scripts/enrich_comparator_disagreement.py
```

## Deposits

| File | Contents |
|------|----------|
| `disagreement_ranked.tsv` | Ranked disagreement proteins + flags |
| `METHOD.md` | Short method note (thresholds + honesty) |
| `ui_tables.json` | Compact payload for Node `/analyses` |
| `summary.json` | Run summary |

## Headline

- Tier A∪B **{n_ab:,}**
- Disagreement (either low) **{len(records):,}** · both low **{n_both:,}** · Ser/MAP focus **{n_ser:,}**
- PSAP ≤ {psap_lo:.6g} · LLPhyScore %ile ≤ {llphy_lo:.4g}

## Honesty

{honesty}
""",
        encoding="utf-8",
    )

    meta = {k: summary[k] for k in summary if k != "anchor_F6HET4"}
    if summary["anchor_F6HET4"]:
        a = summary["anchor_F6HET4"]
        meta["anchor_F6HET4"] = {
            "rank": a["rank"],
            "tier": a["tier"],
            "score_plantllps": a["score_plantllps"],
            "score_psap": a["score_psap"],
            "score_llphyscore_percentile": a["score_llphyscore_percentile"],
            "low_psap": a["low_psap"],
            "low_llphyscore": a["low_llphyscore"],
            "ser_map_focus": a["ser_map_focus"],
        }
    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("comparator_disagreement", json.dumps(meta)),
    )
    con.commit()
    con.close()

    print(json.dumps({k: summary[k] for k in summary if k != "anchor_F6HET4"}, indent=2))
    print(f"Wrote {out_dir} ({len(records)} rows)", flush=True)


if __name__ == "__main__":
    main()
