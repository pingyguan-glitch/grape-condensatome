#!/usr/bin/env python3
"""Deposit reconstructed grapevine MAP screen + Spearman ρ vs PlantLLPS / PSAP.

The companion 208-protein table with original PSAP columns was never present in
this repository or the project store. This script reconstructs an
annotation-based MAP set from the atlas UniProt joins and recomputes Spearman ρ
against **native proteome PSAP / LLPhyScore** and PlantLLPS scores already in
SQLite.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
from pathlib import Path
from typing import Iterable, List, Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
OUT_DIR = ROOT / "data" / "map_screen"

import sys

sys.path.insert(0, str(ROOT))
from plantllps.metrics import spearman_rho  # noqa: E402


def _finite_pairs(xs: Sequence[float], ys: Sequence[float]):
    a, b = [], []
    for x, y in zip(xs, ys):
        if x is None or y is None:
            continue
        fx, fy = float(x), float(y)
        if fx == fx and fy == fy:
            a.append(fx)
            b.append(fy)
    return a, b


def spearman_with_n(xs, ys):
    a, b = _finite_pairs(xs, ys)
    if len(a) < 3:
        return {"n": len(a), "rho": float("nan")}
    return {"n": len(a), "rho": spearman_rho(a, b)}


def p_value_approx(rho: float, n: int) -> float:
    """Two-sided p via Student-t approximation (large-n; stdlib only)."""
    if n < 4 or rho != rho or abs(rho) >= 1.0:
        return float("nan")
    t = rho * math.sqrt((n - 2) / max(1e-12, 1 - rho * rho))
    # crude survival using erfc for normal approx of t when n large
    # better: incomplete beta; use normal approx for reporting honesty
    z = abs(t)
    # Abramowitz-style erfc approx via math.erfc
    return math.erfc(z / math.sqrt(2.0))


def fetch_rows(con: sqlite3.Connection) -> List[dict]:
    con.row_factory = sqlite3.Row
    sql = """
    SELECT protein_id, gene_name, protein_name, go_terms, keywords, length, tier,
           score, score_companion, score_logistic, score_psap,
           score_llphyscore_percentile, frac_disorder_metapredict,
           af_frac_very_low, af_mean_plddt, cross_crop_family, features_json
    FROM proteins
    WHERE go_terms LIKE '%GO:0008017%'
       OR keywords LIKE '%Microtubule%'
    ORDER BY protein_id
    """
    rows = []
    for r in con.execute(sql):
        d = dict(r)
        memb = []
        if d.get("go_terms") and "GO:0008017" in d["go_terms"]:
            memb.append("GO:0008017")
        if d.get("keywords") and "Microtubule" in d["keywords"]:
            memb.append("KW-Microtubule")
        d["membership"] = ";".join(memb)
        d["set_class"] = "annotation"
        rows.append(d)
    # name-rescue extras (companion expanded 208→256)
    seen = {r["protein_id"] for r in rows}
    rescue_sql = """
    SELECT protein_id, gene_name, protein_name, go_terms, keywords, length, tier,
           score, score_companion, score_logistic, score_psap,
           score_llphyscore_percentile, frac_disorder_metapredict,
           af_frac_very_low, af_mean_plddt, cross_crop_family, features_json
    FROM proteins WHERE
      lower(IFNULL(protein_name,'')) LIKE '%map65%'
      OR lower(IFNULL(protein_name,'')) LIKE '%map70%'
      OR lower(IFNULL(protein_name,'')) LIKE '%wave-dampened%'
      OR lower(IFNULL(protein_name,'')) LIKE '%wvd2%'
      OR lower(IFNULL(protein_name,'')) LIKE '%spiral1%'
      OR lower(IFNULL(protein_name,'')) LIKE '%tpx2%'
      OR lower(IFNULL(protein_name,'')) LIKE '%iq67%'
      OR lower(IFNULL(protein_name,'')) LIKE '%iqd%'
      OR lower(IFNULL(protein_name,'')) LIKE '%kinesin%'
      OR lower(IFNULL(protein_name,'')) LIKE '%microtubule-associated%'
      OR lower(IFNULL(protein_name,'')) LIKE '%microtubule associated%'
      OR IFNULL(cross_crop_family,'') != ''
    """
    for r in con.execute(rescue_sql):
        if r["protein_id"] in seen:
            continue
        d = dict(r)
        d["membership"] = "name_rescue"
        d["set_class"] = "name_rescue"
        rows.append(d)
        seen.add(r["protein_id"])
    return rows


COLS = [
    "protein_id",
    "gene_name",
    "protein_name",
    "set_class",
    "membership",
    "length",
    "tier",
    "score_plantllps",
    "score_companion",
    "score_logistic",
    "score_psap",
    "score_llphyscore_percentile",
    "frac_disorder_metapredict",
    "af_frac_very_low",
    "af_mean_plddt",
    "cross_crop_family",
]


def write_tsv(path: Path, rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            out = {
                "protein_id": r["protein_id"],
                "gene_name": r.get("gene_name") or "",
                "protein_name": r.get("protein_name") or "",
                "set_class": r.get("set_class") or "",
                "membership": r.get("membership") or "",
                "length": r.get("length"),
                "tier": r.get("tier") or "",
                "score_plantllps": r.get("score"),
                "score_companion": r.get("score_companion"),
                "score_logistic": r.get("score_logistic"),
                "score_psap": r.get("score_psap"),
                "score_llphyscore_percentile": r.get("score_llphyscore_percentile"),
                "frac_disorder_metapredict": r.get("frac_disorder_metapredict"),
                "af_frac_very_low": r.get("af_frac_very_low"),
                "af_mean_plddt": r.get("af_mean_plddt"),
                "cross_crop_family": r.get("cross_crop_family") or "",
            }
            w.writerow(out)


def metrics_for(rows: List[dict], label: str) -> dict:
    comp = [r.get("score_companion") for r in rows]
    hyb = [r.get("score") for r in rows]
    psap = [r.get("score_psap") for r in rows]
    llphy = [r.get("score_llphyscore_percentile") for r in rows]
    m = {
        "set": label,
        "n": len(rows),
        "n_go8017": sum(1 for r in rows if "GO:0008017" in (r.get("membership") or "")),
        "n_kw_microtubule": sum(
            1 for r in rows if "KW-Microtubule" in (r.get("membership") or "")
        ),
        "spearman_companion_vs_psap": spearman_with_n(comp, psap),
        "spearman_plantllps_vs_psap": spearman_with_n(hyb, psap),
        "spearman_companion_vs_llphyscore": spearman_with_n(comp, llphy),
        "spearman_plantllps_vs_llphyscore": spearman_with_n(hyb, llphy),
    }
    for key in (
        "spearman_companion_vs_psap",
        "spearman_plantllps_vs_psap",
        "spearman_companion_vs_llphyscore",
        "spearman_plantllps_vs_llphyscore",
    ):
        block = m[key]
        block["p_approx_normal"] = p_value_approx(block["rho"], block["n"])
    # VvMAP65-1 ranks
    with_psap = [r for r in rows if r.get("score_psap") is not None]
    with_psap.sort(key=lambda r: -float(r["score_psap"]))
    with_ll = [r for r in rows if r.get("score_llphyscore_percentile") is not None]
    with_ll.sort(key=lambda r: -float(r["score_llphyscore_percentile"]))
    def rank_of(sorted_rows, pid="F6HET4"):
        for i, r in enumerate(sorted_rows, 1):
            if r["protein_id"] == pid:
                return i
        return None
    m["F6HET4"] = {
        "psap_rank": rank_of(with_psap),
        "psap_n": len(with_psap),
        "llphyscore_rank": rank_of(with_ll),
        "llphyscore_n": len(with_ll),
        "score_psap": next((r["score_psap"] for r in rows if r["protein_id"] == "F6HET4"), None),
        "score_llphyscore_percentile": next(
            (r["score_llphyscore_percentile"] for r in rows if r["protein_id"] == "F6HET4"),
            None,
        ),
        "score_plantllps": next((r["score"] for r in rows if r["protein_id"] == "F6HET4"), None),
    }
    return m


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    con = sqlite3.connect(args.db)
    rows = fetch_rows(con)
    ann = [r for r in rows if r["set_class"] == "annotation"]
    go_only = [r for r in ann if "GO:0008017" in (r.get("membership") or "")]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out_dir / "map_screen_annotation.tsv", ann)
    write_tsv(args.out_dir / "map_screen_extended.tsv", rows)
    write_tsv(args.out_dir / "map_screen_go8017.tsv", go_only)

    report = {
        "version": 1,
        "provenance": {
            "companion_cited": {
                "n_go_annotated": 208,
                "n_name_rescued_extra": 48,
                "n_extended": 256,
                "spearman_composite_vs_psap": 0.368,
                "spearman_p": 4.6e-8,
                "spearman_composite_vs_llphyscore": 0.232,
                "note": (
                    "Companion Supplemental Data S1 table with original PSAP columns "
                    "was never deposited in this repo or project store."
                ),
            },
            "this_deposit": {
                "method": (
                    "Reconstructed from UniProt joins in grape_condensatome.sqlite: "
                    "GO:0008017 ∪ keyword Microtubule (= annotation set), plus modest "
                    "name-rescue / cross_crop_family extras (= extended set). "
                    "PSAP / LLPhyScore columns are **native proteome** scores from "
                    "scripts/score_psap_proteome.py and score_llphyscore_proteome.py, "
                    "not the companion's original MAP-screen PSAP run."
                ),
                "annotation_drift": (
                    "Current UP000009183 annotations yield GO:0008017 n≠208; "
                    "do not treat set membership as identical to the companion 208."
                ),
            },
        },
        "metrics": {
            "annotation": metrics_for(ann, "annotation_GO8017_or_KW_Microtubule"),
            "go8017_only": metrics_for(go_only, "GO:0008017_only"),
            "extended": metrics_for(rows, "annotation_plus_name_rescue"),
        },
        "honesty": (
            "Recomputed ρ uses PlantLLPS companion/hybrid scores (composition-proxy "
            "disorder + MAP65 grammar bonus) versus native proteome PSAP. These are "
            "**not** a byte-for-byte replication of companion ρ=0.368. Cite companion "
            "ρ when discussing the published transfer gap; cite this deposit when "
            "discussing the public atlas reconstruction."
        ),
    }
    (args.out_dir / "map_screen_spearman.json").write_text(json.dumps(report, indent=2))
    readme = f"""# Grapevine MAP screen deposit (reconstructed)

## Why this exists

The companion MAP/LLPS manuscript cites Spearman ρ = **0.368** (PSAP vs plant
composite) on **208** GO-annotated grapevine microtubule-binding proteins.
That table (with original PSAP columns) was **not** present in the public
repository or project store, so ρ could not be recomputed from a deposited file.

## What is deposited here

| File | Contents |
|------|----------|
| `map_screen_go8017.tsv` | Current atlas proteins with GO:0008017 |
| `map_screen_annotation.tsv` | GO:0008017 **∪** UniProt keyword `Microtubule` |
| `map_screen_extended.tsv` | Annotation set + name-rescue / cross-crop hooks |
| `map_screen_spearman.json` | Recomputed Spearman ρ + provenance |

## Recomputed metrics (annotation set)

- n = **{report['metrics']['annotation']['n']}** (GO:0008017 n={report['metrics']['annotation']['n_go8017']})
- ρ(companion vs native PSAP) ≈ **{report['metrics']['annotation']['spearman_companion_vs_psap']['rho']:.3f}**
- ρ(PlantLLPS hybrid vs native PSAP) ≈ **{report['metrics']['annotation']['spearman_plantllps_vs_psap']['rho']:.3f}**
- ρ(companion vs LLPhyScore %ile) ≈ **{report['metrics']['annotation']['spearman_companion_vs_llphyscore']['rho']:.3f}**
- F6HET4 native PSAP rank ≈ **{report['metrics']['annotation']['F6HET4']['psap_rank']}** / {report['metrics']['annotation']['F6HET4']['psap_n']}

## Honesty

- Companion-cited ρ = **0.368** remains the manuscript transfer-gap number.
- This deposit recomputes ρ on an **annotation-reconstructed** set with **native
  proteome PSAP** columns — expected to differ from 0.368 (annotation drift,
  composition-proxy vs AF disorder composite, PSAP run differences).
- Regenerate: `python3 scripts/deposit_map_screen.py`
"""
    (args.out_dir / "README.md").write_text(readme)
    # also mirror summary into benchmarks + meta
    bench = ROOT / "data" / "benchmarks" / "map_screen_spearman.json"
    bench.write_text(json.dumps(report, indent=2))
    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("map_screen_spearman", json.dumps(report["metrics"])),
    )
    con.commit()
    print(json.dumps(report["metrics"], indent=2))
    print(f"Wrote {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
