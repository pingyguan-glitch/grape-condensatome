# Structure–score joint (PlantLLPS × low-pLDDT)

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

- AF-covered **23,421** / 25,972; Tier A∪B with AF **3,423**
- Joint flags (af_frac ≥ 0.25): **1,872**
- Enrichment tests FDR&lt;0.05: **4** / 4

## Honesty

Joint overlap of PlantLLPS Tier A∪B with AlphaFold low-pLDDT (af_frac_very_low ≥ 0.25 for per-protein flags). Hypergeometric tests across discrete fraction thresholds are BH-controlled; the ranked flag table is exploratory co-occurrence, not proof of LLPS. AF coverage incomplete (23421/25972 proteins).
