# Tier A∪B Pfam / InterPro family enrichment

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

- Background **25,972**; Tier A∪B **3,870** (rate 0.1490)
- Pfam tested **424** (FDR&lt;0.05: **144**)
- InterPro tested **984** (FDR&lt;0.05: **323**)

## Honesty

Hypergeometric enrichment of Pfam/InterPro families among PlantLLPS Tier A∪B vs the scored proteome (N=25972, K=3870). BH FDR controlled within each database. FDR<0.05 hits are the primary claim; remaining ranks are exploratory domain bias, not causal condensate membership. Labels from deposited domains_json when present.
