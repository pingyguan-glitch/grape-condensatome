# Comparator disagreement set

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

- Tier A∪B **3,870**
- Disagreement (either low) **1,555** · both low **277** · Ser/MAP focus **448**
- PSAP ≤ 0.0108333 · LLPhyScore %ile ≤ 33.25

## Honesty

Tier A∪B proteins with native PSAP and/or LLPhyScore percentile in the bottom proteome quartile (q=0.25). Ser-rich/MAP focus = frac_ser≥0.12 or C-term Ser≥0.15 or MAP annotation/name or gold (map65_tail_like alone is not sufficient). Exploratory disagreement set — not an FDR-controlled claim.
