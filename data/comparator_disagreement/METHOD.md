# Comparator disagreement method note

## Definition

**High PlantLLPS:** Tier A∪B on the deposited grape condensatome
(`n` = 3870).

**Low comparator:** native proteome columns from
`scripts/score_psap_proteome.py` / `scripts/score_llphyscore_proteome.py`:

- PSAP ≤ proteome quantile **0.25** (threshold = **0.01083333**)
- LLPhyScore percentile ≤ proteome quantile **0.25**
  (threshold = **33.249816**)

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
