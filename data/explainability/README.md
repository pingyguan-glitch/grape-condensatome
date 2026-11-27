# PlantLLPS explainability deposits

SHAP-equivalent **Ser-rich / charge-block / IDR-disorder** attributions for
GrapeCondensatome proteins in the curated LLPS gold overlay **or** score **Tier A**
(always including anchor **F6HET4**).

## Method

- **Logistic head:** exact linear SHAP in logit space (`φ_i = w_i · z_i`).
- **Companion composite:** transparent weighted term decomposition
  (disorder 0.55 / LC 0.25 / blockiness 0.20 + MAP65/Ser bonuses).
- **Hybrid grammar axes:** average of logistic probability-local SHAP (grouped)
  and companion terms mapped onto the same axes.
- **Residue windows:** Gly occlusion ΔPlantLLPS (SHAP-equivalent leave-one-window-out);
  windows tagged `ser_rich` / `charge_block` when composition warrants; joined to
  deposited IDR / domain intervals.

Not KernelSHAP / TreeSHAP — the logistic model admits an exact closed form.
Occlusion is the sequence-space equivalent for local drivers.

## Files

| Path | Notes |
|------|-------|
| `attributions.json.gz` | All gold ∪ Tier A payloads (compact) |
| `F6HET4.json` | Pretty NAR+ example for VvMAP65-1 |
| `summary.json` | Counts + F6HET4 axis snapshot |
| `README.md` | This file |

## Coverage

- Proteins deposited: **1159**
- Gold: **34** · Tier A: **1142** · Union: **1159**
- With occlusion windows: **1159**
- Model: `/agent/repos/grape-condensatome/data/model_formal.json` · window=25 step=12

## Anchor F6HET4 (VvMAP65-1)

- Tier **B** · PlantLLPS ≈ **0.761767** · gold=1
- Grammar hybrid axes: Ser-rich ≈ **0.133** ·
  charge-block ≈ **0.121201** ·
  IDR/disorder ≈ **-0.083465**

Regenerate::

```bash
python3 scripts/deposit_explainability.py
```

UI: protein pages load `attributions.json.gz` at startup (gold ∪ Tier A); live
`/score` uses the JS SHAP port for feature/grammar bars without occlusion.
