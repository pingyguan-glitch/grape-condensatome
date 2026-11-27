# Grapevine MAP screen deposit (reconstructed)

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

- n = **233** (GO:0008017 n=159)
- ρ(companion vs native PSAP) ≈ **0.202**
- ρ(PlantLLPS hybrid vs native PSAP) ≈ **0.153**
- ρ(companion vs LLPhyScore %ile) ≈ **0.146**
- F6HET4 native PSAP rank ≈ **126** / 233

## Honesty

- Companion-cited ρ = **0.368** remains the manuscript transfer-gap number.
- This deposit recomputes ρ on an **annotation-reconstructed** set with **native
  proteome PSAP** columns — expected to differ from 0.368 (annotation drift,
  composition-proxy vs AF disorder composite, PSAP run differences).
- Regenerate: `python3 scripts/deposit_map_screen.py`
