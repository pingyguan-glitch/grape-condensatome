# Data provenance

## LLPSDB

- Downloaded `Phase_separation_unambiguous.zip` and `No_phase_separation_unambiguous.zip` from `http://bio-comp.org.cn/llpsdbv2/download/` (LLPSDB v2).
- Also fetched v1 `LLPS_Natural_protein.zip` from `http://bio-comp.org.cn/llpsdb/download/`.
- Plant taxonomy rows → `llpsdb_v2_plant_entries.json`.

## UniProt proteome

- Query: `proteome:UP000009183` (*Vitis vinifera* reference; 29 861 sequences).
- File: `UP000009183.fasta.gz` via UniProt REST stream API.

## Training / models

- `train_formal.csv` — LLPSDB plant (filtered) + companion gold standards + controls.
- `train_tiny.csv` — earlier 13-row stub (kept for regression).
- `model_formal.json` / `model_formal.metrics.json` — formal train + CV.
- `grape_condensatome.sqlite` — whole-proteome scores (`proteome_scored_v1`).
- `grape_condensatome.prototype.sqlite` — 4-protein seed DB from first prototype.
