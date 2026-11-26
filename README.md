# GrapeCondensatome

Browse PlantLLPS scores for the *Vitis vinifera* UniProt reference proteome (UP000009183).

## Run

```bash
npm install
npm start
# → http://127.0.0.1:8765/
```

Env: `PORT` (default 8765), `HOST` (default `0.0.0.0`), optional `GRAPE_DB` / `GRAPE_MODEL`.

Requires **Node ≥ 22** (`node:sqlite`).

## Features

- Search by UniProt ID / keyword
- Browse tiers A–D with pagination
- Protein detail (scores, percentile, features, map65_tail_like)
- FASTA paste → score with bundled `data/model_formal.json`

## Data

Bundled under `data/`:

- `grape_condensatome.sqlite` (~22 MB, 25 972 proteins)
- `model_formal.json` (PlantLLPS v12 weights)
