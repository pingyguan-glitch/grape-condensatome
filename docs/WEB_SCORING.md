# Web paste-FASTA scoring

Browser UI for PlantLLPS scoring without leaving the GrapeCondensatome site.

## Demo

```bash
npm install   # once
npm start     # http://127.0.0.1:8765/
```

1. Open [/score](http://127.0.0.1:8765/score).
2. Click **Load sample** (or paste multi-FASTA / upload `.fasta`).
3. Click **Score** — results show the CLI feature vector, hybrid scores, grape-background percentile, composition-proxy disorder track, and grammar SHAP bars.

JSON: `POST /api/score` with `{ "fasta": ">id\\nSEQ..." }`.

Sample input: [`public/examples/score_demo.fasta`](../public/examples/score_demo.fasta).

## Parity with CLI

- Weights: `data/model_formal.json` (v13 hybrid).
- Feature vector: the same names as `plantllps.features.FEATURE_NAMES` (JS port in `lib/features.js`).
- Tool version: root `VERSION` (shown in UI footer, `/api/health`, and `plantllps --version` after package sync).
- Live paste does **not** write rows back into `grape_condensatome.sqlite`.

## Percentiles

- **Grape background:** empirical CDF vs deposited proteome PlantLLPS scores in SQLite:  
  `100 × (# scores ≤ query) / N`.
- **Batch (n>1):** within-paste rank percentile (same spirit as CLI `plantllps score`).

## Caps (browser)

| Cap | Default |
|-----|---------|
| Sequences / request | 25 |
| Length / sequence | 20–8000 aa |
| Total residues | 50 000 |
| Paste characters | 1 500 000 |
| Suggested client timeout | 120 s |

Larger jobs: `python3 -m plantllps score …`.

## Tracks

- Live: composition-proxy disorder trajectory + logistic grammar SHAP.
- Atlas protein pages: deposited metapredict + AF pLDDT dual-rails when present under `data/tracks/` (not recomputed in-browser).

## Docs in UI

- Methods: `/methods`
- Limitations / wet-lab caveats: `/limitations`
