# Release notes

## v0.3.1 (2026-11-26)

Submission-readiness follow-on (base: v0.3.0 / PR #1 branch).

### metapredict

- Optional `--use-metapredict` on `score` / `score-proteome` (+ `METAPREDICT_PYTHON`).
- `plantllps/metapredict_bridge.py` + `scripts/score_metapredict_subset.py` (29 sequences scored; mean frac ≈ 0.35).
- Install / pin notes in `docs/EXTERNAL_TOOLS.md`.

### Native PSAP

- Runnable in Python 3.9 env (`numpy<1.24`, `scipy<1.11`); `scripts/run_psap_gold.py`.
- Gold-standard AUC ≈ **0.718** (companion cited 0.692). Side-by-side in `data/benchmarks/gold_standard_comparison.json`.
- MAP-screen ρ still blocked (208-protein table missing).

### Expression Atlas

- `scripts/enrich_expression_atlas.py` joins **E-MTAB-8758** TPMs → SQLite (`expr_*` columns).
- **25 810 / 25 972** proteins matched (99.4%); Node + Python UIs show mean TPM + per-condition.

### Hosting

- `docs/HOSTING.md`: Cursor publish-app on `:8765`; GitHub Pages unsuitable; Railway/Render/Fly alternatives.

## v0.3.0 (2026-11-26)

Submission-oriented resource packaging (still composition-proxy disorder; not metapredict).

### Model / training

- **v13** training set: v12 + **19 curated disordered non-LLPS hard negatives** (Arabidopsis extensins / PRPs / AGPs / LRX / linker histones + human prothymosin-α).
- Canonical files: `data/train_formal_v13.csv` (mirrored as `data/train_formal.csv`), `data/model_formal.json`.
- Metric delta vs v12 (stratified 5-fold, same code path):

| Metric | v12 | v13 | Δ |
|--------|----:|----:|--:|
| n (pos/neg) | 143 (76/67) | 162 (76/86) | +19 neg |
| mean AUC logistic | 0.859 | 0.863 | +0.004 |
| OOF AUC logistic | 0.853 | 0.868 | +0.015 |
| mean AUC companion | 0.861 | 0.774 | −0.088 |

Companion composite AUC drops as expected when highly disordered non-LLPS proteins enter the negative class; logistic holds / slightly improves — the intended hard-negative effect.

### Benchmarks

- Reconstructed companion **gold standard** (13 pos / 15 neg); scored with PlantLLPS + **LLPhyScore** standalone 20220627.
- LLPhyScore gold AUC ≈ **0.769** (companion cited 0.774).
- **PSAP** not runnable on Python 3.12 in this environment (sklearn-json / sklearn~=0.21); companion PSAP gold AUC **0.692** and MAP-screen ρ **0.368** cited, not re-run.
- Grapevine MAP screen table with PSAP columns **absent** from store/repo — MAP ρ not recomputed.

### Atlas / UI

- Proteome re-scored with v13: **25 972** proteins; tiers A=1142, B=2728, C=9653, D=12449.
- UniProt gene names / GO IDs / keywords joined for all scored rows.
- Node + Python web UIs show gene, protein name, GO, keywords; search includes those fields.
- Anchor: VvMAP65-1 `F6HET4` / `VIT_19s0090g01230` remains **Tier B** (score ≈ 0.762).

### Packaging

- Dual deliverable: Python `plantllps` CLI + Express browser.
- `LICENSE` (MIT), `CITATION.cff`, expanded `README.md`.

## v0.2.0 (prior)

Web app publish snapshot with v12 proteome scores (Node-only repo start).
