# PlantLLPS + GrapeCondensatome

Plant-aware LLPS propensity scoring (**PlantLLPS**) and a browsable *Vitis vinifera* condensatome atlas (**GrapeCondensatome**) for UniProt proteome **UP000009183**.

[![Version](https://img.shields.io/badge/version-1.11.0-0B6E4F)](./VERSION)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)
[![Public atlas](https://img.shields.io/badge/public_atlas-live-0B6E4F)](https://grape-condensatome.apps.grok.bot)
[![Uptime](https://img.shields.io/badge/host-grape--condensatome.apps.grok.bot-0B6E4F)](https://grape-condensatome.apps.grok.bot)
[![Zenodo DOI](https://img.shields.io/badge/DOI-pending_Zenodo_publish-lightgrey)](./docs/ZENODO.md)

**Version:** 1.11.0 · **License:** MIT · **Cite:** see [`CITATION.cff`](./CITATION.cff) · **Public atlas:** [https://grape-condensatome.apps.grok.bot](https://grape-condensatome.apps.grok.bot)

> Hosting: live Cursor publish-app URL recorded in [`docs/HOSTING.md`](./docs/HOSTING.md). Zenodo DOI still pending — do not invent one ([`docs/ZENODO.md`](./docs/ZENODO.md)).

> **v1.0.0** combines Wave 2 curated gold/MAP/hold-out, STRING neighborhoods, and the soybean second-crop condensatome.

Author: Pingyin Guan (College of Horticulture, China Agricultural University; ORCID [0000-0003-4890-3412](https://orcid.org/0000-0003-4890-3412)).

## What ships

| Piece | Path | Notes |
|-------|------|-------|
| Python package / CLI | `plantllps/` | train, score, score-proteome, serve; optional `--use-metapredict` |
| Training history | `data/train_formal_v*.csv`, `data/training/` | v1–v13 snapshots + metrics |
| Model weights | `data/model_formal.json` | **v13** hybrid logistic + companion |
| Proteome SQLite | `data/grape_condensatome.sqlite` | 25 972 scored rows + UniProt + Expression Atlas |
| Proteome TSV | Zenodo bundle / local export | flat table (gitignored; rebuild from SQLite or `export_zenodo_bundle.py`) |
| Node web UI | `server.js` + `lib/` + `public/` | Express + `node:sqlite` (Node ≥22); paste-FASTA `/score` + `/api/score` |
| Web scoring docs | [`docs/WEB_SCORING.md`](./docs/WEB_SCORING.md) | Caps, percentiles, methods/limitations pages, sample FASTA |
| Benchmarks | `data/benchmarks/` | gold vs LLPhyScore + **native PSAP** + family hold-out / transfer + **feature ablation** + **fair NAR+ suite** |
| Fair leaderboard | [`docs/BENCHMARK_LEADERBOARD.md`](./docs/BENCHMARK_LEADERBOARD.md) | leakage-controlled PlantLLPS vs PSAP / LLPhy / PLAAC / ProX-style |
| Experimental models | `data/experimental/` | optional metapredict / AF-pLDDT disorder retrains (**not** production default) |
| Curated LLPS gold | `data/curated/` | experimentally anchored overlay + SQLite/UI badge |
| MAP screen deposit | `data/map_screen/` | annotation-reconstructed MAP table + Spearman ρ |
| External tool docs | [`docs/EXTERNAL_TOOLS.md`](./docs/EXTERNAL_TOOLS.md) | metapredict / PSAP env pins |
| Stress condensatome modules | `data/stress_modules/` | DE × Tier A∪B enrichment + response modules; UI `/modules` |
| Atlas analyses 2/6/7 (v1.5.0) | `data/stress_modules/deepdive/`, `data/string/stress_seed_hubs.tsv`, `data/isoforms/` | Botrytis GO/keyword deep dive + STRING/stress hubs + VIT multi-accession deltas |
| Atlas analyses family/comparator/structure (v1.5.1) | `data/family_enrichment/`, `data/comparator_disagreement/`, `data/structure_score/` | Domain enrichment + comparator disagreement + structure–score; UI `/analyses` |
| Hosting notes | [`docs/HOSTING.md`](./docs/HOSTING.md) | **https://grape-condensatome.apps.grok.bot** (+ local / Render / Railway / Fly fallbacks) |
| Zenodo deposit | [`docs/ZENODO.md`](./docs/ZENODO.md) + `scripts/export_zenodo_bundle.py` | Scoring-server snapshot + sample FASTA + weights + SQLite/TSV; **no invented DOI** |
| Explainability | `data/explainability/` + `plantllps/explain.py` | SHAP-equivalent Ser-rich / charge-block / IDR attributions (gold ∪ Tier A; F6HET4 example) |
| Cross-crop conserved (v1.5.2) | `data/cross_crop/` | Grape↔soy↔Ath Tier A∪B OrthoDB overlap + Ath score deposit |
| Tomato third crop (v1.5.2) | `data/second_crop/tomato_*` | UP000004994 scores + grape `ortholog_sly*`; UI `/hooks` |
| Fruit-crop condensatomes (v1.6.0) | `data/fruit_crops/` | Multi-species fruit-tree PlantLLPS deposits + manifest; UI `/hooks` |
| Fruit cross-species high-score (v1.7.0) | `data/fruit_crops/cross_species/` | Conservative vs species-specific Tier A∪B clusters; UI `/hooks#fruit-cross-species` |
| LLPS DB crosswalk (v1.8.0) | `data/llps_crosswalk/` | PhaSePro / LLPSDB / DrLLPS / CD-CODE ID joins; unique vs shared; UI `/crosswalk` · [`docs/LLPS_CROSSWALK.md`](./docs/LLPS_CROSSWALK.md) |
| RNA-binding crosswalk (v1.9.0) | `data/rbp/` | RBP2GO Ath × OrthoDB + ATtRACT + UniProt GO/KW + RG/RGG; UI `/rbp` |
| Per-protein tracks (v1.10.0) | `data/tracks/` | metapredict + AF pLDDT rails + orthogonal cross-check; grape + candidate deposits |
| NAR Database draft | [`docs/MANUSCRIPT.md`](./docs/MANUSCRIPT.md) | Manuscript synced to software **v1.10.0**; Zenodo DOI pending |
| Supplemental figures S1–S5 (v1.5.3) | [`docs/figures/supp/`](./docs/figures/supp/) | PNG/PDF + [`LEGENDS.md`](./docs/figures/supp/LEGENDS.md); regenerate via `scripts/generate_supp_figures.py` (`requirements-figures.txt`) |

## Honesty (read before citing metrics)

- Default disorder is a **composition proxy** (`data/model_formal.json`). **metapredict** mean fractions via `scripts/score_metapredict_proteome.py`; **per-residue trajectories** via `scripts/deposit_residue_tracks.py` → `data/tracks/` (protein-page sparklines). Optional experimental retrains that *replace* the proxy live under `data/experimental/` (`scripts/retrain_disorder_features.py`) and must not silently replace production weights.
- **Feature-ablation panels** for Figure 3: `scripts/eval_feature_ablation.py` → `data/benchmarks/feature_ablation_*`.
- **Explainability (NAR+):** exact logistic SHAP (`φ=w·z`) + companion decomposition + Gly-occlusion windows via `scripts/deposit_explainability.py` → `data/explainability/` (gold ∪ Tier A; protein-page panel). Not KernelSHAP.
- **AlphaFold pLDDT** summaries + **scored-proteome low-pLDDT segments** via `scripts/enrich_alphafold_plddt.py` (network-bound; resumable cache; browser-like User-Agent).
- **Ath/rice OrthoDB orthologs** + curated `data/cross_crop_hooks.json` via `scripts/enrich_orthologs.py`.
- **Soybean / tomato crop condensatomes** (UP000008827 / UP000004994) via `scripts/score_second_crop_proteome.py` + `scripts/enrich_second_crop_orthologs.py` → `data/second_crop/` + grape `ortholog_gma*` / `ortholog_sly*`; UI `/hooks`.
- **Fruit-crop condensatomes** (apple, pear, peach, cherry, citrus, banana, …) via `scripts/score_fruit_crop_proteomes.py` + optional `scripts/enrich_fruit_crop_orthologs.py` → `data/fruit_crops/`; UI `/hooks` (grape-trained v13 — not species-retrained).
- **Fruit cross-species high-score comparison** via `scripts/analyze_fruit_cross_species.py` → `data/fruit_crops/cross_species/` (conservative ≥5-species Tier A∪B clusters + species-specific peaks; OrthoDB / grape-join / gene-proxy); UI `/hooks#fruit-cross-species`.
- **Orthogonal LLPS DB crosswalk** via `scripts/fetch_llps_external_ids.py` + `scripts/enrich_llps_crosswalk.py` → `data/llps_crosswalk/` (PhaSePro / LLPSDB / DrLLPS / CD-CODE; unique-to-us filter `evidence=unique`; UI `/crosswalk`).
- **Conserved cross-crop Tier A∪B** (grape↔soy↔Ath OrthoDB) via `scripts/analyze_cross_crop_conserved.py` → `data/cross_crop/`.
- Side-by-side **PSAP / LLPhyScore / PLAAC-approx** via `scripts/enrich_comparator_scores.py`. Full-proteome native PSAP / LLPhyScore via `scripts/score_psap_proteome.py` and `scripts/score_llphyscore_proteome.py` (external envs; batched + checkpointed).
- **Pfam / InterPro domains** + **IDR segments** (AF pLDDT&lt;50 ∪ UniProt MobiDB-lite disordered) via `scripts/enrich_domains_idr.py`; protein-page track + Ensembl Plants link-outs.
- **Bulk download** UI `/download` and `GET /api/download?tier=A&format=tsv|fasta` (optional `family` / `keyword` / `q`).
- **STRING** Vitis PPI neighborhoods + Tier A∪B enrichment via `scripts/enrich_string_neighborhood.py` (UI `/network`); combined scores are **not** LLPS evidence.
- **Stress condensatome modules** via `scripts/enrich_stress_condensatome_modules.py` (UI `/modules`): DE × Tier A∪B enrichment for E-MTAB-4324 / E-GEOD-67932 + response-pattern modules + deposited gene sets (`data/stress_modules/`).
- **Atlas analyses 2/6/7 (v1.5.0)** via `scripts/analyze_stress_network_isoforms.py`: Botrytis-down GO/keyword deep dive + heatmap matrix; STRING/stress seed hubs TSV; VIT multi-accession score deltas (`data/isoforms/` — nearly 1:1 in UP000009183).
- **Atlas analyses family/comparator/structure (v1.5.1)** (UI `/analyses`): Pfam/InterPro Tier A∪B enrichment (`scripts/enrich_family_enrichment.py`); comparator disagreement (`scripts/enrich_comparator_disagreement.py`); structure–score joint (`scripts/enrich_structure_score_joint.py`). Family/structure threshold tests are FDR-controlled; disagreement flags are exploratory.
- **RNA-binding crosswalk (v1.9.0)** via `scripts/enrich_rbp_crosswalk.py` (UI `/rbp`): RBP2GO Ath OrthoDB joins + ATtRACT Ath motifs (sparse; no Vitis) + UniProt GO/KW + RG/RGG grammar. Ortholog/motif flags are triage annotations, not grape CLIP evidence.
- v13 adds **hard disordered non-LLPS negatives**; companion-composite AUC falls vs folded-heavy v12 by design.
- Gold-standard PlantLLPS AUCs from the production model are **partially in-sample** (IDs overlap the training corpus). Prefer the **leakage-controlled fair suite** (`scripts/eval_fair_benchmark.py` → `data/benchmarks/fair_comparison_metrics.json`, leaderboard [`docs/BENCHMARK_LEADERBOARD.md`](./docs/BENCHMARK_LEADERBOARD.md)): leave-gold-out PlantLLPS hybrid ≈ **0.959**, family-hold-out mean ≈ **0.833**, vs OOB LLPhyScore ≈ **0.769** / PSAP ≈ **0.718** on gold. ProX/PU-style is a feature-family logistic baseline (not the published Plant PS proteins model).
- **PSAP** runs in a pinned Python 3.9 env (`scripts/run_psap_gold.py`); gold AUC ≈ **0.718** (companion cited 0.692). System Python 3.12 cannot import `sklearn-json`.
- Companion-cited MAP-screen ρ = **0.368** used the unpublished Supplemental Data S1 PSAP table. This repo deposits an **annotation reconstruction** (`data/map_screen/`) and recomputes ρ vs **native proteome PSAP** (annotation set ≈ 0.20) — not a byte-for-byte replication of 0.368.

## Install

### Python CLI (stdlib; Python ≥3.10)

```bash
git clone https://github.com/pingyguan-glitch/grape-condensatome.git
cd grape-condensatome
python3 -m plantllps --version
```

Optional extras / external tools: [`requirements.txt`](./requirements.txt), [`docs/EXTERNAL_TOOLS.md`](./docs/EXTERNAL_TOOLS.md).

### Node web UI (Node ≥22)

```bash
npm install
npm start
# → http://127.0.0.1:8765/  (HOST/PORT env; default 0.0.0.0:8765)
# Routes: /compare, /gold, /map, /hooks (Ath/rice + soybean/tomato/fruit crops), /network (STRING),
#         /expression, /modules (stress DE enrichment), /analyses (atlas analyses),
#         /rbp (RNA-binding crosswalk), /download, /protein/:id
# API:    GET /api/download?tier=A&evidence=llps|rbp|rg&format=tsv|fasta
```

## RNA-binding crosswalk (NAR nucleic-acid layer)

```bash
python3 scripts/import_rbp_sources.py --rbp2go-csv RBP_AT.csv --fetch-attract
python3 scripts/enrich_rbp_crosswalk.py
# deposits: data/rbp/; UI /rbp ; methods data/rbp/METHOD.md
```

## Atlas analyses (NAR+)

```bash
python3 scripts/enrich_family_enrichment.py
python3 scripts/enrich_comparator_disagreement.py
python3 scripts/enrich_structure_score_joint.py
# deposits: data/family_enrichment/, data/comparator_disagreement/, data/structure_score/
# UI: /analyses
```

## Stress condensatome modules (NAR+)

```bash
python3 scripts/enrich_stress_rnaseq.py                         # prerequisite
python3 scripts/enrich_stress_condensatome_modules.py
python3 scripts/analyze_stress_network_isoforms.py              # v1.5 analyses 2/6/7
# deposits: data/stress_modules/{de_enrichment,coexpression_modules,gene_sets.gmt.gz,deepdive/…}
#           data/string/stress_seed_hubs.tsv · data/isoforms/
# UI: /modules
```

## Second-crop condensatome (soybean)

```bash
python3 scripts/score_second_crop_proteome.py --crop soybean
python3 scripts/enrich_second_crop_orthologs.py --crop soybean
python3 scripts/score_second_crop_proteome.py --crop tomato
python3 scripts/enrich_second_crop_orthologs.py --crop tomato
python3 scripts/score_second_crop_proteome.py --crop arabidopsis --out-dir data/cross_crop
python3 scripts/analyze_cross_crop_conserved.py
# deposits: data/second_crop/{soybean,tomato}_* ; data/cross_crop/conserved_*
# optional: --crop tomato (UP000004994)
```

## Fruit-crop condensatomes (multi-species)

```bash
python3 scripts/score_fruit_crop_proteomes.py --list
python3 scripts/score_fruit_crop_proteomes.py --all
python3 scripts/enrich_fruit_crop_orthologs.py --all   # fruit→grape OrthoDB; crop DB only
python3 scripts/analyze_fruit_cross_species.py         # conservative vs species-specific layer
# deposits: data/fruit_crops/{slug}_condensatome.* + index.json
#          data/fruit_crops/cross_species/{conservative,species_specific,...}
# UI: /hooks fruit-crop panel + /hooks#fruit-cross-species
```

## Train

```bash
python3 -m plantllps train \
  --train-csv data/train_formal_v13.csv \
  --model-out data/model_formal.json

python3 scripts/expand_v13_hard_negatives.py   # rebuild v13 from v12 + hard negs

# Figure 3 feature ablation (composition proxy; deposits under data/benchmarks/)
python3 scripts/eval_feature_ablation.py

# Optional experimental disorder retrain (does NOT overwrite model_formal.json)
export METAPREDICT_PYTHON="$(micromamba run -n meta311 which python)"  # if needed
python3 scripts/retrain_disorder_features.py --metapredict-python "$METAPREDICT_PYTHON"
```

## Score a proteome

```bash
python3 -m plantllps score-proteome \
  --fasta data/UP000009183.fasta.gz \
  --model data/model_formal.json \
  -o data/grape_condensatome.proteome.tsv \
  --db data/grape_condensatome.sqlite

python3 scripts/enrich_uniprot_annotations.py --db data/grape_condensatome.sqlite
python3 scripts/enrich_expression_atlas.py --db data/grape_condensatome.sqlite
python3 scripts/enrich_stress_rnaseq.py --db data/grape_condensatome.sqlite
python3 scripts/enrich_stress_condensatome_modules.py   # DE × Tier A∪B + modules
python3 scripts/enrich_alphafold_plddt.py --segments-only --segments-tier ALL --workers 24
python3 scripts/enrich_string_neighborhood.py   # Vitis STRING PPI + enrichment
python3 scripts/deposit_explainability.py       # SHAP/grammar attributions (gold ∪ Tier A)
python3 scripts/analyze_stress_network_isoforms.py  # deep dive + hubs + isoform deltas
```

Optional metapredict (subset recommended on CPU):

```bash
export METAPREDICT_PYTHON=/path/to/meta311/bin/python
python3 scripts/score_metapredict_subset.py
# or: python3 -m plantllps score proteins.fasta -o out.tsv --use-metapredict
```

## Serve

```bash
python3 -m plantllps serve --host 0.0.0.0 --port 8765
# Or Node UI (publishable Express app)
npm start
```

Public hosting: **https://grape-condensatome.apps.grok.bot** (see [`docs/HOSTING.md`](./docs/HOSTING.md)). Local `npm start` still works. **GitHub Pages will not work** (needs a SQLite server). Deploy configs: `render.yaml`, `railway.toml`, `fly.toml`, `Dockerfile`.

Code/weights/example deposit for Zenodo (tip **v1.11.0**): [`docs/ZENODO.md`](./docs/ZENODO.md) (`python3 scripts/export_zenodo_bundle.py [--parquet]`). Archives land under gitignored `dist/` — upload to Zenodo (or `scripts/zenodo_publish.py` with `ZENODO_ACCESS_TOKEN`); paste the real DOI into `CITATION.cff` after publish — do not invent one.

## Benchmarks vs PSAP / LLPhyScore / PLAAC / ProX-style

```bash
# Fair leakage-controlled suite (preferred for NAR+ claims)
export PSAP_PYTHON=/path/to/psap39/bin/python
export LLPHYSCORE_HOME=/path/to/LLPhyScore/standalone_package/LLPhyScore
export LLPHYSCORE_PYTHON=/path/to/llphy39/bin/python
python3 scripts/eval_fair_benchmark.py
# → data/benchmarks/fair_comparison_metrics.json + docs/BENCHMARK_LEADERBOARD.md

# Legacy in-sample gold QA + native OOB baselines
python3 scripts/compare_gold_standard.py
python3 scripts/run_psap_gold.py
# → data/benchmarks/gold_standard_comparison.json

# Full proteome (resumable)
python3 scripts/score_psap_proteome.py --db-ids-only --batch-size 500
python3 scripts/score_llphyscore_proteome.py --db-ids-only --batch-size 100
python3 scripts/enrich_comparator_scores.py
```

| Predictor | Gold AUC (13+15) | Notes |
|-----------|----------------:|-------|
| PlantLLPS hybrid (leave-gold-out) | ≈ 0.959 | **Fair** — `eval_fair_benchmark.py` |
| PlantLLPS hybrid (production) | ≈ 0.990 | **In-sample** — do not claim as transfer |
| LLPhyScore | ≈ 0.769 | OOB; companion cited 0.774 |
| PSAP (native) | ≈ 0.718 | OOB in py3.9 env; companion cited 0.692 |
| ProX/PU-style (leave-gold-out) | ≈ 0.467 | Feature-family logistic; not published PU |
| PLAAC-approx | ≈ 0.497 | Stdlib approx; not Java PLAAC |

## Expression Atlas

**Baseline TPM:** Experiment **E-MTAB-8758** (Pinot Noir berry) joined on `gene_name` = Atlas `GeneID` (`VIT_*`): **25 810 / 25 972** (99.4%).

**Stress RNA-seq differentials (v1.1):** **E-MTAB-4324** (downy mildew / Trichoderma) and **E-GEOD-67932** (Botrytis noble rot) — both **25 810 / 25 972** (99.4%). UI `/expression` + protein-page contrasts (log2FC / p; sig when p&lt;0.05 and |log2FC|≥1). Abiotic drought/salt Atlas grape entries remain microarray and are not claimed as RNA-seq.


## Version tag notes

Tag **`v0.3.0`** = hard negatives + LLPhyScore + UniProt joins + packaging.  
**`v0.3.1`** = metapredict wiring, native PSAP gold AUC, Expression Atlas joins, hosting docs.  
**`v0.4.0`** = AlphaFold pLDDT joins, batched proteome metapredict, Ath/rice OrthoDB + cross-crop hooks, side-by-side PSAP/LLPhyScore/PLAAC-approx UI.  
**`v0.5.0`** = full-proteome **native** PSAP + LLPhyScore (batched/checkpointed) joined into SQLite + `/compare`.  
**`v0.6.0`** = Pfam/InterPro domains + IDR tracks, bulk download API, Ensembl Plants link-outs.  
**`v1.1.1`** = feature ablation + experimental disorder retrains.  
**`v1.2.0`** = stress DE × condensatome modules (parallel PR).
**`v1.3.0`** = fair NAR+ benchmark suite + leaderboard (leave-gold-out / family / expanded species).  
**`v1.4.0`** = explainability SHAP/grammar attributions on protein pages (gold ∪ Tier A).  
**`v1.5.0`** = stress GO/keyword deep dive + STRING/stress hubs + VIT multi-accession score deltas.
**`v1.5.1`** = family enrichment + comparator disagreement + structure–score (UI `/analyses`).
**`v1.5.2`** = conserved grape↔soy↔Ath Tier A∪B + tomato third-crop condensatome.
**`v1.5.3`** = Supplemental Figures S1–S5 (PNG/PDF + generator).
**`v1.6.0`** = multi-species fruit-crop condensatome pipeline (`data/fruit_crops/` + `/hooks`).  
**`v1.7.0`** = fruit cross-species high-score comparison (`data/fruit_crops/cross_species/` + `/hooks#fruit-cross-species`).

## Still open (not claimed done)

- Exact RBGD4 AGI  
- Zenodo DOI (bundle ready for tip **v1.11.0** under `dist/`; mint on Zenodo then paste DOI — `docs/ZENODO.md`)  
- Host hardening follow-up (publish review: resource exhaustion / data exposure — noted in `docs/HOSTING.md`)  
- Broader splice-isoform annotation (UP000009183 is nearly 1:1 VIT↔accession; v1.5 deposits the rare multi-accession deltas); coexpression-only track  



