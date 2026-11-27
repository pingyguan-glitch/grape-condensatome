# External tools & atlas enrichments

PlantLLPS **core** scoring is stdlib-only. Optional tracks (metapredict, AlphaFold,
PSAP, LLPhyScore) and OrthoDB ortholog joins use network downloads and/or separate
environments. All enrichment scripts are **resumable** (disk caches under `data/`).

## metapredict (disorder)

```bash
# Dedicated env recommended (pulls PyTorch); system Python 3.12 also works if
# python3-dev is present so the cython extension builds.
micromamba create -n meta311 -c conda-forge python=3.11 pip -y
micromamba run -n meta311 pip install metapredict
export METAPREDICT_PYTHON="$(micromamba run -n meta311 which python)"

# Gold + F6HET4 subset
python3 scripts/score_metapredict_subset.py --metapredict-python "$METAPREDICT_PYTHON"

# Full scored proteome (batched + checkpointed; ~5–15 min on a 4-core CPU once
# metapredict is warm — document longer wall-times on cold / memory-constrained hosts)
python3 scripts/score_metapredict_proteome.py --db-ids-only --batch-size 200 \
  --metapredict-python "$METAPREDICT_PYTHON"
# → data/benchmarks/metapredict_proteome.tsv + SQLite frac_disorder_metapredict

# CLI path (stores frac_disorder_metapredict in TSV / features_json)
python3 -m plantllps score seed.fasta -o out.tsv --use-metapredict \
  --metapredict-python "$METAPREDICT_PYTHON"
```

Notes:

- v13 **production** logistic weights were trained on the **composition** `frac_disorder_proxy`
  (`data/model_formal.json`).
- `--metapredict-replace-proxy` overwrites that field for companion/tiers only.
- Optional experimental retrain (comparable CV; does **not** overwrite production)::

    export METAPREDICT_PYTHON="$(micromamba run -n meta311 which python)"
    python3 scripts/retrain_disorder_features.py --metapredict-python "$METAPREDICT_PYTHON"
    # → data/experimental/model_formal.metapredict.json (+ af_plddt / metapredict_af_mean)

- Feature-ablation panels (Figure 3)::

    python3 scripts/eval_feature_ablation.py
    # → data/benchmarks/feature_ablation_metrics.json + feature_ablation_panel.svg

- Checkpoint: `data/benchmarks/metapredict_proteome.checkpoint.json` (gitignored).

## AlphaFold pLDDT

```bash
# Mean pLDDT + fraction bins from alphafold.ebi.ac.uk API; low-pLDDT segments from
# confidence JSON (browser-like User-Agent required; file CDN 403 otherwise).
python3 scripts/enrich_alphafold_plddt.py --workers 16 --with-segments --segments-tier ALL
# Resume segments only (skips IDs that already have af_disorder_segments_json):
python3 scripts/enrich_alphafold_plddt.py --segments-only --segments-tier ALL --workers 24
# Full proteome is network-bound (~15–45+ min depending on rate limits / workers).
# v1.1 default segment scope = ALL scored proteins with an AFDB model.
```

Columns: `af_mean_plddt`, `af_frac_very_low` (pLDDT&lt;50 fraction; matches companion
“pLDDT-based disorder”), `af_disorder_segments_json`, …

Cache: `data/alphafold/*.json` (gitignored). Coverage: `data/alphafold/coverage.json` + `meta.alphafold_plddt`.

## Stress RNA-seq (Expression Atlas)

```bash
python3 scripts/enrich_stress_rnaseq.py --db data/grape_condensatome.sqlite
# Default: E-MTAB-4324 (downy mildew / Trichoderma) + E-GEOD-67932 (Botrytis noble rot)
```

Writes `data/expression_atlas/stress_by_protein.json.gz` (UI sidecar) +
`stress_rnaseq_join_report.json` / `meta.stress_rnaseq`. Optional
`--also-sqlite-column` is off by default (GH 100 MiB SQLite limit).

## Orthologs (Arabidopsis / rice) + cross-crop hooks

```bash
python3 scripts/enrich_orthologs.py
# Streams UniProt OrthoDB xrefs for grape UP000009183, Ath UP000006548,
# rice japonica UP000059680; joins on shared OrthoDB group; overlays
# data/cross_crop_hooks.json (MAP65-1, WVD2/WDL, …).
```

Honest limit: OrthoDB eukaryote-level groups can be broad (e.g. large enzyme
families). Curated hooks seed high-value families even when xref coverage is thin
(F6HET4 itself often lacks an OrthoDB xref in UniProt).

## Second-crop condensatome (soybean / tomato)

```bash
# Soybean Glycine max UP000008827 (default; full proteome ~71k sequences, ~3–5 min)
python3 scripts/score_second_crop_proteome.py --crop soybean
python3 scripts/enrich_second_crop_orthologs.py --crop soybean

# Optional tomato Solanum lycopersicum UP000004994
python3 scripts/score_second_crop_proteome.py --crop tomato
python3 scripts/enrich_second_crop_orthologs.py --crop tomato
```

Deposits under `data/second_crop/` (TSV.gz + slim SQLite). Grape atlas gains thin
`ortholog_gma*` (or `ortholog_sly*`) columns only — keeps `grape_condensatome.sqlite`
under GitHub’s ~100 MB limit. UI: `/hooks` soybean panel. Uses grape-trained PlantLLPS
v13 with composition-proxy disorder; tiers recalibrated per crop proteome.

## Fair NAR+ benchmark suite

Leakage-controlled same-set comparison (PlantLLPS vs PSAP / LLPhyScore / PLAAC-approx /
ProX-style). See [`BENCHMARK_LEADERBOARD.md`](./BENCHMARK_LEADERBOARD.md).

```bash
export PSAP_PYTHON="$(micromamba run -n psap39 which python)"
export LLPHYSCORE_HOME=/tmp/LLPhyScore/standalone_package/LLPhyScore
export LLPHYSCORE_PYTHON="$(micromamba run -n llphy39 which python)"  # optional
python3 scripts/eval_fair_benchmark.py
# → data/benchmarks/fair_comparison_*.json|tsv + docs/BENCHMARK_LEADERBOARD.md
# Reuse deposited OOB cache: python3 scripts/eval_fair_benchmark.py --skip-oob-run
```

## Side-by-side comparators (PSAP / LLPhyScore / PLAAC)

```bash
python3 scripts/enrich_comparator_scores.py
# Joins gold + optional proteome native TSVs; PLAAC-approx for full proteome
# → data/benchmarks/comparator_side_by_side.tsv + UI /compare
```

### PSAP (van Mierlo et al. 2021)

PyPI `psap` 1.0.7 needs **Python ≤3.9**, `scikit-learn~=0.21`, and older NumPy/SciPy.
It does **not** import on system Python 3.12.

```bash
micromamba create -n psap39 -c conda-forge python=3.9 pip -y
micromamba run -n psap39 pip install 'psap' 'numpy<1.24' 'scipy<1.11'
export PSAP_PYTHON="$(micromamba run -n psap39 which python)"

# Gold-standard AUC
python3 scripts/run_psap_gold.py --psap-python "$PSAP_PYTHON"

# Full scored proteome (batched + checkpointed; ~5–15 min on a 4-core CPU)
python3 scripts/score_psap_proteome.py --db-ids-only --batch-size 500 \
  --psap-python "$PSAP_PYTHON"
# → data/benchmarks/psap_proteome.tsv + SQLite score_psap
```

Checkpoint: `data/benchmarks/psap_proteome.checkpoint.json` (gitignored).

**MAP-screen Spearman ρ:** companion Supplemental Data S1 (208 IDs + original PSAP
columns) was never recovered. `scripts/deposit_map_screen.py` deposits an
annotation reconstruction under `data/map_screen/` and recomputes ρ vs **native
proteome PSAP** (annotation set companion-vs-PSAP ≈ 0.20). Cite companion ρ =
0.368 for the published transfer gap; cite the deposit for the public atlas.

### LLPhyScore (Cai / Forman-Kay lab)

```bash
git clone https://github.com/julie-forman-kay-lab/LLPhyScore.git /tmp/LLPhyScore
cd /tmp/LLPhyScore/standalone_package && tar -xf LLPhyScore-standalone-20220627.tar.gz
export LLPHYSCORE_HOME=/tmp/LLPhyScore/standalone_package/LLPhyScore

# Pinned deps (standalone README pins older wheels; py3.9 + numpy 1.19 works)
micromamba create -n llphy39 -c conda-forge python=3.9 pip numpy=1.19.5 \
  pandas=1.1.5 scipy=1.5.3 tqdm scikit-learn -y
export LLPHYSCORE_PYTHON="$(micromamba run -n llphy39 which python)"

# Gold-standard AUC
python3 scripts/compare_gold_standard.py --llphyscore-home "$LLPHYSCORE_HOME"

# Full scored proteome (batched + checkpointed; ~2–4 h single-process on a 4-core CPU)
python3 scripts/score_llphyscore_proteome.py --db-ids-only --batch-size 100 \
  --llphyscore-home "$LLPHYSCORE_HOME" --llphyscore-python "$LLPHYSCORE_PYTHON"
# → data/benchmarks/llphyscore_proteome.tsv + SQLite score_llphyscore_percentile

# Optional: fan-out across CPU cores with separate shard checkpoints, then merge
# (each shard uses --fasta of a partitioned ID set + its own --checkpoint/--out)
python3 scripts/merge_llphyscore_shards.py
```

Checkpoint: `data/benchmarks/llphyscore_proteome.checkpoint.json` (gitignored; shard
files `llphyscore_proteome.shard*.checkpoint.json` likewise). Score column is the
**8-feature-sum human+PDB percentile** (same as gold compare).

### PLAAC-approx (stdlib)

`plantllps/plaac_approx.py` — windowed Q/N enrichment + NLLR-style score. **Not** the
Java PLAAC HMM; intended as a comparable, always-on column. Companion MAP65 NLLR
values remain near zero / negative under both.

## Expression Atlas

```bash
python3 scripts/enrich_expression_atlas.py --db data/grape_condensatome.sqlite
# Default experiment E-MTAB-8758 (Pinot Noir berry baseline); GeneID = VIT_* loci
```

## Domains / IDR tracks

```bash
python3 scripts/enrich_domains_idr.py
# UniProt FT Domain + Pfam/InterPro xrefs; InterPro API backfill for Pfam coords
# IDR = AF pLDDT<50 segments (when present) ∪ UniProt REGION Disordered (MobiDB-lite)
# → SQLite domains_json / idr_segments_json; UI protein-page track
```

Cache: `data/domains/uniprot_up000009183_domains.tsv.gz` (committed). InterPro checkpoint
`data/domains/interpro_pfam_backfill.checkpoint.json` is gitignored.

## STRING neighborhood (Vitis PPI)

```bash
python3 scripts/enrich_string_neighborhood.py
# Downloads taxid 29760 (Vitis) aliases / info / links from stringdb-downloads.org (v12.0 files)
# Maps STRING ↔ UniProt (~99.6% of atlas); builds 1-hop neighborhoods at combined_score ≥ 400
# → SQLite string_neighborhood (Tier A seeds only) + TSV.gz sidecars; UI /network + protein page
# Deposits: enrichment_summary.json, neighborhood_stats.tsv.gz, top_neighbors.tsv.gz,
#           seed_edges.tsv.gz, README.md
```

Notes:

- **Combined score** (0–1000) mixes experiments, databases, co-expression, textmining,
  neighborhood, fusion, and co-occurrence — **not** LLPS-specific evidence.
  Conventional cutoffs: 150 low / 400 medium / 700 high / 900 highest.
- Species is **native grape** STRING (29760), not Arabidopsis transfer. A small set of
  atlas IDs remain unmapped (~100).
- Full `*.protein.links*.txt.gz` caches are gitignored (~47 MB); re-download on demand.
- SQLite keeps **Tier A seed** neighborhood rows only (GitHub 100 MB SQLite limit).
  Full-proteome stats / top-N / seed edges are committed as `data/string/*.tsv.gz`;
  the Node UI loads those sidecars at startup.
- Enrichment fold = (neighbor Tier A∪B fraction) / (proteome A∪B rate); hypergeometric
  p-value is descriptive adjacency bias, not condensate membership.
