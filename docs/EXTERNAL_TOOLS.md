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

- v13 logistic weights were trained on the **composition** `frac_disorder_proxy`.
- `--metapredict-replace-proxy` overwrites that field for companion/tiers only; retrain
  if you want logistic weights to use metapredict.
- Checkpoint: `data/benchmarks/metapredict_proteome.checkpoint.json` (gitignored).

## AlphaFold pLDDT

```bash
# Mean pLDDT + fraction bins from alphafold.ebi.ac.uk API; optional low-pLDDT
# segments from confidence JSON (User-Agent required; files 403 without it).
python3 scripts/enrich_alphafold_plddt.py --workers 12 --with-segments --segments-tier A
# Full proteome summaries are network-bound (~20–60+ min depending on rate limits).
# Segments default to gold + anchors + Tier A; use --segments-tier ALL for every ID.
```

Columns: `af_mean_plddt`, `af_frac_very_low` (pLDDT&lt;50 fraction; matches companion
“pLDDT-based disorder”), `af_disorder_segments_json`, …

Cache: `data/alphafold/*.json` (gitignored).

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

## Side-by-side comparators (PSAP / LLPhyScore / PLAAC)

```bash
python3 scripts/enrich_comparator_scores.py
# Joins gold PSAP + LLPhyScore from data/benchmarks/gold_standard_scores.tsv
# Computes stdlib PLAAC-approx for the full scored proteome
# → data/benchmarks/comparator_side_by_side.tsv + UI /compare
```

### PSAP (van Mierlo et al. 2021)

PyPI `psap` 1.0.7 needs **Python ≤3.9**, `scikit-learn~=0.21`, and older NumPy/SciPy.
It does **not** import on system Python 3.12.

```bash
micromamba create -n psap39 -c conda-forge python=3.9 pip -y
micromamba run -n psap39 pip install 'psap' 'numpy<1.24' 'scipy<1.11'
export PSAP_PYTHON="$(micromamba run -n psap39 which python)"
python3 scripts/run_psap_gold.py --psap-python "$PSAP_PYTHON"
```

**Full-proteome PSAP** is not shipped (CPU + env pin). Gold columns are joined into
SQLite for the overlapping proteome IDs; UI `/compare` shows the gold panel.

**MAP-screen Spearman ρ** still cannot be recomputed until the companion 208-protein
table (with PSAP columns) is deposited.

### LLPhyScore

Set `LLPHYSCORE_HOME` to an unpacked standalone tree, then
`python3 scripts/compare_gold_standard.py`. Gold percentiles are joined by
`enrich_comparator_scores.py`. Proteome-wide LLPhyScore is not run in CI.

### PLAAC-approx (stdlib)

`plantllps/plaac_approx.py` — windowed Q/N enrichment + NLLR-style score. **Not** the
Java PLAAC HMM; intended as a comparable, always-on column. Companion MAP65 NLLR
values remain near zero / negative under both.

## Expression Atlas

```bash
python3 scripts/enrich_expression_atlas.py --db data/grape_condensatome.sqlite
# Default experiment E-MTAB-8758 (Pinot Noir berry baseline); GeneID = VIT_* loci
```
