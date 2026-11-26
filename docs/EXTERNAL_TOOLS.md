# External tools: metapredict & PSAP

PlantLLPS core scoring is **stdlib-only**. Disorder prediction (metapredict) and the
PSAP baseline need separate Python environments.

## metapredict

```bash
# Recommended: dedicated env (pulls PyTorch)
micromamba create -n meta311 -c conda-forge python=3.11 pip -y
micromamba run -n meta311 pip install metapredict
export METAPREDICT_PYTHON="$(micromamba run -n meta311 which python)"

# Score gold standards (+ F6HET4) → data/benchmarks/metapredict_subset.tsv
python3 scripts/score_metapredict_subset.py --metapredict-python "$METAPREDICT_PYTHON"

# Or wire into CLI scoring (stores frac_disorder_metapredict in TSV / features_json)
python3 -m plantllps score seed.fasta -o out.tsv --use-metapredict \
  --metapredict-python "$METAPREDICT_PYTHON"

# Full proteome is feasible but slow on CPU (~hours for 26k). Prefer a subset first.
# python3 -m plantllps score-proteome --use-metapredict --metapredict-python "$METAPREDICT_PYTHON"
```

Notes:

- v13 logistic weights were trained on the **composition** `frac_disorder_proxy`.
- `--metapredict-replace-proxy` overwrites that field for companion/tiers only; retrain
  if you want logistic weights to use metapredict.

## PSAP (van Mierlo et al. 2021)

PyPI `psap` 1.0.7 needs **Python ≤3.9**, `scikit-learn~=0.21`, and older NumPy/SciPy.
It does **not** import on system Python 3.12 (`sklearn-json` / `_gb_losses`).

```bash
micromamba create -n psap39 -c conda-forge python=3.9 pip -y
micromamba run -n psap39 pip install 'psap' 'numpy<1.24' 'scipy<1.11'
export PSAP_PYTHON="$(micromamba run -n psap39 which python)"

# Gold-standard AUC vs PlantLLPS / LLPhyScore
python3 scripts/run_psap_gold.py --psap-python "$PSAP_PYTHON"
# → updates data/benchmarks/gold_standard_comparison.json (auc.psap)
# → data/benchmarks/psap_gold_prediction.csv
```

FASTA headers must be UniProt-style (`>sp|ACCESSION|NAME`); the helper script writes them.

**MAP-screen Spearman ρ** still cannot be recomputed until the companion 208-protein
table (with PSAP columns) is deposited.

## LLPhyScore

See README: set `LLPHYSCORE_HOME` to an unpacked standalone tree, then
`python3 scripts/compare_gold_standard.py`.

## Expression Atlas

```bash
python3 scripts/enrich_expression_atlas.py --db data/grape_condensatome.sqlite
# Default experiment E-MTAB-8758 (Pinot Noir berry baseline); GeneID = VIT_* loci
```
