# PlantLLPS + GrapeCondensatome

Plant-aware LLPS propensity scoring (**PlantLLPS**) and a browsable *Vitis vinifera* condensatome atlas (**GrapeCondensatome**) for UniProt proteome **UP000009183**.

**Version:** 0.7.0 · **License:** MIT · **Cite:** see [`CITATION.cff`](./CITATION.cff)

Author: Pingyin Guan (College of Horticulture, China Agricultural University; ORCID [0000-0003-4890-3412](https://orcid.org/0000-0003-4890-3412)).

## What ships

| Piece | Path | Notes |
|-------|------|-------|
| Python package / CLI | `plantllps/` | train, score, score-proteome, serve; optional `--use-metapredict` |
| Training history | `data/train_formal_v*.csv`, `data/training/` | v1–v13 snapshots + metrics |
| Model weights | `data/model_formal.json` | **v13** hybrid logistic + companion |
| Proteome SQLite | `data/grape_condensatome.sqlite` | 25 972 scored rows + UniProt + Expression Atlas |
| Proteome TSV | `data/grape_condensatome.proteome.tsv` | flat table |
| Node web UI | `server.js` + `lib/` + `public/` | Express + `node:sqlite` (Node ≥22) |
| Benchmarks | `data/benchmarks/` | gold vs LLPhyScore + **native PSAP** + metapredict subset |
| External tool docs | [`docs/EXTERNAL_TOOLS.md`](./docs/EXTERNAL_TOOLS.md) | metapredict / PSAP env pins |
| Hosting notes | [`docs/HOSTING.md`](./docs/HOSTING.md) | publish-app + Railway/Render/Fly; not GitHub Pages |

## Honesty (read before citing metrics)

- Default disorder is a **composition proxy**. **metapredict** can be joined proteome-wide via `scripts/score_metapredict_proteome.py` (batched/checkpointed; v13 logistic weights still expect the proxy unless you retrain).
- **AlphaFold pLDDT** summaries + **scored-proteome low-pLDDT segments** via `scripts/enrich_alphafold_plddt.py` (network-bound; resumable cache; browser-like User-Agent).
- **Ath/rice OrthoDB orthologs** + curated `data/cross_crop_hooks.json` via `scripts/enrich_orthologs.py`.
- Side-by-side **PSAP / LLPhyScore / PLAAC-approx** via `scripts/enrich_comparator_scores.py`. Full-proteome native PSAP / LLPhyScore via `scripts/score_psap_proteome.py` and `scripts/score_llphyscore_proteome.py` (external envs; batched + checkpointed).
- **Pfam / InterPro domains** + **IDR segments** (AF pLDDT&lt;50 ∪ UniProt MobiDB-lite disordered) via `scripts/enrich_domains_idr.py`; protein-page track + Ensembl Plants link-outs.
- **Bulk download** UI `/download` and `GET /api/download?tier=A&format=tsv|fasta` (optional `family` / `keyword` / `q`).
- v13 adds **hard disordered non-LLPS negatives**; companion-composite AUC falls vs folded-heavy v12 by design.
- Gold-standard PlantLLPS AUCs are **partially in-sample** (IDs overlap the training corpus). LLPhyScore / PSAP are out-of-the-box.
- **PSAP** runs in a pinned Python 3.9 env (`scripts/run_psap_gold.py`); gold AUC ≈ **0.718** (companion cited 0.692). System Python 3.12 cannot import `sklearn-json`.
- The companion **208-protein MAP screen table** (with PSAP columns) is not in this repo; MAP ρ cannot be recomputed until that table is deposited (companion ρ = 0.368).

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
# Routes: /compare, /hooks, /download, /protein/:id
# API:    GET /api/download?tier=A&format=tsv|fasta
```

## Train

```bash
python3 -m plantllps train \
  --train-csv data/train_formal_v13.csv \
  --model-out data/model_formal.json

python3 scripts/expand_v13_hard_negatives.py   # rebuild v13 from v12 + hard negs
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
python3 scripts/enrich_alphafold_plddt.py --segments-only --segments-tier ALL --workers 24
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

Public hosting: see [`docs/HOSTING.md`](./docs/HOSTING.md). **GitHub Pages will not work** (needs a SQLite server).

## Benchmarks vs PSAP / LLPhyScore

```bash
export LLPHYSCORE_HOME=/path/to/LLPhyScore/standalone_package/LLPhyScore
export LLPHYSCORE_PYTHON=/path/to/llphy39/bin/python
python3 scripts/compare_gold_standard.py

export PSAP_PYTHON=/path/to/psap39/bin/python
python3 scripts/run_psap_gold.py
# → data/benchmarks/gold_standard_comparison.json

# Full proteome (resumable)
python3 scripts/score_psap_proteome.py --db-ids-only --batch-size 500
python3 scripts/score_llphyscore_proteome.py --db-ids-only --batch-size 100
python3 scripts/enrich_comparator_scores.py
```

| Predictor | Gold AUC (13+15) | Notes |
|-----------|----------------:|-------|
| PlantLLPS hybrid | ≈ 0.990 | **In-sample** — do not claim as transfer |
| LLPhyScore | ≈ 0.769 | Runnable; companion cited 0.774 |
| PSAP (native) | ≈ 0.718 | Runnable in py3.9 env; companion cited 0.692 |

## Expression Atlas

**Baseline TPM:** Experiment **E-MTAB-8758** (Pinot Noir berry) joined on `gene_name` = Atlas `GeneID` (`VIT_*`): **25 810 / 25 972** (99.4%).

**Stress RNA-seq differentials (v0.7):** **E-MTAB-4324** (downy mildew / Trichoderma) and **E-GEOD-67932** (Botrytis noble rot) — both **25 810 / 25 972** (99.4%). UI `/expression` + protein-page contrasts (log2FC / p; sig when p&lt;0.05 and |log2FC|≥1). Abiotic drought/salt Atlas grape entries remain microarray and are not claimed as RNA-seq.

## Version tag notes

Tag **`v0.3.0`** = hard negatives + LLPhyScore + UniProt joins + packaging.  
**`v0.3.1`** = metapredict wiring, native PSAP gold AUC, Expression Atlas joins, hosting docs.  
**`v0.4.0`** = AlphaFold pLDDT joins, batched proteome metapredict, Ath/rice OrthoDB + cross-crop hooks, side-by-side PSAP/LLPhyScore/PLAAC-approx UI.  
**`v0.5.0`** = full-proteome **native** PSAP + LLPhyScore (batched/checkpointed) joined into SQLite + `/compare`.  
**`v0.6.0`** = Pfam/InterPro domains + IDR tracks, bulk download API, Ensembl Plants link-outs.  
**`v0.7.0`** = scored-proteome AF low-pLDDT segments + multi-stress RNA-seq differentials.

## Still open (not claimed done)

- MAP-screen table deposit (ρ vs PSAP)  
- Exact RBGD4 AGI; family-held-out PlantLLPS transfer AUC  
- Confirmed public URL after publish-app / Railway / Render / Fly deploy  
- Splice isoforms; STRING / coexpression layers; user score queue  

