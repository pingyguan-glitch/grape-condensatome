# PlantLLPS + GrapeCondensatome

Plant-aware LLPS propensity scoring (**PlantLLPS**) and a browsable *Vitis vinifera* condensatome atlas (**GrapeCondensatome**) for UniProt proteome **UP000009183**.

**Version:** 0.3.0 · **License:** MIT · **Cite:** see [`CITATION.cff`](./CITATION.cff)

Author: Pingyin Guan (College of Horticulture, China Agricultural University; ORCID [0000-0003-4890-3412](https://orcid.org/0000-0003-4890-3412)).

## What ships

| Piece | Path | Notes |
|-------|------|-------|
| Python package / CLI | `plantllps/` | train, score, score-proteome, serve |
| Training history | `data/train_formal_v*.csv`, `data/training/` | v1–v13 snapshots + metrics |
| Model weights | `data/model_formal.json` | **v13** hybrid logistic + companion |
| Proteome SQLite | `data/grape_condensatome.sqlite` | 25 972 scored rows + UniProt joins |
| Proteome TSV | `data/grape_condensatome.proteome.tsv` | flat table |
| Node web UI | `server.js` + `lib/` + `public/` | Express + `node:sqlite` (Node ≥22) |
| Benchmarks | `data/benchmarks/` | gold-standard vs LLPhyScore |
| Hard-neg curation | `data/hard_neg_v13/`, `scripts/expand_v13_hard_negatives.py` | disordered non-LLPS negatives |

## Honesty (read before citing metrics)

- Disorder is a **composition proxy**, not metapredict / AlphaFold pLDDT.
- v13 adds **hard disordered non-LLPS negatives**; companion-composite AUC falls vs folded-heavy v12 by design.
- Gold-standard PlantLLPS AUCs are **partially in-sample** (IDs overlap the training corpus). LLPhyScore / PSAP are out-of-the-box.
- **PSAP** could not be executed here (PyPI `psap` / sklearn-json incompatible with Python 3.12). Companion-reported PSAP gold AUC = 0.692 and MAP-screen Spearman ρ = 0.368 are cited, not re-run.
- The companion **208-protein MAP screen table** (with PSAP columns) is not in this repo; MAP ρ cannot be recomputed until that table is deposited.

## Install

### Python CLI (stdlib; Python ≥3.10)

```bash
git clone https://github.com/pingyguan-glitch/grape-condensatome.git
cd grape-condensatome
# optional: pip install -r requirements.txt  # only needed for optional extras
python3 -m plantllps --version
```

### Node web UI (Node ≥22)

```bash
npm install
npm start
# → http://127.0.0.1:8765/  (HOST/PORT env; default 0.0.0.0:8765)
```

## Train

```bash
# Canonical v13 set (already trained weights in data/model_formal.json)
python3 -m plantllps train \
  --train-csv data/train_formal_v13.csv \
  --model-out data/model_formal.json

# Rebuild v13 from v12 + hard negatives (needs data/hard_neg_v13/seqs/)
python3 scripts/expand_v13_hard_negatives.py
```

Historical rounds v3–v12: `python3 scripts/expand_rounds_v3_v12.py` (needs cached FASTA under `/tmp/plantllps_*/seqs`).

## Score a proteome

```bash
# Uses bundled UP000009183.fasta.gz when present
python3 -m plantllps score-proteome \
  --fasta data/UP000009183.fasta.gz \
  --model data/model_formal.json \
  -o data/grape_condensatome.proteome.tsv \
  --db data/grape_condensatome.sqlite

# Join UniProt gene names / GO / keywords into SQLite
python3 scripts/enrich_uniprot_annotations.py --db data/grape_condensatome.sqlite
```

Single FASTA:

```bash
python3 -m plantllps score path/to/proteins.fasta -o scores.tsv --model data/model_formal.json
```

## Serve

```bash
# Python stdlib UI
python3 -m plantllps serve --host 0.0.0.0 --port 8765

# Or Node UI (publishable Express app)
npm start
```

## Benchmarks vs PSAP / LLPhyScore

```bash
# Requires LLPhyScore standalone unpacked; set LLPHYSCORE_HOME if not at the default path
export LLPHYSCORE_HOME=/path/to/LLPhyScore/standalone_package/LLPhyScore
python3 scripts/compare_gold_standard.py
# → data/benchmarks/gold_standard_comparison.json
```

On this release machine: LLPhyScore gold AUC ≈ **0.769**; PlantLLPS hybrid AUC ≈ **0.990** (**in-sample** — do not claim as independent transfer). PSAP not runnable — see JSON `psap_status`.

## Re-fetch base data

```bash
curl -L -o data/UP000009183.fasta.gz \
  'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28proteome%3AUP000009183%29'
```

## Version tag notes

Tag **`v0.3.0`** marks the submission-packaging cut: hard negatives, LLPhyScore gold compare, UniProt annotation joins, dual UI, MIT + CITATION. See [`RELEASE_NOTES.md`](./RELEASE_NOTES.md).

```bash
git tag -a v0.3.0 -m "PlantLLPS/GrapeCondensatome v0.3.0"
git push origin v0.3.0
```

## Still open (not claimed done)

- metapredict / AlphaFold pLDDT features  
- Native PSAP re-run (older sklearn env) and MAP-screen table deposit  
- Expression Atlas joins; public hosted deployment beyond local/`publish-app`  
- Exact RBGD4 AGI resolution  
