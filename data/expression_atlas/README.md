# Expression Atlas cache

## Baseline TPM

**E-MTAB-8758** (*Vitis vinifera* Pinot Noir berry RNA-seq baseline).

- `E-MTAB-8758-tpms.tsv` — gene × assay-group TPMs (GeneID = `VIT_*`)
- `E-MTAB-8758-configuration.xml` — assay-group labels

```bash
python3 scripts/enrich_expression_atlas.py --db data/grape_condensatome.sqlite
```

Source: https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-8758/

## Stress RNA-seq differentials (v0.7)

Two public Atlas **RNA-seq** differentials with `VIT_*` GeneIDs (same join key):

| Accession | Stress | Why chosen |
|-----------|--------|------------|
| **E-MTAB-4324** | Downy mildew (*Plasmopara*) ± *Trichoderma* T39 (leaves) | Clear biotic-stress contrasts; Atlas analytics.tsv |
| **E-GEOD-67932** | *Botrytis* noble-rot stages vs control (berries) | Independent biotic stress; berry tissue; analytics.tsv |

```bash
python3 scripts/enrich_stress_rnaseq.py --db data/grape_condensatome.sqlite
```

Join report: `stress_rnaseq_join_report.json` (+ SQLite `meta.stress_rnaseq`).
Per-protein DE blobs: `stress_by_protein.json.gz` (sidecar; keeps main SQLite under GitHub 100 MiB).

**Not joined as RNA-seq:** abiotic drought/salt Atlas grape entries **E-GEOD-31594** / **E-GEOD-31677** are microarray (documented honesty; no invented RNA-seq TPMs).

## Condensatome × stress enrichment (v1.2)

```bash
python3 scripts/enrich_stress_condensatome_modules.py
```

Deposits DE ∩ Tier A∪B enrichment, response-pattern modules, seed neighborhoods, and GMT gene sets under `data/stress_modules/` (UI `/modules`). See that folder’s README for parameters and honesty notes.
