# Expression Atlas cache

Default experiment: **E-MTAB-8758** (*Vitis vinifera* Pinot Noir berry RNA-seq baseline).

- `E-MTAB-8758-tpms.tsv` — gene × assay-group TPMs (GeneID = `VIT_*`)
- `E-MTAB-8758-configuration.xml` — assay-group labels

Re-download / join into SQLite:

```bash
python3 scripts/enrich_expression_atlas.py --db data/grape_condensatome.sqlite
```

Source: https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-8758/
