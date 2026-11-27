# Second-/third-crop condensatome (soybean + tomato)

Additive crop atlases on top of the grapevine condensatome. Scores use
grape-trained PlantLLPS **v13** with composition-proxy disorder; tiers are
recalibrated per proteome. **Not** crop-retrained models.

| File | Role |
|------|------|
| `soybean_condensatome.sqlite` | *Glycine max* UP000008827 slim browse DB |
| `soybean_condensatome.summary.json` | Soy tier counts / honesty |
| `soybean_grape_ortholog_join.json` | OrthoDB join coverage vs grape |
| `tomato_condensatome.sqlite` | *Solanum lycopersicum* UP000004994 slim browse DB |
| `tomato_condensatome.summary.json` | Tomato tier counts / honesty |
| `tomato_grape_ortholog_join.json` | OrthoDB join coverage vs grape |
| `*_condensatome.proteome.tsv.gz` | Full score tables (gitignored; rebuild / Zenodo) |

Grape atlas stores thin `ortholog_gma*` / `ortholog_sly*` columns only (keeps
`grape_condensatome.sqlite` under GitHub’s ~100 MiB limit).

## Reproduce

```bash
python3 scripts/score_second_crop_proteome.py --crop soybean
python3 scripts/enrich_second_crop_orthologs.py --crop soybean

python3 scripts/score_second_crop_proteome.py --crop tomato
python3 scripts/enrich_second_crop_orthologs.py --crop tomato
```

If UniProt proteome download fails after retries, document the blocker and optionally:

```bash
python3 scripts/write_tomato_schema_stub.py --blocker "…"
```

FASTA caches (`UP00000*.fasta.gz`) are gitignored and re-fetched by the score script.

## Honesty

- OrthoDB eukaryote-level groups can be broad; `data/cross_crop_hooks.json` seeds
  MAP65-1 soybean / tomato IDs when grape F6HET4 lacks an OrthoDB xref.
- Conserved grape↔soy↔Ath Tier A/B tables live under `data/cross_crop/`.
