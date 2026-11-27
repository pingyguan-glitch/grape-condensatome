# Second-crop condensatome (soybean)

Additive layer (**v0.9.1**) on top of the grapevine atlas. Does **not** stack on open
Wave 2 gold (**v0.8.0** / PR #7) or STRING (**v0.9.0** / PR #8).

| File | Role |
|------|------|
| `soybean_condensatome.proteome.tsv.gz` | Full PlantLLPS scores (composition-proxy disorder) |
| `soybean_condensatome.sqlite` | Slim browse DB (scores + grape ortholog pointers) |
| `soybean_condensatome.summary.json` | Tier counts / honesty / wall-time |
| `soybean_grape_ortholog_join.json` | OrthoDB join coverage vs grape atlas |

**Proteome:** UniProt **UP000008827** (*Glycine max*), NCBI taxid **3847**.  
**Model:** grape-trained PlantLLPS **v13** (`data/model_formal.json`); tiers recalibrated on this proteome. **Not** a soybean-retrained model.

## Reproduce

```bash
python3 scripts/score_second_crop_proteome.py --crop soybean
python3 scripts/enrich_second_crop_orthologs.py --crop soybean
# optional tomato path (same scripts):
# python3 scripts/score_second_crop_proteome.py --crop tomato
# python3 scripts/enrich_second_crop_orthologs.py --crop tomato
```

FASTA downloads land here as `UP000008827.fasta.gz` (gitignored; re-fetched by the score script).

## Honesty

- Full proteome scored in this deposit (see `summary.json` for `n_scored`).
- OrthoDB eukaryote-level groups can be broad; curated `data/cross_crop_hooks.json` seeds MAP65-1 soybean IDs when grape lacks an OrthoDB xref (e.g. F6HET4).
- Grape SQLite only stores thin `ortholog_gma*` columns + join meta (keeps the atlas under GitHub’s ~100 MB file limit). Full crop scores live in this folder.
