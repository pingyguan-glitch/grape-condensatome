# Fruit-crop condensatomes

Multi-species fruit-tree / fruit-crop PlantLLPS scores (grape-trained **v13**,
composition-proxy disorder). Tier cutoffs recalibrated per proteome.
**Not** species-retrained models.

| File | Role |
|------|------|
| `species_manifest.yaml` | Pinned UniProt proteome IDs + unavailable list |
| `index.json` | Aggregate n / tier counts after scoring |
| `failures.json` | Download/score failures + no-proteome species |
| `{slug}_condensatome.sqlite` | Slim browse DB (scores + optional grape/Ath ortholog pointers) |
| `{slug}_condensatome.summary.json` | Per-species tier counts / honesty |
| `{slug}_grape_ortholog_join.json` | OrthoDB fruit→grape join coverage |
| `*_condensatome.proteome.tsv.gz` | Full score tables (gitignored; rebuild / Zenodo) |
| `*.fasta.gz` | UniProt FASTA caches (gitignored) |

Linked deposits (not re-scored here): grape `data/grape_condensatome.sqlite`,
soybean/tomato under `data/second_crop/`.

## Reproduce

```bash
python3 scripts/score_fruit_crop_proteomes.py --list
python3 scripts/score_fruit_crop_proteomes.py --all
# optional OrthoDB fruit→grape (+ Ath via grape) joins — crop DB only:
python3 scripts/enrich_fruit_crop_orthologs.py --all
```

UI: `/hooks` → Fruit-crop condensatomes panel. Health: `/api/health` → `second_crop.fruit_crops`.

## Honesty

- Cross-species transfer of grape-trained weights; do not treat ranks as experimental LLPS.
- OrthoDB eukaryote-level groups can be broad; no new columns are added to the grape atlas
  for fruit species (keeps `grape_condensatome.sqlite` ≤100 MiB).
- Species without a UniProt reference/near-complete proteome are listed in the manifest
  `unavailable` block — no invented deposits.
