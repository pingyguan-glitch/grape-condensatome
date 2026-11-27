# Curated experimentally anchored LLPS layer

Structured positive / mechanism-excluded entries used as a **resource-quality** gold overlay for GrapeCondensatome (not a new training corpus by itself).

| File | Role |
|------|------|
| `llps_experimental_gold.tsv` | Master table: UniProt ID, species, condensate type, assay/evidence, citation, grape ortholog when known |
| `llps_experimental_gold.json` | Same rows + join summary (written by enrich script) |

## Provenance

Prefer literature already cited in the PlantLLPS / companion MAP manuscripts and deposited train notes:

- Companion mechanism-matched gold (AtMAP65-1, ELF3, FLOE1, RBGD2, SEUSS, FCA, animal controls; CRY2 excluded)
- Modest literature / LLPSDB plant expansions already present in `train_formal_v13.csv`
- Grape ortholog column filled only when curated (`F6HET4`, `A0ACT1LGI9`) or when the atlas OrthoDB join resolves a grape ID at enrich time

## Honesty

- This is a **curated overlay**, not an exhaustive plant LLPS database.
- Several train-round entries are mechanism-associated rather than droplet-assay gold; `assay_evidence` and `citation` fields state that.
- CRY2 remains `excluded_mechanism_mismatch` (folded-domain oligomerization).
- PlantLLPS gold AUCs on the 13+15 reconstructed set remain **partially in-sample**; see `data/benchmarks/family_holdout_metrics.json` for family / transfer metrics.

## Join

```bash
python3 scripts/enrich_curated_llps_gold.py
```

Adds SQLite columns `llps_gold` / `llps_gold_json` on matching grape proteins and stores the curated table under `meta`.
