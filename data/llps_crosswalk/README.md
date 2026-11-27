# Orthogonal LLPS database crosswalk

Cross-annotation of the PlantLLPS gold-standard panel and GrapeCondensatome Tier A/B/C candidates against **PhaSePro**, **LLPSDB**, **DrLLPS**, and **CD-CODE**.

## Rebuild

```bash
python3 scripts/fetch_llps_external_ids.py          # cache raw dumps → compact ID lists
python3 scripts/enrich_llps_crosswalk.py            # gold + candidates TSV/JSON + SQLite columns
```

Raw dumps stay in gitignored `cache/`. Committed artifacts are UniProt/gene ID tables and join results only.

## Deposits

| Path | Role |
|------|------|
| `ids/phasepro_uniprot.txt` | PhaSePro UniProt accessions |
| `ids/llpsdb_uniprot.txt` | LLPSDB v2 unambiguous PS UniProt IDs |
| `ids/drllps_experimental_uniprot.txt` | DrLLPS experimental UniProt IDs |
| `ids/drllps_vitis_genes.txt` | DrLLPS *Vitis* ortholog-expansion gene IDs |
| `ids/cdcode_uniprot.txt` | CD-CODE UniProt accessions (from Edmond dump; CC BY-NC) |
| `ids/sources.json` | URLs, citations, license notes, freeze counts |
| `gold_crosswalk.tsv` / `.json` | PlantLLPS gold panel × external DBs |
| `candidates_crosswalk.tsv.gz` | Tier A∪B∪C atlas rows × external DBs |
| `summary.json` | Coverage counts (unique vs shared) |

SQLite (`data/grape_condensatome.sqlite`) columns: `ext_llps_dbs`, `ext_llps_unique`, `ext_llps_json`.

## Match rules

1. Direct UniProt ID in PhaSePro / LLPSDB / DrLLPS-experimental / CD-CODE  
2. `ortholog_ath` UniProt in those sets  
3. Grape `gene_name` (`VIT_…`) in DrLLPS Vitis gene list → counted as **DrLLPS**

`ext_llps_unique=1` means no match under any rule (unique to GrapeCondensatome / PlantLLPS relative to these catalogs).

## Honesty

- DrLLPS *Vitis* hits are **ortholog expansions**, not plant droplet assays.
- CD-CODE dump is **CC BY-NC**; we redistribute accession lists for cross-annotation, not the full dump.
- “Unique” ≠ experimentally proven novelty.

See [`docs/LLPS_CROSSWALK.md`](../../docs/LLPS_CROSSWALK.md) for the reviewer-facing methods note.
