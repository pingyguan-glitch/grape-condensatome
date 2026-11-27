# Methods note — how GrapeCondensatome / PlantLLPS differs from existing LLPS databases

**Audience:** NAR Database reviewers asking *“how do you differ from LLPSDB / PhaSePro / DrLLPS / CD-CODE?”*  
**Release:** v1.8.0 · deposits under [`data/llps_crosswalk/`](../data/llps_crosswalk/)

## Short answer

Those resources are **community knowledgebases** of (mostly animal/yeast) experimental condensates, with DrLLPS further **propagating orthologs** into plant proteomes. GrapeCondensatome is a **specialty-crop scored atlas**: whole *Vitis vinifera* UP000009183 PlantLLPS tiers plus physiology/structure joins. Orthogonal cross-annotation makes the difference concrete: proteins already in those catalogs are badged, and **Tier A/B/C proteins absent from all four** are filterable as a *unique-to-us* candidate set (hypotheses, not wet-lab claims).

## What we cross

| Panel | Definition | Freeze *n* |
|-------|------------|----------:|
| PlantLLPS gold-standard | `train_formal_v13` ∪ curated experimental gold ∪ reconstructed gold-standard CSV (unique UniProt IDs) | **164** |
| Browsable candidates | GrapeCondensatome Tier **A∪B∪C** | **13 523** |

## External sources (UniProt-mapped / gene ID tables)

| Database | Artifact used | License / redistribution note |
|----------|---------------|-------------------------------|
| **PhaSePro** | `download_full.json` (UniProt keys) | Cite Mészáros et al. 2020 *NAR*; ID list only |
| **LLPSDB** v2 | Unambiguous `Phase_separation` `protein.xls` UniProt column | Free download; ID list only (Li et al. 2020; v2 update) |
| **DrLLPS** | `LLPS.txt` experimental UniProt + `Vitis_vinifera.zip` gene headers | Cite Ning et al. 2020 *NAR*; no FASTA redistributed |
| **CD-CODE** | Edmond dump `protein2cdcode` / protein CSV (doi:10.17617/3.WDM3PG) | **CC BY-NC 4.0** — accession list for annotation only |

Rebuild: `python3 scripts/fetch_llps_external_ids.py` then `python3 scripts/enrich_llps_crosswalk.py`.

## Match rules

1. Direct UniProt accession overlap.  
2. Grape `ortholog_ath` UniProt overlap.  
3. Grape `gene_name` in DrLLPS *Vitis* ortholog-expansion gene list → **DrLLPS** match.

## Coverage (v1.8.0 freeze)

| Set | *n* | Unique to us | Shared with ≥1 DB |
|-----|----:|-------------:|------------------:|
| Gold panel | 164 | 100 | 64 |
| Tier A | 1 142 | 974 | 168 |
| Tier B | 2 728 | 2 422 | 306 |
| Tier C | 9 653 | 8 458 | 1 195 |
| Tier A∪B | 3 870 | **3 396** | 474 |
| Tier A∪B∪C | 13 523 | 11 854 | 1 669 |

Shared grape hits are dominated by **DrLLPS Vitis ortholog expansion** (not PhaSePro/LLPSDB direct plant curation). Tier **B/C unique** rows are the visible “new resource” slice for reviewers: plant-scored candidates not already catalogued under these joins.

## Limits (state explicitly)

- Crosswalk is **ID/ortholog/gene** based — not sequence identity clustering beyond OrthoDB Ath joins already in the atlas.  
- DrLLPS plant rows ≠ experimental plant LLPS.  
- CD-CODE NC license constrains commercial reuse of their dump; our accession table is for transparent overlap counting.  
- Gold panel *n*=164 is the deposited union (train + curated + reconstructed); it is not an exhaustive plant LLPS census.  
- Unique flags are **catalog-relative**, not experimental priority scores.

## UI / API

- Protein pages: external DB badges + unique/shared status.  
- `/crosswalk` summary + unique-candidate table.  
- `/download?evidence=unique` (optional `tier=B`) and `GET /api/crosswalk`.
