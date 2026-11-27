# RNA-binding annotation crosswalk (NAR layer A3)

PlantLLPS / GrapeCondensatome nucleic-acid scope layer: cross LLPS-scored grape
proteins against public RBP resources and RG/RGG sequence grammar.

## Sources

| Resource | Role | Plant / grape coverage |
|---|---|---|
| **RBP2GO 2.0** Ath CSV | Primary experimental+GO RBP catalog | Ath only (2 352 RBPs); **no Vitis** |
| **ATtRACT** | RBP-associated RNA motifs | Animal-heavy; Ath genes ≈ 3 (sparse); **no Vitis** |
| UniProt UP000009183 | Native GO:0003723 + KW RNA-binding | Direct grape annotations |
| Grape FASTA | RG / RGG motif scan | Whole scored proteome |

Import: `scripts/import_rbp_sources.py` → `data/rbp/cache/`.  
Crosswalk: `scripts/enrich_rbp_crosswalk.py`.

## Join logic

1. Map each grape protein’s `ortholog_ath` to RBP2GO Ath `Uniprot_ID`.
2. Attach ATtRACT Ath motif rows when the Ath ortholog / gene resolves.
3. Flag native UniProt GO:0003723 and Keyword RNA-binding.
4. Scan sequences for overlapping `RGG` and `RG`; call when
   `n_RGG ≥ 2` **or** `n_RG ≥ 8` **or** max consecutive RGG-run ≥ 2.

## Headline counts (this deposit)

- Atlas proteins: **25972**
- Any RBP-related annotation: **1943** (7.48% proteome)
- Via RBP2GO Ath ortholog: **598**
- Via ATtRACT Ath motifs: **0**
- RG/RGG motif call: **395**
- Tier A∪B with any RBP annotation: **348**
- RG motif ∩ RBP2GO: **31**

## Honesty

- OrthoDB groups can be broad; an Ath RBP ortholog is **not** proof the grape
  protein binds RNA.
- ATtRACT motif strings are transferred labels for triage, not grape CLIP peaks.
- RG/RGG thresholds are conventional grammar filters used to prioritize
  RNA-related LLPS candidates (stress-granule / RGG-box biology), not binders.
- No new species condensatome expansion in this layer; no Zenodo touch.
