# Fruit-crop cross-species high-score comparison

Deposit produced by `scripts/analyze_fruit_cross_species.py` from existing
`data/fruit_crops/*_condensatome.sqlite` scores (**no re-scoring**).

## Panel

- **18** scored fruit-crop proteomes under `data/fruit_crops/` (status=`score`
  in `index.json`). Linked grape / soybean / tomato atlases are **not** folded into
  these tables (they remain on `/hooks` separately).
- Unavailable taxa stay in the fruit-crop manifest `unavailable` block — not invented.

## High-score definition

- **High** = PlantLLPS Tier **A∪B** on that species' recalibrated cutoffs
  (grape-trained v13 + composition-proxy disorder; not species-retrained).

## Orthology / grouping (approximation)

Cluster keys, in priority order:

1. **`orthodb:…`** — UniProt OrthoDB group id written by
   `scripts/enrich_fruit_crop_orthologs.py` (eukaryote-level; can be broad).
2. **`grape:…`** — fruit→grape OrthoDB join accession when OrthoDB id is missing
   but a grape pointer exists.
3. **`gene:…`** — normalized UniProt primary gene symbol proxy when OrthoDB xrefs
   are absent (common for several fruit proteomes). Gene symbols are **not** true
   orthology; cross-species collisions are rare in this panel.
4. **`acc:…@slug`** — singleton accession when no joinable key exists.

Species with dense OrthoDB xrefs dominate multi-species clusters
(~12 / 18 species contribute ≥1 OrthoDB-linked Tier A/B
row in this deposit). Proteomes lacking OrthoDB xrefs mostly enter as gene-proxy
or singleton clusters and therefore enrich the **species-specific** table — that is
expected, not a bug.

## Conservative LLPS candidates

- Cluster has high-score members in **≥ 5** panel species.
- Ranked by `n_species_high` ↓, then `median_high_score` ↓, then `mean_high_score` ↓.
- Per-species representative = highest-score member in that species (Tier A preferred
  over B only via score/tier sort already used for “best member”).

## Species-specific peaks

A cluster is a **species-specific peak** when:

- high-score in **≤ 2** panel species, **and**
- either:
  - **panel-contrast peak**: present (any tier) in **≥ 3** species, with
    the non-peak species not Tier A/B; or
  - **orphan high**: high in 1–2 species and present in fewer than
    3 species (typical for gene-proxy / singleton / OrthoDB-sparse taxa).

Peak species = the species (or species list) carrying the high-score member(s).

## Files

| File | Role |
|------|------|
| `conservative_candidates.tsv.gz` | Conservative cluster table |
| `species_specific_peaks.tsv.gz` | Species-specific peak table |
| `cluster_members.tsv.gz` | Per-member long table for clusters in either class |
| `summary.json` | Counts, thresholds, honesty, species list |
| `browse.json` | UI-sized slice (top conserved + per-species peaks) |
| `README.md` | This methods note |

## Honesty

- Ranks are **hypotheses** from cross-species transfer of grape-trained weights.
- OrthoDB eukaryote groups and gene-name proxies are **near-ortholog approximations**,
  not curated 1:1 orthologs.
- Do **not** treat conservation as experimental LLPS evidence.

## Reproduce

```bash
python3 scripts/analyze_fruit_cross_species.py
# optional threshold knobs:
python3 scripts/analyze_fruit_cross_species.py --min-conserved 5 --max-specific 2
```

This deposit: **75544** scored clusters · conservative **1438** ·
species-specific **73336**.
