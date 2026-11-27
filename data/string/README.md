# STRING neighborhood deposit (Vitis)

- **Source:** STRING v12.0 downloads for NCBI taxid **29760** (*Vitis vinifera*)
  from `stringdb-downloads.org`.
- **Threshold:** combined_score ≥ **400** (medium by default).
- **Mapping:** UniProt_AC (+ preferred name / OrderedLocusName fallbacks) → atlas IDs.
  Coverage: **25,867 / 25,972**.
- **SQLite seeds:** Tier **A** → **833** rows in `string_neighborhood`
  (top-N JSON included).
- **Full proteome (TSV.gz):** `neighborhood_stats.tsv.gz` (19,523 rows),
  `top_neighbors.tsv.gz`, `seed_edges.tsv.gz` (44,932 edges).
- **Enrichment:** Tier A∪B fraction among 1-hop neighbors vs proteome background
  (0.1490); hypergeometric p-value. See `enrichment_summary.json`.

## Honesty

Combined scores mix experimental, database, co-expression, and **textmining** channels —
they are **not** LLPS evidence. This pass uses **native grape** STRING, not Arabidopsis
transfer. Large raw `*.protein.links*.txt.gz` caches are gitignored; re-download with
`python3 scripts/enrich_string_neighborhood.py`.

SQLite seed-only storage keeps `grape_condensatome.sqlite` under GitHub's 100 MB limit;
the Node UI loads the TSV.gz sidecars for non-seed protein pages.
