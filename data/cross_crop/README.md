# Cross-crop conserved condensatome

Grape ↔ soybean ↔ *Arabidopsis* Tier **A∪B** overlap via OrthoDB groups.

| File | Role |
|------|------|
| `conserved_tier_ab_overlap.tsv.gz` | Per grape Tier A/B protein: soy/Ath IDs, tiers, `overlap_class` |
| `conserved_orthodb_clusters.tsv.gz` | OrthoDB groups with ≥1 grape Tier A/B member |
| `conserved_summary.json` | Protein- and cluster-level counts + honesty |
| `arabidopsis_condensatome.sqlite` | Ath PlantLLPS scores (slim; tiers recalibrated on Ath) |
| `arabidopsis_condensatome.summary.json` | Ath scoring summary |

**Ath scoring** uses the same grape-trained PlantLLPS v13 path as second-crop scripts
(`--crop arabidopsis --out-dir data/cross_crop`). Not an Ath-retrained model.

## Reproduce

```bash
# Ath scores (once)
python3 scripts/score_second_crop_proteome.py --crop arabidopsis --out-dir data/cross_crop

# Overlap tables (requires deposited soybean SQLite + grape OrthoDB columns)
python3 scripts/analyze_cross_crop_conserved.py
```

## Headline counts (v1.5.0 deposit)

See `conserved_summary.json`. Typical classes:

- `triple_tier_ab` — grape ∩ soybean ∩ Ath all Tier A/B
- `grape_soy_ab_ath_ortholog` — grape∩soy Tier A/B + Ath OrthoDB link (Ath not Tier A/B)
- `grape_specific` — grape Tier A/B without soy Tier A/B and without Ath Tier A/B

Honest limit: OrthoDB eukaryote-level groups can be broad; curated hooks in
`data/cross_crop_hooks.json` cover MAP65-1 when OrthoDB is missing on F6HET4.
