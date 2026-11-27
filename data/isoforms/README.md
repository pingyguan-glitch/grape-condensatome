# Isoform / multi-protein gene score deltas (analysis 7)

PlantLLPS score spreads for UniProt accessions that share the same ``VIT_*`` locus.

## Headline

- Proteins with a parseable VIT ID: **25,894**
- Distinct VIT loci: **25,939**
- Multi-accession VIT loci: **1**
- Deposited rows (delta ≥ 0.05 or all multi if none): **1**
- Max observed delta: **0.069686**

## Deposit

| File | Contents |
|------|----------|
| `multi_accession_vit_score_deltas.tsv` | High-delta (or all multi-accession) VIT genes |
| `summary.json` | Counts + honesty |

## Reproduce

```bash
python3 scripts/analyze_stress_network_isoforms.py
```

## Honesty

Grouping key = VIT_* locus IDs parsed from UniProt gene_name / gene_names. UP000009183 is nearly one UniProt accession per VIT locus, so multi-accession isoform sets are rare. Score deltas are PlantLLPS hybrid score spreads across accessions of the same VIT — not experimentally confirmed splice-isoform LLPS differences. Shared gene symbols that map to distinct VIT IDs (paralogs) are intentionally excluded from this deposit.
