# Stress condensatome modules (NAR+ omics layer)

Condensatome × stress DE enrichment for **E-MTAB-4324** and **E-GEOD-67932**, plus
coexpression modules among high-tier (PlantLLPS Tier A∪B) proteins.

## Reproduce

```bash
python3 scripts/enrich_stress_rnaseq.py          # prerequisite sidecar
python3 scripts/enrich_stress_condensatome_modules.py
```

## Deposits

| File | Contents |
|------|----------|
| `de_enrichment.tsv` | Per-contrast × direction hypergeometric enrichment vs Tier A∪B |
| `coexpression_modules.tsv` | Response-pattern modules (shared signed DE calls) |
| `seed_neighborhoods.tsv.gz` | Tier A seed top-15 partners (r≥0.7) |
| `module_members.tsv.gz` | Membership for every gene set |
| `gene_sets.gmt.gz` | GMT deposit (gene names; UniProt fallback) |
| `ui_tables.json` | Compact payload for Node `/modules` |
| `summary.json` | Run summary + honesty note |

## Parameters

- Background: stress-matched atlas proteins (**25,810**); Tier A∪B in background **3,857** (rate 0.1494)
- DE call: p < 0.05 and |log2FC| ≥ 1.0
- Pattern modules: shared −1/0/+1 signature across contrasts; size ≥ 5
- Seed neighborhoods: Pearson r on stress log2FC (≥4 shared); top-15 at r≥0.7
- Profiled Tier A∪B: **3,857**; seed partner links r≥0.85: **307,885**; modules: **42**; gene sets: **1198**

## Honesty

DE call rule matches stress join (p<0.05, |log2FC|≥1). Hypergeometric enrichment of Tier A∪B among DE genes is descriptive. Response-pattern modules group Tier A∪B proteins by shared signed significant DE calls across Atlas contrasts; seed neighborhoods use Pearson r on the same short log2FC vectors — not a dedicated multi-condition coexpression atlas. Berry TPM E-MTAB-8758 has only 6 assay groups and is not the primary module signal. Not WGCNA.
