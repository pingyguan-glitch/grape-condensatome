# Stress module deep dive (analysis 2)

GO-ID + UniProt **keyword** enrichment for condensatome × stress DE gene sets, plus a
heatmap-ready log2FC matrix for the FDR-significant Botrytis down set.

## Primary set

- **DE_E-GEOD-67932_g1_g4_down_TierAB** — E-GEOD-67932 stage-3 Botrytis down ∩ Tier A∪B
- Members: **493**
- DE∩TierAB FDR: **0.001** (fold ≈ 1.167796)

## Deposits

| File | Contents |
|------|----------|
| `go_keyword_enrichment.tsv` | GO + keyword hypergeometric tests; `scope=primary\|exploratory` |
| `heatmap_matrix.tsv.gz` | Primary members × contrast log2FC (+ significance flags) |
| `summary.json` | Run summary + honesty |

## Reproduce

```bash
python3 scripts/analyze_stress_network_isoforms.py
```

## Honesty

Primary scope = FDR-significant Botrytis stage-3 down ∩ Tier A∪B only. All other contrast×direction GO/keyword tests are exploratory and must not be cited as confirmed pathway hits. Enrichment is descriptive hypergeometric overlap on UniProt annotations — not g:Profiler/clusterProfiler with ontology propagation, and not KEGG pathway analysis.

GO IDs and UniProt keywords from proteins.go_terms / keywords. No KEGG pathway table is deposited in this atlas; keyword enrichment is the KEGG-style functional-class proxy.
