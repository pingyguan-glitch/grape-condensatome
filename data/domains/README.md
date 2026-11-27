# Domain / IDR enrichment caches

- `uniprot_up000009183_domains.tsv.gz` — UniProt stream (`ft_domain`, `ft_region`, Pfam/InterPro xrefs) for proteome UP000009183.
- `interpro_pfam_backfill.checkpoint.json` — resumable InterPro API Pfam coordinates for proteins with Pfam xref but empty UniProt Domain FT (gitignored).

Rebuild / join into SQLite:

```bash
python3 scripts/enrich_domains_idr.py
```

IDR intervals use AlphaFold pLDDT&lt;50 segments (when deposited) plus UniProt REGION “Disordered” (MobiDB-lite).
Per-residue metapredict / pLDDT trajectories live under `data/tracks/` (see that README);
entry name cache via `scripts/deposit_interpro_entries.py`.
