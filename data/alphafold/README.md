# AlphaFold cache

Per-accession `*.prediction.json` / `*.confidence.json` from alphafold.ebi.ac.uk.
Gitignored (large). Browser-like User-Agent required (file CDN returns HTTP 403 otherwise).

Rebuild (v0.7 — full scored-proteome low-pLDDT segments):

```bash
# Summaries (if needed) + segments for every AFDB model in the scored proteome
python3 scripts/enrich_alphafold_plddt.py --workers 16 --with-segments --segments-tier ALL

# Or resume segments only (skips IDs that already have af_disorder_segments_json)
python3 scripts/enrich_alphafold_plddt.py --segments-only --segments-tier ALL --workers 24
```

Coverage summary is written to `coverage.json` and SQLite `meta.alphafold_plddt`.
