# AlphaFold cache

Per-accession `*.prediction.json` / `*.confidence.json` from alphafold.ebi.ac.uk.
Gitignored (large). Rebuild with:

```bash
python3 scripts/enrich_alphafold_plddt.py --workers 12 --with-segments --segments-tier A
```
