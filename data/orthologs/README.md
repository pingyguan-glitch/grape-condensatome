# Ortholog UniProt OrthoDB streams

Cached `*_orthodb.tsv.gz` for grape / Arabidopsis / rice / soybean / tomato proteomes. Gitignored.

```bash
python3 scripts/enrich_orthologs.py                 # Ath + rice
python3 scripts/enrich_second_crop_orthologs.py --crop soybean
python3 scripts/enrich_second_crop_orthologs.py --crop tomato
python3 scripts/analyze_cross_crop_conserved.py     # grape↔soy↔Ath Tier A/B overlap
```
