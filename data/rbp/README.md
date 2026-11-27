# data/rbp — RNA-binding crosswalk deposits

NAR-priority nucleic-acid annotation layer (v1.9.0).

| File | Contents |
|---|---|
| `rbp2go_arabidopsis.tsv.gz` | Slim RBP2GO 2.0 Ath RBP table |
| `attract_arabidopsis_motifs.tsv.gz` | ATtRACT Ath motif rows |
| `grape_rbp_annotations.tsv.gz` | Per-grape-protein crosswalk |
| `annotations_by_protein.json.gz` | UI sidecar |
| `ui_tables.json` | Browse slices (Tier A∪B + RG hits) |
| `summary.json` | Hit rates / species coverage |
| `METHOD.md` | Methods blurb |
| `cache/` | Raw RBP_AT.csv + ATtRACT.zip (gitignored) |

Rebuild::

    python3 scripts/import_rbp_sources.py --rbp2go-csv RBP_AT.csv --fetch-attract
    python3 scripts/enrich_rbp_crosswalk.py
