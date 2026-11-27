# Per-residue annotation tracks (NAR A2)

Offline deposits for **AlphaFold pLDDT** as an orthogonal disorder rail next to
**metapredict** trajectories. Prefer cache over live API for NAR reproducibility.

| File | Contents |
|------|----------|
| `grape_residue_tracks.sqlite` | Grape UP000009183 — uint8 **metapredict disorder** + **AFDB pLDDT** |
| `grape_residue_tracks.coverage.json` | Grape deposit report |
| `crop_candidates_residue_tracks.sqlite` | High-score fruit/crop **conservative candidates** (+ grape ortholog pointers) |
| `crop_candidates_residue_tracks.ids.json` | Species → accession lists used for the candidate deposit |
| `candidate_sequences.fasta` | UniProt FASTA cache for candidate disorder (gitignored if large) |
| `interpro_entries.json.gz` | Optional InterPro/Pfam name cache (non-blocking for A2) |

Schema (`residue_tracks`): `protein_id`, `length`, `disorder_u8` (0–255 ← [0,1]),
`plddt_u8` (0–100), sources/means, `updated_at`. Helpers: `scripts/tracks_lib.py`.

## Critical path (A2)

```bash
# 1) Grape proteome disorder + pLDDT (run separately to avoid SQLite writer contention)
python3 scripts/deposit_residue_tracks.py --species grape --disorder --batch-size 50
python3 scripts/deposit_residue_tracks.py --species grape --plddt --workers 20

# 2) Fruit/crop high-score candidates from cross-species deposit (no new species)
#    Default = all species in conservative_candidates.tsv.gz; use --species-limit 8
#    for grape + top 7 by candidate count.
python3 scripts/deposit_candidate_tracks.py --plddt --disorder --species-limit 8 --workers 16
```

Candidate source: `data/fruit_crops/cross_species/conservative_candidates.tsv.gz`
(Tier A∪B in ≥5 panel species; produced by PR #23 / v1.7.0). **No new species.**

## UI

Protein page (`/protein/:id`): metapredict + pLDDT sparklines on the same residue
axis, plus Jaccard cross-check (disorder≥0.5 vs pLDDT&lt;50). JSON:
`/api/protein/:id/tracks`.

## Non-blocking extras

```bash
# InterPro entry names / Tier A∪B locations (optional)
python3 scripts/deposit_interpro_entries.py --protein-locations --tier-ab --workers 12
```

Grape Expression Atlas TPM bars use existing `expr_tpm_json` (E-MTAB-8758) — not
re-fetched here.

## Licenses / sources

| Track | Source | Notes |
|-------|--------|-------|
| metapredict | [metapredict](https://github.com/idptools/metapredict) | Local deposit-time only |
| AlphaFold pLDDT | [AlphaFold DB](https://alphafold.ebi.ac.uk/) | Confidence JSON; AFDB terms |
| InterPro | [InterPro](https://www.ebi.ac.uk/interpro/) | Optional name cache |
| Expression Atlas | E-MTAB-8758 | Grape only; see `data/expression_atlas/README.md` |

Do **not** invent tracks for accessions without a successful deposit. Full fruit
proteome InterPro/AF is **out of scope** — candidates only.
