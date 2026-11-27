# Zenodo deposit — PlantLLPS + GrapeCondensatome

Prepare a **versioned data+model deposit** for Zenodo (weights, SQLite / tabular
exports, checksums, citation metadata). This document does **not** invent a DOI.

Until you publish on Zenodo and paste the real DOI back into the repo, cite the
GitHub tag (e.g. `v1.1.0`) via [`CITATION.cff`](../CITATION.cff).

## What goes in the deposit

| Artifact | Repo path | Role |
|----------|-----------|------|
| PlantLLPS v13 weights | `data/model_formal.json` | Scorable model |
| Training set + metrics | `data/train_formal_v13.csv`, `data/model_formal.metrics.v13.json` | Reproducibility |
| GrapeCondensatome SQLite | `data/grape_condensatome.sqlite` | Full atlas (~100 MiB) |
| Slim proteins table | generated TSV.gz / optional Parquet | Analyst-friendly export |
| Soybean second crop | `data/second_crop/soybean_condensatome.sqlite` (+ summary/join JSON) | Cross-crop layer |
| Citation metadata | `CITATION.cff` | Software citation (**DOI empty until publish**) |
| License | `LICENSE` | MIT |
| Checksums | `SHA256SUMS.txt` (generated) | Integrity |

Optional sidecars (included when present): MAP-screen ρ JSON, curated gold TSV,
family hold-out metrics, STRING enrichment summary.

## Build the bundle locally

```bash
# From repo root (v1.1.0 tree)
python3 scripts/export_zenodo_bundle.py
# → dist/zenodo_1.1.0/  +  dist/grape-condensatome-zenodo-v1.1.0.{zip,tar.gz}

# Optional Parquet (needs: pip install pyarrow)
python3 scripts/export_zenodo_bundle.py --parquet

# Preview only
python3 scripts/export_zenodo_bundle.py --dry-run
```

The bundle directory is gitignored (`dist/`). Do not commit the archive to GitHub
(GitHub soft limit ~100 MiB; the SQLite alone is already near that ceiling).

Verify checksums after unpack:

```bash
cd dist/zenodo_1.1.0
sha256sum -c SHA256SUMS.txt
```

## Exact steps — create the Zenodo upload (you do this)

Agents / scripts **must not** mint or invent `10.5281/zenodo.…` identifiers.

1. **Log in** at [https://zenodo.org](https://zenodo.org) (ORCID recommended).
2. **New upload** → upload either:
   - the folder contents under `dist/zenodo_1.1.0/`, or
   - `dist/grape-condensatome-zenodo-v1.1.0.zip` (or `.tar.gz`).
3. **Metadata** (suggested):
   - **Title:** PlantLLPS and GrapeCondensatome (v1.1.0)
   - **Upload type:** Dataset (or Software)
   - **Creators:** Pingyin Guan (ORCID `0000-0003-4890-3412`)
   - **Description:** paste the abstract from `CITATION.cff` + note that gold AUCs are partially in-sample
   - **Version:** `1.1.0`
   - **License:** MIT
   - **Related identifiers:** GitHub repo  
     `https://github.com/pingyguan-glitch/grape-condensatome` (isSupplementTo / isIdenticalTo as appropriate)  
     and git tag `v1.1.0`
   - **Keywords:** liquid-liquid phase separation; condensatome; Vitis vinifera; PlantLLPS
4. **Save draft** → review files + `SHA256SUMS.txt`.
5. **Publish.** Zenodo assigns a DOI of the form `10.5281/zenodo.<digits>`.
   Copy that DOI from the published record page. **Only this value is real.**

### Concept DOI vs version DOI

Zenodo also shows a **Concept DOI** that always resolves to the latest version.
Prefer citing the **version DOI** for `v1.1.0` in papers; keep the Concept DOI
for “latest release” badges if desired.

## Paste the DOI back into the repo

After publish, edit [`CITATION.cff`](../CITATION.cff) and add:

```yaml
identifiers:
  - type: doi
    value: "10.5281/zenodo.PASTE_REAL_ID_HERE"
    description: "Zenodo deposit for PlantLLPS + GrapeCondensatome v1.1.0"
```

Also update:

1. `scripts/export_zenodo_bundle.py` output is regenerated only if you re-deposit.
2. Project store honesty line in `docs/grape-condensatome-prototype.md` (DOI pending → DOI set).
3. Optionally add the DOI badge / link to the root `README.md`.
4. Commit, tag a patch release if you want the DOI in a git tag (e.g. `v1.1.1`), or amend release notes.

**Placeholder DOIs, example DOIs, or `10.5281/zenodo.XXXXXXX` must never be committed as real identifiers.**

## What this prep PR does *not* do

- Does not create a Zenodo account or publish a record.
- Does not claim a public hosted atlas URL (see [`HOSTING.md`](./HOSTING.md)).
- Does not enlarge GitHub with the 100 MiB+ archive.

## Related

- Hosting status: [`HOSTING.md`](./HOSTING.md) (local server unless a real URL is live)
- External tool pins: [`EXTERNAL_TOOLS.md`](./EXTERNAL_TOOLS.md)
