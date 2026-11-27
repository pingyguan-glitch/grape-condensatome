# Zenodo deposit — PlantLLPS + GrapeCondensatome

Prepare a **versioned data+model deposit** for Zenodo (weights, SQLite / tabular
exports, checksums, citation metadata). This document does **not** invent a DOI.

Until you publish on Zenodo and paste the real DOI back into the repo, cite the
GitHub tag (e.g. `v1.1.1`) via [`CITATION.cff`](../CITATION.cff).

## What goes in the deposit

| Artifact | Repo path | Role |
|----------|-----------|------|
| PlantLLPS v13 weights | `data/model_formal.json` | Scorable model |
| Training set + metrics | `data/train_formal_v13.csv`, `data/model_formal.metrics.v13.json` | Reproducibility |
| GrapeCondensatome SQLite | `data/grape_condensatome.sqlite` | Full atlas (~100 MiB) |
| Slim proteins table | generated TSV.gz / optional Parquet | Analyst-friendly export (preferred over the gitignored full proteome TSV) |
| Soybean second crop | `data/second_crop/soybean_condensatome.sqlite` (+ summary/join JSON) | Cross-crop layer |
| Citation metadata | `CITATION.cff` | Software citation (**DOI empty until publish**) |
| License | `LICENSE` | MIT |
| Checksums | `SHA256SUMS.txt` (generated) | Integrity |

Optional sidecars (included when present): MAP-screen ρ JSON, curated gold TSV,
family hold-out metrics, STRING enrichment summary. Feature-ablation / experimental
disorder weights stay in the GitHub repo under `data/benchmarks/` and
`data/experimental/` (not required for the core DOI deposit).

## Build the bundle locally

```bash
# From repo root (reads VERSION; currently 1.1.1)
python3 scripts/export_zenodo_bundle.py
# → dist/zenodo_1.1.1/  +  dist/grape-condensatome-zenodo-v1.1.1.{zip,tar.gz}

# Optional Parquet (needs: pip install 'numpy>=2' pyarrow)
python3 scripts/export_zenodo_bundle.py --parquet

# Preview only
python3 scripts/export_zenodo_bundle.py --dry-run
```

### Artifact policy (do not commit archives)

The bundle directory is gitignored (`dist/`). **Do not** commit the archive to GitHub
(GitHub soft limit ~100 MiB; the SQLite alone is already near that ceiling).

Verified locally for **v1.1.1** (Cloud Agent run):

| Artifact | Approx. size | Notes |
|----------|-------------:|-------|
| `dist/zenodo_1.1.1/` | folder | includes SQLite, TSV.gz, Parquet, checksums |
| `dist/grape-condensatome-zenodo-v1.1.1.zip` | ~33 MiB | upload this **or** the folder |
| `dist/grape-condensatome-zenodo-v1.1.1.tar.gz` | ~32 MiB | alternate archive |
| `grape_condensatome.proteins.parquet` | ~3 MiB | present when `--parquet` succeeds |
| Proteins exported | **25 972** | matches atlas |

Re-run the script on your machine before upload if your tree moved. Verify:

```bash
cd dist/zenodo_1.1.1
sha256sum -c SHA256SUMS.txt
```

## Exact steps — create the Zenodo upload (you do this)

Agents / scripts **must not** mint or invent `10.5281/zenodo.…` identifiers.

1. **Log in** at [https://zenodo.org](https://zenodo.org) (ORCID recommended).
2. **New upload** → upload either:
   - the folder contents under `dist/zenodo_1.1.1/`, or
   - `dist/grape-condensatome-zenodo-v1.1.1.zip` (or `.tar.gz`).
3. **Metadata** (suggested):
   - **Title:** PlantLLPS and GrapeCondensatome (v1.1.1)
   - **Upload type:** Dataset (or Software)
   - **Creators:** Pingyin Guan (ORCID `0000-0003-4890-3412`)
   - **Description:** paste the abstract from `CITATION.cff` + note that gold AUCs are partially in-sample; composition-proxy disorder is the production default
   - **Version:** `1.1.1`
   - **License:** MIT
   - **Related identifiers:** GitHub repo  
     `https://github.com/pingyguan-glitch/grape-condensatome` (isSupplementTo / isIdenticalTo as appropriate)  
     and git tag `v1.1.1`
   - **Keywords:** liquid-liquid phase separation; condensatome; Vitis vinifera; PlantLLPS
4. **Save draft** → review files + `SHA256SUMS.txt`.
5. **Publish.** Zenodo assigns a DOI of the form `10.5281/zenodo.<digits>`.
   Copy that DOI from the published record page. **Only this value is real.**

### Concept DOI vs version DOI

Zenodo also shows a **Concept DOI** that always resolves to the latest version.
Prefer citing the **version DOI** for `v1.1.1` in papers; keep the Concept DOI
for “latest release” badges if desired.

## Paste the DOI back into the repo

After publish, edit [`CITATION.cff`](../CITATION.cff) and uncomment / set:

```yaml
identifiers:
  - type: doi
    value: "10.5281/zenodo.PASTE_REAL_ID_HERE"
    description: "Zenodo deposit for PlantLLPS + GrapeCondensatome v1.1.1"
```

Also update:

1. Project store honesty line in `docs/grape-condensatome-prototype.md` (DOI pending → DOI set).
2. Optionally add the DOI badge / link to the root `README.md`.
3. Commit, tag a patch release if you want the DOI in a git tag (e.g. `v1.1.2`), or amend release notes.

**Placeholder DOIs, example DOIs, or `10.5281/zenodo.XXXXXXX` must never be committed as real identifiers.**

## What this prep does *not* do

- Does not create a Zenodo account or publish a record.
- Does not invent a DOI (public atlas URL is separate — see [`HOSTING.md`](./HOSTING.md): https://grape-condensatome.apps.grok.bot).
- Does not enlarge GitHub with the 100 MiB+ archive.

## Related

- Hosting status: [`HOSTING.md`](./HOSTING.md) (local server unless a real URL is live)
- External tool pins: [`EXTERNAL_TOOLS.md`](./EXTERNAL_TOOLS.md)
