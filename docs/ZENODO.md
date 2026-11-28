# Zenodo deposit — PlantLLPS + GrapeCondensatome

Prepare a **versioned code + weights + atlas deposit** for Zenodo (scoring-server
snapshot, sample FASTA, model weights, SQLite / tabular exports, checksums,
citation metadata). This document does **not** invent a DOI.

Until you publish on Zenodo and paste the real DOI back into the repo, cite the
GitHub tip (e.g. `v1.11.0`) via [`CITATION.cff`](../CITATION.cff).

## What goes in the deposit

| Artifact | Repo path | Role |
|----------|-----------|------|
| Scoring-server snapshot | `server.js`, `lib/`, `plantllps/`, `package.json`, … | NAR Web Server runnable code |
| Sample FASTA | `public/examples/score_demo.fasta` | Paste-FASTA `/score` demo |
| Web scoring methods | `docs/WEB_SCORING.md` | Caps / parity notes |
| PlantLLPS v13 weights | `data/model_formal.json` | Scorable model |
| Training set + metrics | `data/train_formal_v13.csv`, `data/model_formal.metrics.v13.json` | Reproducibility |
| GrapeCondensatome SQLite | `data/grape_condensatome.sqlite` | Full atlas (~100 MiB) |
| Slim proteins table | generated TSV.gz / optional Parquet | Analyst-friendly export |
| Soybean / tomato / fruit crops | `data/second_crop/`, `data/fruit_crops/` (when present) | Cross-crop layer |
| Citation metadata | `CITATION.cff`, `VERSION` | Software citation (**DOI empty until publish**) |
| License | `LICENSE` | MIT |
| Checksums | `SHA256SUMS.txt` (generated) | Integrity |

Optional sidecars (included when present): MAP-screen ρ JSON, curated gold TSV,
family hold-out metrics, STRING enrichment summary, LLPS / RBP crosswalks.
Feature-ablation / experimental disorder weights stay in the GitHub repo under
`data/benchmarks/` and `data/experimental/` (not required for the core DOI deposit).

## Build the bundle locally

```bash
# From repo root (reads VERSION; tip 1.11.0)
python3 scripts/export_zenodo_bundle.py
# → dist/zenodo_1.11.0/  +  dist/grape-condensatome-zenodo-v1.11.0.{zip,tar.gz}

# Optional Parquet (needs: pip install 'numpy>=2' pyarrow)
python3 scripts/export_zenodo_bundle.py --parquet

# Preview only
python3 scripts/export_zenodo_bundle.py --dry-run
```

### Artifact policy (do not commit archives)

The bundle directory is gitignored (`dist/`). **Do not** commit the archive to GitHub
(GitHub soft limit ~100 MiB; the SQLite alone is already near that ceiling).

Verified locally for **v1.11.0** (Cloud Agent run):

| Artifact | Approx. size | Notes |
|----------|-------------:|-------|
| `dist/zenodo_1.11.0/` | folder | weights + SQLite + `code/` + sample FASTA + checksums |
| `dist/grape-condensatome-zenodo-v1.11.0.zip` | ~92 MiB | upload this **or** the folder |
| `dist/grape-condensatome-zenodo-v1.11.0.tar.gz` | ~91 MiB | alternate archive |
| `examples/score_demo.fasta` | ~3 KiB | NAR Web Server sample input |
| `code/` | scoring snapshot | Node Express + `plantllps` (copy weights/SQLite into `code/data/` to run) |
| Proteins exported | **25 972** | matches atlas |

Re-run the script on your machine before upload if your tree moved. Verify:

```bash
cd dist/zenodo_1.11.0
sha256sum -c SHA256SUMS.txt
```

## Exact steps — mint the DOI (you do this once)

Agents **must not** invent `10.5281/zenodo.…` identifiers. Minting needs **your**
Zenodo account (or a personal access token you never paste into chat).

### Option A — browser (one click after upload)

1. **Log in** at [https://zenodo.org](https://zenodo.org) (ORCID recommended).
2. **New upload** → upload either:
   - the folder contents under `dist/zenodo_1.11.0/`, or
   - `dist/grape-condensatome-zenodo-v1.11.0.zip` (or `.tar.gz`).
3. **Metadata** (suggested):
   - **Title:** PlantLLPS and GrapeCondensatome (v1.11.0)
   - **Upload type:** Software (or Dataset)
   - **Creators:** Pingyin Guan (ORCID `0000-0003-4890-3412`)
   - **Description:** paste the abstract from `CITATION.cff` + note paste-FASTA `/score`, sample FASTA, and that gold AUCs are partially in-sample; composition-proxy disorder is the production default
   - **Version:** `1.11.0`
   - **License:** MIT
   - **Related identifiers:** GitHub repo  
     `https://github.com/pingyguan-glitch/grape-condensatome`  
     and tip / tag `v1.11.0`
   - **Keywords:** liquid-liquid phase separation; condensatome; Vitis vinifera; PlantLLPS; web server
4. **Save draft** → review files + `SHA256SUMS.txt`.
5. **Publish** (the one click). Zenodo assigns a DOI of the form `10.5281/zenodo.<digits>`.
   Copy that DOI from the published record page. **Only this value is real.**

### Option B — API helper (token required)

Create a personal access token at
[Zenodo applications](https://zenodo.org/account/settings/applications/)
with scopes `deposit:write` and `deposit:actions`. Then:

```bash
export ZENODO_ACCESS_TOKEN=…   # or ZENODO_TOKEN — never commit this
python3 scripts/export_zenodo_bundle.py
# Draft only (review on Zenodo, then Publish in the UI):
python3 scripts/zenodo_publish.py --archive dist/grape-condensatome-zenodo-v1.11.0.zip
# Or mint immediately:
python3 scripts/zenodo_publish.py --archive dist/grape-condensatome-zenodo-v1.11.0.zip --publish
```

If the token is missing, the script exits with code 2 and prints the path of the
ready archive — use Option A.

Cloud Agents in this Project currently have **no** `ZENODO_ACCESS_TOKEN` /
`ZENODO_TOKEN` secret, so they stop after building the bundle and documenting
these steps.

### Concept DOI vs version DOI

Zenodo also shows a **Concept DOI** that always resolves to the latest version.
Prefer citing the **version DOI** for `v1.11.0` in papers; keep the Concept DOI
for “latest release” badges if desired.

## Paste the DOI back into the repo

After publish, edit [`CITATION.cff`](../CITATION.cff) and uncomment / set:

```yaml
identifiers:
  - type: doi
    value: "10.5281/zenodo.PASTE_REAL_ID_HERE"
    description: "Zenodo deposit for PlantLLPS + GrapeCondensatome v1.11.0"
```

Also update:

1. Project store honesty line in `docs/grape-condensatome-prototype.md` (DOI pending → DOI set), if used.
2. Root `README.md` Zenodo badge (pending → real DOI).
3. `docs/MANUSCRIPT.md` `[ZENODO_DOI — pending]` placeholders.
4. Commit, tag a patch release if you want the DOI in a git tag (e.g. `v1.11.1`), or amend release notes.

**Placeholder DOIs, example DOIs, or `10.5281/zenodo.XXXXXXX` must never be committed as real identifiers.**

## What this prep does *not* do

- Does not create a Zenodo account or publish a record without your token / click.
- Does not invent a DOI (public atlas URL is separate — see [`HOSTING.md`](./HOSTING.md): https://grape-condensatome.apps.grok.bot).
- Does not enlarge GitHub with the 100 MiB+ archive.

## Related

- Paste-FASTA scoring: [`WEB_SCORING.md`](./WEB_SCORING.md)
- Hosting status: [`HOSTING.md`](./HOSTING.md)
- External tool pins: [`EXTERNAL_TOOLS.md`](./EXTERNAL_TOOLS.md)
