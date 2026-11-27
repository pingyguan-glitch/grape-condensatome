# Release notes

## v1.1.1 (2026-11-27)

Figure 3 feature-ablation panels + optional experimental disorder retrains (composition production weights unchanged).

### Feature ablation (Figure 3)

- `scripts/eval_feature_ablation.py` — leave-one-group-out + group-only stratified 5-fold CV on v13 (same seed as production).
- Deposits: `data/benchmarks/feature_ablation_metrics.{json,tsv}`, `feature_ablation_panel.svg`, `FEATURE_ABLATION_README.md`.
- Headline (composition proxy): full-model mean logistic AUC ≈ **0.863** / hybrid ≈ **0.859**; largest LOGO drop **map65_cterm** (Δlogistic ≈ **−0.036**), then **disorder** (≈ **−0.021**). Family hold-out remains the transfer probe (`family_holdout_metrics.json`).

### Experimental disorder retrains (optional)

- `scripts/retrain_disorder_features.py` — comparable CV with metapredict and/or AFDB `fractionPlddtVeryLow` replacing `frac_disorder_proxy`.
- **Does not overwrite** `data/model_formal.json`. Experimental weights under `data/experimental/`.
- Side-by-side (mean logistic AUC / Δ vs composition): composition **0.863**; metapredict **0.882** (+0.020); AF pLDDT **0.862** (−0.000); metapredict+AF mean **0.877** (+0.014). AF coverage **156 / 162** (6 rows fall back to composition).
- Install notes unchanged in `docs/EXTERNAL_TOOLS.md` (`METAPREDICT_PYTHON`).

## Deposit / hosting prep (post-v1.1.0)

- `docs/ZENODO.md` + `scripts/export_zenodo_bundle.py` — build Zenodo bundle (weights, SQLite, TSV.gz, optional Parquet, checksums). **No invented DOI**; user publishes on Zenodo and pastes DOI into `CITATION.cff`.
- `docs/HOSTING.md` — clarified **local-server-only** until a real public URL exists; Railway/Render one-liners; SQLite size note (~100 MiB).

## v1.1.0 (2026-11-27)

P1 atlas layers on v1.0.0: full scored-proteome AlphaFold low-pLDDT segments + multi-stress RNA-seq.

### AlphaFold low-pLDDT segments (scored proteome)

- `scripts/enrich_alphafold_plddt.py` — `--segments-tier ALL` (default) + `--segments-only` resume; browser-like User-Agent for AFDB file CDN (403 otherwise); cache under `data/alphafold/`.
- Coverage on scored proteome (**25 972**): **23 421** with mean pLDDT; **23 421 / 23 421** AFDB models with segment residue counts (**22 449** nonempty low-pLDDT). Remaining **2 551** proteins have no AFDB model.
- Protein page shows segments alongside domain / IDR tracks; IDR refresh → **23 590** proteins with IDR intervals (AF∪MobiDB; was 9 401 when segments were gold+Tier A only).
- Coverage written to `meta.alphafold_plddt`. AF interval geometry ships via the IDR track (`source=af_plddt_lt50`) plus `af_n_disorder_residues` so SQLite stays under GitHub's 100 MiB limit.

### Multi-stress transcriptome

- `scripts/enrich_stress_rnaseq.py` joins Expression Atlas RNA-seq differentials:
  - **E-MTAB-4324** — downy mildew / Trichoderma induced resistance (leaves)
  - **E-GEOD-67932** — Botrytis noble-rot stages vs control (berries)
- Both match **25 810 / 25 972** (99.4%) on `VIT_*` gene_name (same key as baseline E-MTAB-8758).
- UI `/expression` + protein-page contrast cards (log2FC / p; significant when p&lt;0.05 and |log2FC|≥1).
- DE blobs ship as `data/expression_atlas/stress_by_protein.json.gz` (Node loads at startup).
- Honesty: abiotic drought/salt Atlas grape experiments (E-GEOD-31594 / E-GEOD-31677) are microarray — not claimed as RNA-seq.

## v1.0.0 (2026-11-27)

Combined Wave 2 landing on `main`: curated LLPS gold / MAP hold-out (**v0.8.0**), STRING neighborhoods (**v0.9.0**), and soybean second-crop condensatome (**v0.9.1** content). Single coherent atlas release — not a claim of production-complete PlantLLPS.

## v0.9.1 (2026-11-27)

Second-crop condensatome: **soybean** (*Glycine max*, UniProt **UP000008827**).

Originally a parallel additive layer; now merged into **v1.0.0** with gold/MAP + STRING.

- `scripts/score_second_crop_proteome.py` — download + PlantLLPS v13 score (composition-proxy); writes `data/second_crop/soybean_condensatome.{proteome.tsv.gz,sqlite,summary.json}`. Tomato (`--crop tomato`, UP000004994) supported by the same script.
- **Full proteome scored:** **71 064** proteins (min length 80); tiers A=3 137 / B=7 983 / C=26 235 / D=33 709; ~3.5 min wall on a 4-core CPU.
- `scripts/enrich_second_crop_orthologs.py` — OrthoDB join soybean↔grape; grape columns `ortholog_gma*` + crop DB reverse grape pointers.
- Join coverage: **16 388 / 25 972** grape proteins with a soybean OrthoDB hit (**16 304** with crop scores); curated MAP65-1 soybean IDs overlay F6HET4 → `A0A0R0GGZ9` (Tier B, score ≈ 0.805).
- UI `/hooks` expanded: soybean condensatome table (tier filter + pager) + soybean columns on family hooks / protein pages; `/api/health` reports `second_crop`.
- Honesty: grape-trained weights, not soybean-retrained; OrthoDB groups can be broad; FASTA cache gitignored.

## v0.9.0 (2026-11-27)

Additive **STRING neighborhood + condensatome enrichment** on top of Wave 2 **v0.8.0** (gold/MAP/hold-out).

- `scripts/enrich_string_neighborhood.py` — native *Vitis* STRING (taxid 29760) v12 downloads;
  UniProt mapping **25 867 / 25 972** (99.6%); 1-hop neighborhoods at combined_score ≥ 400.
- SQLite: `string_neighborhood` for **Tier A seeds** (~833 deposited seed rows) + `string_id` on mapped proteins
  (full-proteome neighbor JSON would exceed GitHub's 100 MB SQLite limit); combined DB also keeps curated LLPS gold columns from v0.8.0.
- UI: protein-page neighborhood snippet + `/network` enrichment table; `/api/health` reports STRING meta.
- Deposits: `data/string/enrichment_summary.json`, `neighborhood_stats.tsv.gz` (19 523),
  `top_neighbors.tsv.gz`, `seed_edges.tsv.gz` (44 932), `README.md`.
- Honest limits: combined score ≠ LLPS evidence (includes textmining); no Ath PPI transfer;
  seeds = score Tier A (gold overlay is separate).

## v0.8.0 (2026-11-27)

Wave 2 resource-quality: curated LLPS gold, MAP-screen deposit + Spearman ρ, family / transfer AUCs.

- **Curated LLPS gold** — `data/curated/llps_experimental_gold.tsv` (protein ID, species, condensate type, assay/evidence, citation, grape ortholog); `scripts/enrich_curated_llps_gold.py` joins SQLite `llps_gold` / `llps_gold_json`; UI `/gold` + protein-page badge + download `evidence=llps`.
- **MAP screen deposit** — reconstructed GO:0008017 ∪ KW-Microtubule table under `data/map_screen/` with provenance README; recomputed Spearman ρ vs PlantLLPS / native PSAP / LLPhyScore (`map_screen_spearman.json`). Companion-cited ρ=0.368 retained as manuscript number (original 208-ID PSAP table still not recovered).
- **Family hold-out / transfer** — `scripts/eval_family_holdout.py` → `data/benchmarks/family_holdout_metrics.json` + README; mean hybrid family-hold-out AUC ≈ **0.833** (13 families); Ath→rice hybrid AUC = **1.00** on tiny n=6 (honesty: small rice slice). Contrasts in-sample gold hybrid AUC ≈ 0.990.
- Version **0.8.0** (skips 0.7.x to avoid collision with in-flight P1 AF-segments / stress RNA-seq work).

## v0.6.0 (2026-11-27)

P0 atlas layers: domains / IDR tracks, bulk download, genome link-outs.

- `scripts/enrich_domains_idr.py` — UniProt Domain FT + InterPro Pfam coordinate backfill; IDR from AF pLDDT&lt;50 ∪ UniProt REGION Disordered (MobiDB-lite). Coverage on scored proteome: **21 856** with domain coords, **21 670** with Pfam xref, **9 401** with IDR segments.
- Protein page: domain/IDR track visualization + Ensembl Plants / location (JBrowse) links from `VIT_*` gene names.
- `GET /api/download` + UI `/download` — filter by tier / family / keyword / q; TSV or FASTA.
- Honest limit: metapredict remains mean-fraction only (no per-residue IDR); splice isoforms not annotated.

## v0.5.0 (2026-11-26)

Full-proteome **native** PSAP and LLPhyScore (user priority: comparator columns).

- `scripts/score_psap_proteome.py` — batched + checkpointed native PSAP (Python ≤3.9 env).
- `scripts/score_llphyscore_proteome.py` — batched + checkpointed LLPhyScore standalone (8-feature-sum human+PDB %ile).
- `scripts/enrich_comparator_scores.py` joins proteome TSVs into SQLite; UI `/compare` shows coverage + proteome top-40 with real PSAP / LLPhyScore values.
- Honest docs in `docs/EXTERNAL_TOOLS.md` (env pins, wall-time expectations, resumable checkpoints).

## v0.4.0 (2026-11-26)

Atlas enrichments for user priorities 1–3 + 6 (honest compute scope).

### AlphaFold pLDDT

- `scripts/enrich_alphafold_plddt.py` — API mean pLDDT + fraction bins into SQLite; optional low-pLDDT segments (gold + Tier A by default).
- UI protein pages show AF mean / frac&lt;50 and segment lists when present.
- Network-bound; cache under `data/alphafold/` (gitignored). Full proteome summaries are resumable.

### metapredict proteome

- `scripts/score_metapredict_proteome.py` — batched, checkpointed scoring of DB IDs; joins `frac_disorder_metapredict`.
- On a warm 4-core CPU, ~5–15 min for ~26k sequences once metapredict is installed; cold starts / memory pressure can be longer — use `--limit` / `--batch-size`.

### Orthologs + cross-crop hooks

- `scripts/enrich_orthologs.py` — UniProt OrthoDB streams for grape / Arabidopsis / rice; curated `data/cross_crop_hooks.json`.
- UI `/hooks` + Ath/rice columns in tables.

### Side-by-side comparators

- `scripts/enrich_comparator_scores.py` + `plantllps/plaac_approx.py`.
- Gold PSAP / LLPhyScore joined; PLAAC-approx for full proteome; UI `/compare`.
- Full-proteome native PSAP/LLPhyScore deferred to **v0.5.0** (env + CPU).

## v0.3.1 (2026-11-26)

Submission-readiness follow-on (base: v0.3.0 / PR #1 branch).

### metapredict

- Optional `--use-metapredict` on `score` / `score-proteome` (+ `METAPREDICT_PYTHON`).
- `plantllps/metapredict_bridge.py` + `scripts/score_metapredict_subset.py` (29 sequences scored; mean frac ≈ 0.35).
- Install / pin notes in `docs/EXTERNAL_TOOLS.md`.

### Native PSAP

- Runnable in Python 3.9 env (`numpy<1.24`, `scipy<1.11`); `scripts/run_psap_gold.py`.
- Gold-standard AUC ≈ **0.718** (companion cited 0.692). Side-by-side in `data/benchmarks/gold_standard_comparison.json`.
- MAP-screen ρ still blocked (208-protein table missing).

### Expression Atlas

- `scripts/enrich_expression_atlas.py` joins **E-MTAB-8758** TPMs → SQLite (`expr_*` columns).
- **25 810 / 25 972** proteins matched (99.4%); Node + Python UIs show mean TPM + per-condition.

### Hosting

- `docs/HOSTING.md`: Cursor publish-app on `:8765`; GitHub Pages unsuitable; Railway/Render/Fly alternatives.

## v0.3.0 (2026-11-26)

Submission-oriented resource packaging (still composition-proxy disorder; not metapredict).

### Model / training

- **v13** training set: v12 + **19 curated disordered non-LLPS hard negatives** (Arabidopsis extensins / PRPs / AGPs / LRX / linker histones + human prothymosin-α).
- Canonical files: `data/train_formal_v13.csv` (mirrored as `data/train_formal.csv`), `data/model_formal.json`.
- Metric delta vs v12 (stratified 5-fold, same code path):

| Metric | v12 | v13 | Δ |
|--------|----:|----:|--:|
| n (pos/neg) | 143 (76/67) | 162 (76/86) | +19 neg |
| mean AUC logistic | 0.859 | 0.863 | +0.004 |
| OOF AUC logistic | 0.853 | 0.868 | +0.015 |
| mean AUC companion | 0.861 | 0.774 | −0.088 |

Companion composite AUC drops as expected when highly disordered non-LLPS proteins enter the negative class; logistic holds / slightly improves — the intended hard-negative effect.

### Benchmarks

- Reconstructed companion **gold standard** (13 pos / 15 neg); scored with PlantLLPS + **LLPhyScore** standalone 20220627.
- LLPhyScore gold AUC ≈ **0.769** (companion cited 0.774).
- **PSAP** not runnable on Python 3.12 in this environment (sklearn-json / sklearn~=0.21); companion PSAP gold AUC **0.692** and MAP-screen ρ **0.368** cited, not re-run.
- Grapevine MAP screen table with PSAP columns **absent** from store/repo — MAP ρ not recomputed.

### Atlas / UI

- Proteome re-scored with v13: **25 972** proteins; tiers A=1142, B=2728, C=9653, D=12449.
- UniProt gene names / GO IDs / keywords joined for all scored rows.
- Node + Python web UIs show gene, protein name, GO, keywords; search includes those fields.
- Anchor: VvMAP65-1 `F6HET4` / `VIT_19s0090g01230` remains **Tier B** (score ≈ 0.762).

### Packaging

- Dual deliverable: Python `plantllps` CLI + Express browser.
- `LICENSE` (MIT), `CITATION.cff`, expanded `README.md`.

## v0.2.0 (prior)

Web app publish snapshot with v12 proteome scores (Node-only repo start).
