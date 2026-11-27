# Release notes

## v1.8.0 (2026-11-27)

Orthogonal LLPS database cross-annotation (NAR content layer A1).

- `scripts/fetch_llps_external_ids.py` — fetch/cache PhaSePro JSON, LLPSDB unambiguous PS zip, DrLLPS `LLPS.txt` + *Vitis* gene FASTA, CD-CODE Edmond dump; write compact UniProt/gene ID tables under `data/llps_crosswalk/ids/` with `sources.json` (URLs, citations, license notes).
- `scripts/enrich_llps_crosswalk.py` — cross PlantLLPS gold panel (train ∪ curated ∪ reconstructed; **n=164**) + Tier A/B/C candidates; deposits `gold_crosswalk.tsv|json`, `candidates_crosswalk.tsv.gz`, `summary.json`; SQLite columns `ext_llps_dbs` / `ext_llps_unique` / `ext_llps_json`.
- UI `/crosswalk` + protein-page external DB badges; download/API filter `evidence=unique|shared|phasepro|llpsdb|drllps|cdcode`; `/api/crosswalk`; `/api/health.llps_crosswalk`.
- Methods note [`docs/LLPS_CROSSWALK.md`](./docs/LLPS_CROSSWALK.md) for reviewer Q “how do you differ from existing LLPS DBs”.
- Headline: Tier A∪B **3870** → unique **3396** / shared **474**; gold **164** → unique **100** / shared **64**. Shared dominated by DrLLPS *Vitis* ortholog expansion.
- Honesty: unique = catalog-relative under deposited match rules; CD-CODE IDs are CC BY-NC accession extracts; no Zenodo; no new species.

## v1.7.0 (2026-11-27)

Fruit-crop cross-species high-score comparison layer (no re-scoring; no Zenodo touch).

- `scripts/analyze_fruit_cross_species.py` — reads the 18 scored `data/fruit_crops/*_condensatome.sqlite` deposits; groups near-orthologs by OrthoDB → grape join → gene-symbol proxy → singleton.
- Deposits under `data/fruit_crops/cross_species/`: `conservative_candidates.tsv.gz`, `species_specific_peaks.tsv.gz`, `cluster_members.tsv.gz`, `summary.json`, `browse.json`, methods `README.md`.
- Thresholds (documented): **conservative** = Tier A∪B in ≥5 panel species; **species-specific** = high in ≤2 species (panel-contrast when present in ≥3, else orphan-high).
- UI `/hooks#fruit-cross-species` — filterable conservative / species-specific table (species, conservation rank / peak, min score, search); `/api/fruit-cross-species`; `/api/health.second_crop.fruit_crops.cross_species`.
- Honesty: OrthoDB eukaryote groups + gene proxies are approximations; OrthoDB-sparse proteomes enrich orphan-high peaks; ranks are hypotheses from grape-trained transfer.
- Headline: conservative clusters **1438** (almost all OrthoDB); species-specific **73 336** (panel-contrast **1675**, orphan-high **71 661**); OrthoDB-linked high species in panel **12** / 18.

## v1.6.0 (2026-11-27)

Multi-species fruit-tree / fruit-crop condensatome scoring pipeline.

- `data/fruit_crops/species_manifest.yaml` — pinned UniProt proteome IDs for apple, pear, peach, sweet cherry, citrus (sinensis + clementina), banana, strawberry (*F. ananassa*), blueberry (*V. darrowii*), kiwifruit, avocado, pineapple, melon, date palm, olive, fig, pomegranate, blackberry (*R. argutus*), plus linked grape/soybean/tomato; documented unavailable species (plum, raspberry, *F. vesca*, *V. corymbosum*, mango, papaya, watermelon, pure *P. communis*).
- `scripts/score_fruit_crop_proteomes.py` — bounded parallel FASTA download + PlantLLPS v13 composition-proxy scoring; deposits slim SQLite + summary JSON under `data/fruit_crops/`; `index.json` / `failures.json`.
- `scripts/enrich_fruit_crop_orthologs.py` — OrthoDB fruit→grape joins (+ Ath via grape) written to **crop SQLite only** (no new grape atlas columns; keeps DB ≤100 MiB).
- UI `/hooks` fruit-crop browser (species switcher); `/api/health.second_crop.fruit_crops`.
- Honesty: grape-trained weights — **not** species-retrained; OrthoDB groups can be broad; no invented deposits for missing proteomes.
- Headline (full proteomes scored): apple **39 929** (A/B/C/D 1728/4440/14846/18915); pear **39 683**; peach **37 173**; sweet cherry **30 295**; sweet orange **41 603**; clementine **30 104**; banana **37 776**; strawberry **76 950**; blueberry **34 273**; kiwifruit **32 478**; avocado **40 746**; pineapple **29 488**; melon **28 113**; date palm **39 147**; olive **75 223**; fig **43 146**; pomegranate **47 028**; blackberry **38 432**. Linked: grape **25 972**, soybean **71 064**, tomato **31 194**.
- OrthoDB fruit→grape coverage varies by UniProt xref density (e.g. olive **56 710**, kiwifruit **29 090**, pear **28 696**; several proteomes lack OrthoDB xrefs entirely — joins recorded as 0, not invented).

## v1.5.3 (2026-11-27)

Supplemental Figures S1–S5 for the NAR Database manuscript.

- Generator `scripts/generate_supp_figures.py` (+ pinned `requirements-figures.txt`) deposits PNG/PDF under `docs/figures/supp/` with `LEGENDS.md`.
- S1 stratified CV/gold ROC (in-sample caveat); S2 fair leave-gold-out / family / Ath→non_ath_crop; S3 feature-ablation + F6HET4 SHAP; S4 proteome score histogram; S5 stress DE × Tier A∪B enrichment + modules.
- Linked from `docs/MANUSCRIPT.md` and indexed in `README.md`.

## v1.5.2 (2026-11-27)

Cross-crop conserved condensatome + tomato third-crop atlas.

- `scripts/analyze_cross_crop_conserved.py` — grape ↔ soybean ↔ Ath Tier A∪B / OrthoDB overlap; deposits `data/cross_crop/conserved_tier_ab_overlap.tsv.gz`, `conserved_orthodb_clusters.tsv.gz`, `conserved_summary.json`.
- Ath PlantLLPS deposit (`--crop arabidopsis --out-dir data/cross_crop`) for Ath Tier labels in the overlap.
- Tomato third crop: `scripts/score_second_crop_proteome.py --crop tomato` + `enrich_second_crop_orthologs.py` → `data/second_crop/tomato_condensatome.*`; grape `ortholog_sly*` columns; curated MAP65-1a hook `A0A3Q7HI29` (`cross_crop_hooks.json` v3).
- UI `/hooks` tomato panel + protein-page tomato ortholog; `/api/health.second_crop` reports tomato + conserved summary.
- Stub helper `scripts/write_tomato_schema_stub.py` if proteome download is blocked (not needed this release — UniProt succeeded).
- Headline: grape Tier A∪B **3870**; grape∩soy Tier A/B **968**; triple Tier A/B **271**; grape-specific Tier A/B **1645**; tomato scored **31 194** (A/B/C/D = 1387/3424/11566/14817); tomato OrthoDB matches **14 547 / 25 972**.

## v1.5.1 (2026-11-27)

NAR+ atlas analyses: Tier A∪B domain enrichment, comparator disagreement, structure–score joint.

- `scripts/enrich_family_enrichment.py` — hypergeometric Pfam/InterPro enrichment vs scored proteome; BH FDR within each DB; deposits `data/family_enrichment/` + figure `public/figures/top_families_pfam.svg`.
- `scripts/enrich_comparator_disagreement.py` — PlantLLPS Tier A∪B with native PSAP and/or LLPhyScore in the proteome bottom quartile; Ser-rich/MAP focus flags; `data/comparator_disagreement/` + `METHOD.md`.
- `scripts/enrich_structure_score_joint.py` — Tier A∪B ∩ AlphaFold low-pLDDT (`af_frac_very_low`); threshold hypergeometric + BH FDR; per-protein flags under `data/structure_score/`.
- UI `/analyses` (family / disagreement / structure tabs); `/api/health.atlas_analyses`; SQLite meta pointers only (DB stays ≤100 MiB).
- Headline: Pfam FDR&lt;0.05 **144** / 424 tested; InterPro FDR&lt;0.05 **323** / 984; disagreement **1 555** (both-low **277**; Ser/MAP focus **448**); structure joint flags **1 872** / 3 423 Tier A∪B with AF (all 4 af_frac thresholds FDR&lt;0.05).
- Honesty: family + structure threshold tests are FDR-controlled; disagreement set and per-protein structure flags are exploratory.

## v1.5.0 (2026-11-27)

NAR+ atlas analyses on existing stress_modules / STRING / UniProt joins (analyses 2, 6, 7).

- `scripts/analyze_stress_network_isoforms.py` — one pass depositing:
  1. **Stress module deep dive** — GO-ID + UniProt keyword enrichment for the FDR-significant Botrytis stage-3 down ∩ Tier A∪B set (`E-GEOD-67932` `g1_g4` down); exploratory tests for other contrasts; heatmap-ready log2FC matrix → `data/stress_modules/deepdive/`
  2. **STRING / stress seed hubs** — rank Tier A∪B STRING neighborhood enrichment (min 5 neighbors) and stress log2FC coexpression seed hubs; top-100 TSV → `data/string/stress_seed_hubs.tsv`
  3. **Isoform / multi-protein genes** — PlantLLPS score deltas for UniProt accessions sharing a `VIT_*` locus → `data/isoforms/`
- Honesty: no KEGG table in DB (keywords = functional-class proxy); exploratory contrasts labeled; STRING ≠ LLPS evidence; UP000009183 is nearly 1:1 VIT↔accession (1 multi-accession VIT deposited).
- Meta pointer `atlas_analyses_stress_network_isoforms`; `/api/health` surfaces summary counts. SQLite stays at the GitHub ~100 MiB limit (sidecars).

## v1.4.0 (2026-11-27)

NAR+ explainability: Ser-rich / charge-block / IDR-domain attributions on protein pages.

- `plantllps/explain.py` — exact linear SHAP for the logistic head (logit-space `φ_i = w_i·z_i`), companion term decomposition, grammar axes (Ser-rich, charge-block, IDR/disorder), Gly window occlusion joined to deposited IDR/domain intervals.
- `scripts/deposit_explainability.py` — deposits gold ∪ Tier A sidecar `data/explainability/attributions.json.gz` (**1159** proteins) + pretty `F6HET4.json` example; SQLite meta pointer only (keeps DB ≤100 MiB).
- Node UI: protein pages render grammar SHAP bars + top features + occlusion/IDR attributions when deposited; `/score` shows live feature SHAP via `lib/explain.js`.
- Anchor F6HET4 (Tier B, curated gold): Ser-rich SHAP leads the MAP65 grammar story; occlusion windows tagged ser_rich/charge_block.
- Honesty: logistic SHAP is exact for the linear head; window occlusion is SHAP-equivalent locally, not KernelSHAP/TreeSHAP.

## v1.3.0 (2026-11-27)

NAR+ **fair comparison suite**: leakage-controlled PlantLLPS vs native PSAP / LLPhyScore / PLAAC-approx / ProX-style feature baseline on the same gold, family hold-out, and expanded cross-species panels.

- `scripts/eval_fair_benchmark.py` + `plantllps/prox_features.py`
- Deposits: `data/benchmarks/fair_comparison_metrics.json`, `fair_comparison_scores.tsv`, `fair_comparison_family.tsv`, `train_v13_oob_comparator_scores.tsv`, `FAIR_COMPARISON_README.md`
- Leaderboard: [`docs/BENCHMARK_LEADERBOARD.md`](./docs/BENCHMARK_LEADERBOARD.md)
- Headline (fair): gold leave-out PlantLLPS hybrid ≈ **0.959** (production in-sample ≈ 0.990 caveated); family-hold-out mean hybrid ≈ **0.833**; ProX-style family mean ≈ **0.836** but gold leave-out ≈ **0.467**; OOB LLPhyScore gold ≈ **0.769**, PSAP ≈ **0.718**
- Cross-species expansion: Ath→cereal (rice+maize, n=9), Ath→non_ath_crop (n=38), non_ath_crop→Ath / cereal→Ath; Ath→maize / Ath→soybean skipped (insufficient labels)
- Honesty: ProX/PU-style ≠ published Plant PS proteins LightGBM/PU; PLAAC-approx ≠ Java PLAAC; shared negatives in family protocol documented

## v1.2.0 (2026-11-27)

NAR+ omics layer: condensatome × stress DE enrichment + high-tier response modules.

- `scripts/enrich_stress_condensatome_modules.py` — hypergeometric enrichment of PlantLLPS Tier A∪B among significant DE genes from **E-MTAB-4324** and **E-GEOD-67932** (per contrast × direction; BH FDR); response-pattern modules among Tier A∪B (shared signed DE calls); Tier A seed coexpression neighborhoods (Pearson r on stress log2FC profiles).
- Deposits under `data/stress_modules/`: `de_enrichment.tsv`, `coexpression_modules.tsv`, `seed_neighborhoods.tsv.gz`, `module_members.tsv.gz`, `gene_sets.gmt.gz` (**1 198** gene sets), `ui_tables.json`, `summary.json`.
- UI `/modules` table (enrichment + pattern modules); `/expression` links through; `/api/health.stress_modules`.
- Headline: **18** enrichment tests; **1** at FDR &lt; 0.05 (E-GEOD-67932 stage-3 Botrytis down; fold ≈ 1.17); **42** pattern modules (size ≥ 5); **1 120** seed neighborhoods.
- Honesty: descriptive overlap / short stress log2FC vectors — not WGCNA; berry TPM (E-MTAB-8758) has only 6 assay groups and is not the primary module signal. SQLite stays under GitHub 100 MiB (sidecars).

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

### Deposit / hosting prep (retained from post-v1.1.0)

- `docs/ZENODO.md` + `scripts/export_zenodo_bundle.py` — build Zenodo bundle (weights, SQLite, TSV.gz, optional Parquet, checksums). **No invented DOI**; user publishes on Zenodo and pastes DOI into `CITATION.cff`.
- `docs/HOSTING.md` — clarified **local-server-only** until a real public URL exists; Railway/Render one-liners; SQLite size note (~100 MiB).

## Public host + Zenodo automation (post-v1.1.1)

- Deploy configs: `render.yaml`, `railway.toml`, `fly.toml`, `.dockerignore` (+ existing `Dockerfile`); health path `/api/health`.
- **Public URL live:** https://grape-condensatome.apps.grok.bot (Cursor publish-app). README badges + `CITATION.cff` `url` updated; Zenodo DOI still pending.
- `docs/HOSTING.md` / `docs/ZENODO.md` for **v1.1.1**; export script reads `VERSION`. Security review: resource_exhaustion / data_exposure noted for later hardening.
- Bundle built locally with `--parquet` into gitignored `dist/` (~33 MiB zip); not committed.
- Untracked `data/grape_condensatome.proteome.tsv` + `data/second_crop/*.proteome.tsv.gz` (gitignored) so Cursor publish-app stays under the 50 MiB pack limit; tables remain in the Zenodo bundle / regenerable from SQLite.

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
