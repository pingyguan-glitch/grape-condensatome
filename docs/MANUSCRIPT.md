# PlantLLPS and GrapeCondensatome: a plant-trained LLPS predictor and grapevine whole-proteome condensate atlas

Pingyin Guan*

College of Horticulture, China Agricultural University, Beijing 100193, China.

\* To whom correspondence should be addressed. Email: pyguan@cau.edu.cn; ORCID: 0000-0003-4890-3412

---

> **Draft status (remove before submission).** Software/atlas facts reflect GitHub **v1.10.0** (https://github.com/pingyguan-glitch/grape-condensatome; `main` tip `2818e2e319fe587b44c4803042bdcfcea05e273a`). Companion MAP/LLPS manuscript: Guan, working title *A conserved Ser-rich tail marks the microtubule-bundling protein VvMAP65-1 as a phase-separation candidate* — **final journal, year, and DOI pending** (do not invent). **Public URL:** https://grape-condensatome.apps.grok.bot. **Zenodo DOI:** `[ZENODO_DOI — pending]`.

---

## Abstract

Sequence-based predictors of liquid–liquid phase separation (LLPS) are dominated by animal training sets and transfer only partially to plant proteins. Human-trained PSAP and the biophysical scorer LLPhyScore under-rank Ser-rich, aromatic-poor MAP65-1 tails that drive plant microtubule-bundling condensates, while existing plant portals (Plant PS proteins / ProX) omit grapevine (*Vitis vinifera*). Here we present **PlantLLPS**, a plant-enriched LLPS propensity model (formal corpus **v13**; *n* = 162), and **GrapeCondensatome**, a dual software-plus-atlas resource that scores UniProt proteome UP000009183 (**25 972** proteins; tiers A–D) and joins expression, stress RNA-seq modules, atlas analyses, fruit-crop and tomato/soybean condensatomes, orthogonal LLPS-database and RBP/RG crosswalks, and dual-rail residue tracks (metapredict disorder + AFDB pLDDT). On a leakage-controlled fair suite, leave-gold-out PlantLLPS hybrid AUC ≈ **0.959** versus out-of-the-box LLPhyScore ≈ **0.769**, PSAP ≈ **0.718**, PLAAC-approx ≈ **0.497** and a ProX-style feature baseline ≈ **0.467**; family hold-out mean hybrid ≈ **0.833**. Exact logistic SHAP attributes VvMAP65-1 (F6HET4, **Tier B**) primarily to Ser-rich / MAP65 and charge-block grammars. The MIT-licensed CLI and SQLite-backed UI ship as **v1.10.0**, with public atlas https://grape-condensatome.apps.grok.bot. Together, PlantLLPS and GrapeCondensatome fill a specialty-crop niche that general LLPS databases and Plant PS/ProX do not occupy.

## INTRODUCTION

Liquid–liquid phase separation (LLPS) organizes membrane-less compartments that shape plant stress sensing, development and cytoskeletal dynamics (Banani et al. 2017; Shin and Brangwynne 2017; Jung et al. 2020; Dorone et al. 2021; Zhu et al. 2022; Wang et al. 2022; Liu et al. 2023). Sequence grammars—aromatic stickers, charged blocks, low-complexity and related features—are computationally accessible (Das and Pappu 2013; Vernon et al. 2018; Saar et al. 2021; Cai et al. 2022), and community databases such as LLPSDB, PhaSePro, DrLLPS, PhaSepDB and CD-CODE curate experimental condensates across eukaryotes (Li et al. 2020; Mészáros et al. 2020; Ning et al. 2020). Machine-learning tools including PSAP, FuzDrop, catGRANULE, DeePhase and LLPhyScore score proteomes for condensate likelihood (Bolognesi et al. 2016; Hardenberg et al. 2020; van Mierlo et al. 2021; Saar et al. 2021; Cai et al. 2022; Chu et al. 2022). Most of these models were trained or calibrated primarily on animal and yeast examples.

Plant condensate resources have advanced in parallel. **Plant PS proteins / ProX** (Zhang et al. 2023) provides the strongest plant-native competitor: β-isox enrichment lists and predicted proteomes across multiple model and crop species, with a public web database. General LLPS knowledgebases and DrLLPS-style ortholog expansion further annotate plant proteins, but typically by propagating animal/yeast scaffolds rather than plant-trained propensity. Critically, **grapevine (*Vitis*) is absent from Plant PS/ProX species tables**, and no existing portal combines a whole-*Vitis* scored condensatome with a scorer tuned to plant MAP grammars, stress-module enrichment, multi-fruit-crop transfer layers, orthogonal catalog crosswalks and explainable Ser-rich / charge-block attributions.

That gap is biologically consequential. A companion MAP/LLPS screen (Guan, companion manuscript — final journal/year/DOI pending) showed that PSAP transfers only moderately to grapevine microtubule-associated proteins (companion-cited Spearman ρ = 0.368) and that LLPhyScore misses experimentally proven Ser-rich MAP65-1 tails (human-proteome percentiles ~26–30), despite AtMAP65-1 forming droplets at nanomolar concentration (Sahu et al. 2023). Plant-specific scaffolds such as WVD2/WDL are likewise absent from animal training sets. These observations motivate a dual deliverable: (i) **PlantLLPS**, a plant-enriched predictor sensitive to MAP65-1-class composition, and (ii) **GrapeCondensatome**, a genome-wide grapevine atlas packaged as reusable software plus a browsable database—positioned as a specialty-crop resource complementary to, not a replacement for, Plant PS/ProX.

## RESOURCE DESCRIPTION

### Design goals and positioning

GrapeCondensatome is built around design goals that define its niche relative to Plant PS/ProX and NAR-class LLPS databases:

1. **Grapevine-first proteome coverage** — score and tier the full UniProt *V. vinifera* reference proteome (UP000009183), a specialty crop absent from Plant PS/ProX.
2. **Plant/MAP-grammar-aware scoring** — train and ship PlantLLPS on plant-enriched labels with features that recover Ser-rich, aromatic-poor drivers that PSAP/LLPhyScore bury.
3. **Physiology- and structure-linked atlas** — join Expression Atlas TPM, multi-stress RNA-seq with condensatome enrichment modules, atlas analyses (family enrichment, comparator disagreement, structure–score), metapredict / AlphaFold tracks, Pfam/InterPro + IDR intervals, OrthoDB Ath/rice/soybean/tomato hooks, native STRING neighborhoods, and side-by-side PSAP/LLPhyScore/PLAAC-approx columns.
4. **Fair benchmarks and explainability** — deposit a leakage-controlled leaderboard versus ProX-style / PSAP / LLPhy / PLAAC baselines, plus SHAP-equivalent Ser-rich / charge-block / IDR grammar attributions on priority proteins.
5. **Cross-crop and fruit-crop transfer** — deposit soybean and tomato second/third-crop condensatomes, grape↔soy↔Ath conserved Tier A∪B overlap, 18 fruit-crop scored proteomes (grape-trained weights), and a conservative vs species-specific high-score comparison layer.
6. **Orthogonal catalog and nucleic-acid scope** — crosswalk Tier A/B/C against PhaSePro / LLPSDB / DrLLPS / CD-CODE (unique vs shared); annotate RBP2GO / UniProt GO–KW / RG–RGG nucleic-acid triage flags; ship dual-rail residue trajectories (metapredict + AFDB pLDDT) with honest coverage.

The resource does **not** claim to displace Plant PS/ProX as the multi-species plant LLPS portal, nor to replace curated knowledgebases (LLPSDB, PhaSePro, PhaSepDB, CD-CODE, DrLLPS). Its value is the *Vitis* whole-proteome condensatome plus an open, versioned plant-trained scorer with comparator honesty.

### PlantLLPS training corpus and model

Plant positives were compiled from LLPSDB plant entries (Li et al. 2020; freeze: **11** plant PS + **8** plant no-PS rows; Supplemental Table S1) and companion plant drivers (AtMAP65-1, ELF3, FLOE1, RBGD2, SEUSS, FCA, and related mechanism-matched cases). Mechanism-mismatched folded-domain condensates (prototype: CRY2) are held out of IDR-driver training. Negatives comprise folded globular controls, length-matched unlabeled sequences, and a **v13** hard set of disordered non-LLPS proteins (extensins, PRPs, AGPs, LRX, linker histones, prothymosin-α-class). After ~90% identity deduplication, the released **v13** corpus contains **162** sequences (**76** positives / **86** negatives).

Features combine companion-style disorder / low-complexity / IDR charge blockiness with plant-motivated descriptors: Ser/Thr and SP/TP density in IDRs, aromatic sticker proxies (retained but not exclusive), composition-normalized charge residuals, PLAAC-related Q/N grammar, and MAP65-1-like ortholog-tail conservation. Production model weights use a **composition-proxy** disorder feature for CPU-light scoring; full-proteome metapredict mean fractions and AlphaFold pLDDT summaries enter the atlas as annotation tracks, with per-residue dual rails deposited under `data/tracks/` (v1.10.0). Optional experimental disorder retrains (metapredict / AF pLDDT) are deposited under `data/experimental/` and **do not** overwrite production weights. PlantLLPS fits an interpretable logistic / elastic-net model and a hybrid score that blends logistic output with the transparent companion composite (weights 0.55 / 0.25 / 0.20).

**Leave-one-group-out feature ablation** (v1.1.1): full-model mean logistic AUC ≈ **0.863** / hybrid ≈ **0.859**; largest drops from removing **map65_cterm** (Δ ≈ −0.036) then **disorder** (≈ −0.021), consistent with MAP65-grammar sensitivity.

### Fair benchmark suite versus ProX-style / PSAP / LLPhy / PLAAC

Preferred performance claims use the **leakage-controlled fair suite** (`scripts/eval_fair_benchmark.py`; leaderboard `docs/BENCHMARK_LEADERBOARD.md`, v1.3.0). Supervised methods retrain excluding every test ID; PSAP, LLPhyScore and PLAAC-approx are out-of-the-box; the ProX/PU-style baseline is a stdlib logistic on IDR/LCR/hydropathy/FCR/PLAAC/aromatic/κ features and is **not** the published Plant PS proteins LightGBM/PU model.

| Panel | PlantLLPS hybrid | ProX-style | LLPhy | PSAP | PLAAC-approx |
|---|---:|---:|---:|---:|---:|
| Gold leave-out (*n* = 28; train *n* = 134) | **0.959** | 0.467 | 0.769 | 0.718 | 0.497 |
| Family hold-out mean (13 families) | **0.833** | **0.836** | 0.700 | 0.689 | 0.399 |
| Ath→cereal (*n* = 9) | 1.000 | 0.833 | 0.111 | 0.633 | 0.278 |
| Ath→non_ath_crop (*n* = 38) | 0.962 | 0.855 | 0.278 | 0.573 | 0.449 |

**Benchmark honesty (slim).** Production-model gold hybrid ≈ **0.990** / logistic ≈ **0.979** remain **partially in-sample** and are listed only as a caveat—not as independent transfer. Ath→rice (*n* = 6) and Ath→cereal (*n* = 9) are exploratory. Family hold-out shares negatives by design (positives strictly held out). A MAP-screen reconstruction (*n* = 233) yields recomputed Spearman ρ(companion composite vs native PSAP) ≈ **0.202**, honestly distinct from the companion-cited ρ = **0.368** (exact 208-protein table unrecovered). Stratified five-fold logistic AUC on v13 averages ≈ **0.863** (out-of-fold ≈ 0.868).

### Explainable Ser-rich / charge-block / IDR grammar (SHAP)

For NAR+ interpretability (v1.4.0), the logistic head admits **exact linear SHAP** in logit space (`φ_i = w_i · z_i`), with probability-local scaling via `σ'(logit)`. Attributions aggregate into three grammar axes—**Ser-rich** (`ser_arg` + `map65_cterm`), **charge-block**, and **IDR/disorder** (+ low-complexity)—and join Gly-window occlusion (leave-one-window-out) to deposited IDR/domain intervals. Companion composite scores expose a transparent weighted-term decomposition (not SHAP). This is **not** KernelSHAP / TreeSHAP. Sidecar attributions for curated-gold ∪ Tier A (**1 159** proteins) live in `data/explainability/` (gzip; SQLite meta pointer only to respect the GitHub 100 MiB DB limit); protein pages and live `/score` render grammar bars.

**Anchor example.** VvMAP65-1 (F6HET4 / VIT_19s0090g01230) is **Tier B** (PlantLLPS ≈ 0.762; AF mean pLDDT 75.75; fraction pLDDT < 50 = 0.182; Ath ortholog Q9FLP0; soybean A0A0R0GGZ9 Tier B ≈ 0.805; tomato A0A3Q7HI29; native PSAP ≈ 0.0167; LLPhyScore percentile ≈ 29.8; Pfam PF03999). Deposited SHAP axes attribute the call primarily to **Ser-rich / MAP65 grammar** (logistic SHAP ≈ +0.133) and **charge blockiness** (≈ +0.121), not global IDR enrichment (IDR/disorder SHAP ≈ −0.083)—matching the companion Ser-rich tail biology that PSAP/LLPhyScore bury (Sahu et al. 2023; Guan, companion manuscript).

### GrapeCondensatome atlas content

The UniProt grapevine reference proteome UP000009183 was filtered and scored with PlantLLPS v13, yielding **25 972** proteins with hybrid score, percentile, companion-style features and confidence tiers:

| Tier | Definition | Count |
|---|---|---|
| **A** | Global high propensity (top PlantLLPS + high disorder/LC) | 1 142 |
| **B** | Ortholog-/grammar-rescued (moderate global score, strong local drivers; MAP65-1 class) | 2 728 |
| **C** | Context-dependent / mechanism-special or intermediate | 9 653 |
| **D** | Low-propensity baseline | 12 449 |

Joined annotation layers (through **v1.10.0**):

- **Expression / stress:** Expression Atlas **E-MTAB-8758** TPMs (25 810 / 25 972); multi-stress RNA-seq differentials **E-MTAB-4324** and **E-GEOD-67932** (both 25 810 matched).
- **Stress condensatome modules (v1.2.0):** hypergeometric enrichment of Tier A∪B among significant DE genes (18 contrast×direction tests; **1** FDR < 0.05 — E-GEOD-67932 stage-3 Botrytis down, fold ≈ **1.17**); **42** Tier A∪B response-pattern modules (size ≥ 5); **1 120** Tier A seed neighborhoods (Pearson *r* on stress log2FC); **1 198** deposited gene sets (`data/stress_modules/`; UI `/modules`). Descriptive overlap on short Atlas DE vectors—**not** WGCNA.
- **Atlas analyses (v1.5.0–v1.5.1):** Botrytis-down GO/keyword deep dive + STRING/stress seed hubs + VIT multi-accession score deltas (`data/stress_modules/deepdive/`, `data/string/stress_seed_hubs.tsv`, `data/isoforms/`); Pfam/InterPro Tier A∪B enrichment (Pfam FDR < 0.05: **144** / 424 tested); comparator disagreement (Tier A∪B ∩ low PSAP or low LLPhy: **1 555**; exploratory); structure–score joint flags (`af_frac_very_low` ≥ 0.25: **1 872**; threshold enrichments FDR-controlled). UI `/analyses`.
- **Disorder / structure summaries:** full-proteome metapredict means (25 972; mean fraction ≈ 0.256); AlphaFold mean pLDDT (*n* = 23 421) with proteome-wide low-pLDDT segments (22 449 nonempty; 2 551 atlas IDs lack AFDB).
- **Residue dual rails (v1.10.0):** per-residue metapredict disorder + AFDB pLDDT trajectories in `data/tracks/` (uint8 sidecars; protein-page sparklines + Jaccard cross-check). **Coverage honesty:** grape disorder **25 972 / 25 972**; grape pLDDT **3 956 / 25 972** (partial; AFDB-limited / resumable); eight-species conservative high-score candidates (grape + banana, date palm, kiwifruit, olive, peach, pear, pineapple) dual-rail **8 819 / 9 633** (disorder complete **9 633 / 9 633**). Full fruit-proteome residue deposits are out of scope.
- **Orthologs / cross-crop (through v1.5.2):** Ath/rice OrthoDB (12 706 / 14 546 matched) plus curated cross-crop hooks; **soybean** second-crop condensatome (UP000008827; 71 064 scored with grape-trained weights; 16 388 grape joins); **tomato** third-crop condensatome (UP000004994; **31 194** scored; tiers A/B/C/D = 1 387 / 3 424 / 11 566 / 14 817; 14 547 grape OrthoDB matches); conserved grape↔soy↔Ath Tier A∪B overlap (**271** triple Tier A/B; **1 645** grape-specific Tier A/B; deposits `data/cross_crop/`).
- **Fruit-crop scoring (v1.6.0) + cross-species high-score layer (v1.7.0):** **18** fruit UniProt proteomes scored with grape-trained v13 (tiers recalibrated per proteome — **not** species-retrained); eight unavailable taxa documented without invented deposits. Conservative clusters (Tier A∪B in ≥5 panel species): **1 438** rows; species-specific / orphan-high peaks deposited separately (`data/fruit_crops/cross_species/`; UI `/hooks#fruit-cross-species`). Orthology approx: OrthoDB → grape join → gene proxy → singleton.
- **Comparators:** full-proteome native PSAP and LLPhyScore (25 972 each) plus proteome **PLAAC-approx** (stdlib approx, not the Java PLAAC HMM).
- **Domains / IDR:** Pfam/InterPro coordinates (21 856 / 21 670) and IDR tracks (23 590 after full-AF expand).
- **Evidence overlays:** curated grape LLPS gold (*n* = 34; `/gold`); MAP-screen reconstruction deposit (`/map`); native *Vitis* STRING neighborhoods (taxid 29760; 25 867 mapped; 19 523 neighborhoods; `/network` — combined_score is descriptive adjacency, **not** LLPS evidence).
- **Orthogonal LLPS DB crosswalk (v1.8.0):** PhaSePro / LLPSDB / DrLLPS / CD-CODE ID joins against gold-panel union (*n* = **164**) and Tier A∪B∪C candidates; Tier A∪B **unique 3 396 / shared 474** (shared dominated by DrLLPS *Vitis* ortholog expansion). UI `/crosswalk`; methods `docs/LLPS_CROSSWALK.md`. Unique = catalog-relative under deposited match rules—not wet-lab novelty.
- **RNA-binding / nucleic-acid scope (v1.9.0):** RBP2GO 2.0 Ath RBPs joined via OrthoDB (**598** grape proteins); ATtRACT Ath motifs deposited but OrthoDB-sparse (**0** grape joins; plant-motif gap); native UniProt GO:0003723 (**817**) / KW RNA-binding (**613**); RG/RGG sequence grammar calls (**395**); any RBP-related annotation **1 943 / 25 972** (7.48%); Tier A∪B with RBP **348**. UI `/rbp`. Ortholog transfers and motif grammar are triage annotations—not grape CLIP proof (`data/rbp/METHOD.md`).
- **Explainability:** gold ∪ Tier A SHAP grammar sidecars (**1 159**; `/protein/:id` panel).

**Table 1. Resource contents at a glance (v1.10.0).**

| Layer | Content | Scale / note |
|---|---|---|
| PlantLLPS model | v13 logistic + hybrid; CLI `train` / `score` / `score-proteome` / `serve` | Train *n* = 162 |
| Grapevine atlas | UP000009183 scored + tiers A–D | 25 972 proteins |
| Expression / stress | E-MTAB-8758; E-MTAB-4324; E-GEOD-67932 | 25 810 matched |
| Stress modules | DE enrichment + pattern modules + GMT | 18 tests; 42 modules; 1 198 sets |
| Atlas analyses | Family / disagreement / structure + stress hubs | UI `/analyses` |
| Fair benchmarks | vs ProX-style / PSAP / LLPhy / PLAAC | Leave-gold-out hyb 0.959 |
| Explainability | Exact logistic SHAP + Gly occlusion | 1 159 gold ∪ Tier A |
| Residue dual rails | metapredict + AFDB pLDDT trajectories | Grape disorder 25 972; pLDDT 3 956; candidates dual-rail 8 819 / 9 633 |
| Structure summaries | metapredict means; AF pLDDT + low-pLDDT segments | 25 972 / 23 421 |
| Cross-crop | Ath/rice OrthoDB; soy + tomato; conserved AB | 71 064 soy; 31 194 tomato; 271 triple AB |
| Fruit crops | 18 scored proteomes + high-score clusters | Conservative 1 438; UI `/hooks` |
| LLPS DB crosswalk | PhaSePro / LLPSDB / DrLLPS / CD-CODE | Tier AB unique 3 396 / shared 474 |
| RBP / RG | RBP2GO Ath + UniProt GO/KW + RG/RGG | 1 943 annotated; UI `/rbp` |
| Comparators | native PSAP, LLPhyScore, PLAAC-approx | 25 972 each |
| Network / gold | STRING *Vitis*; curated LLPS gold | 34 gold proteins |
| Delivery | GitHub + Node/Python UI; public host | MIT; **v1.10.0** |

## IMPLEMENTATION AND USAGE

### Software package

PlantLLPS is distributed as an open-source Python package (stdlib-first CLI; Biopython/SciPy-compatible feature stack) with entry points:

```
plantllps train --corpus train_formal.csv
plantllps score --fasta query.faa --out scores.tsv
plantllps score-proteome --proteome UP000009183.fasta.gz
plantllps serve --db grape_condensatome.sqlite --port 8080
```

Optional `--use-metapredict` enables metapredict as an annotation track. Model weights (v13), metrics JSON, enrichment scripts, fair-benchmark, explainability, crosswalk and residue-track deposit scripts ship in the repository. Native PSAP/LLPhyScore reproduction requires a pinned Python ≤3.9 environment and `LLPHYSCORE_HOME` (documented in `docs/EXTERNAL_TOOLS.md`).

### Web interface and API

A lightweight web database (Node Express UI and Python `plantllps serve`) supports search by gene ID (`VIT_*`), UniProt accession, family name, tier and score threshold. Dedicated routes:

| Route | Purpose |
|---|---|
| `/` | Search / browse scored proteome |
| `/compare` | PlantLLPS vs native PSAP / LLPhyScore / PLAAC-approx |
| `/gold` | Curated experimentally anchored LLPS overlay |
| `/map` | MAP-screen reconstruction deposit |
| `/hooks` | Ath/rice + soybean/tomato + fruit-crop condensatomes; `#fruit-cross-species` high-score layer |
| `/network` | STRING *Vitis* neighborhoods |
| `/expression` | TPM and multi-stress RNA-seq contrasts |
| `/modules` | Stress DE × condensatome enrichment and response modules |
| `/analyses` | Family enrichment, comparator disagreement, structure–score |
| `/crosswalk` | Orthogonal LLPS DB unique vs shared candidates |
| `/rbp` | RNA-binding crosswalk (RBP2GO Ath / ATtRACT / UniProt / RG–RGG) |
| `/protein/:id` | Domain/IDR/AF summaries + dual-rail residue tracks + grammar SHAP (gold ∪ Tier A) + RBP/RG + external-DB badges |
| `/download` | Tier-/family-/keyword-/evidence-filtered TSV or FASTA (`evidence=unique|rbp|rg|…`) |

Bulk export: `GET /api/download?tier=&family=&keyword=&evidence=&format=tsv|fasta`. The public application is at **https://grape-condensatome.apps.grok.bot**; local `plantllps serve` / `npm start` and alternate hosting paths remain in `docs/HOSTING.md`. Fair leaderboard: `docs/BENCHMARK_LEADERBOARD.md`. Crosswalk methods: `docs/LLPS_CROSSWALK.md`.

### Typical workflows

1. **Prioritize grapevine condensate candidates.** Browse Tier A/B, filter by family (e.g. MAP65, WVD2/WDL, IQD), export via `/download`, and inspect expression + `/modules` stress enrichments and `/analyses` family/structure panels for wet-lab triage.
2. **Rescue MAP65-class drivers.** Use `/compare` and `/map` to contrast PlantLLPS placement of VvMAP65-1 against PSAP/LLPhyScore; open `/protein/F6HET4` for Ser-rich / charge-block SHAP bars and dual-rail residue tracks; follow Ath/soybean/tomato orthologs via `/hooks`.
3. **Cite fair benchmarks.** Prefer leave-gold-out / family / species-transfer AUCs from `docs/BENCHMARK_LEADERBOARD.md` over production in-sample gold AUCs.
4. **Score a custom FASTA.** Run `plantllps score` or `score-proteome` with released v13 weights; `/score` returns live feature SHAP for the logistic head.
5. **Cross-crop and fruit transfer.** Use deposited soybean and **tomato** condensatomes (grape-trained weights—not species-retrained), conserved grape↔soy↔Ath Tier A∪B tables, and the 18 fruit-crop deposits plus conservative high-score clusters on `/hooks`.
6. **Catalog uniqueness and nucleic-acid triage.** Filter `/crosswalk` or `/download?evidence=unique` for Tier B/C candidates absent from PhaSePro/LLPSDB/DrLLPS/CD-CODE under deposited rules; use `/rbp` for RBP2GO / UniProt / RG–RGG flags (triage only).

## DISCUSSION

### Specialty-crop uniqueness versus Plant PS/ProX and general LLPS databases

Plant PS proteins / ProX already owns the multi-species plant experimental-plus-prediction niche with a public site. LLPSDB, PhaSePro, PhaSepDB, DrLLPS and CD-CODE supply curated community standards. GrapeCondensatome’s differentiators are narrower but concrete: (i) **whole-proteome *Vitis* coverage**; (ii) a **plant/MAP-grammar-aware scorer** with leakage-controlled gains over PSAP/LLPhy/PLAAC on gold leave-out and explainable Ser-rich rescue of MAP65-1; (iii) an **integrated atlas** (expression, stress modules, atlas analyses, AF/metapredict dual rails, domains, STRING, native comparator columns, fruit-crop transfer, orthogonal catalog and RBP crosswalks) under an MIT-reproducible CLI. Orthogonal cross-annotation makes the niche measurable: of Tier A∪B proteins, **3 396 / 3 870** are absent from PhaSePro / LLPSDB / DrLLPS / CD-CODE under deposited ID/ortholog/gene match rules (methods: `docs/LLPS_CROSSWALK.md`)—catalog-relative uniqueness, not experimental proof. That combination is a specialty-crop resource, not a field-wide database update—and is accordingly framed for NAR Database Issue submission while remaining suitable for a Database / plant-resource venue if NAR scope or maturity thresholds are not met.

### Nucleic-acid / RBP scope (methods note)

Many plant condensate scaffolds are RNA-binding or RG/RGG-rich (stress granules, glycine-rich RBPs, nuclear RGG repeats). To place PlantLLPS candidates in that nucleic-acid narrative without inventing grape CLIP data, GrapeCondensatome (v1.9.0) deposits an RNA-binding crosswalk: (i) RBP2GO 2.0 *Arabidopsis* RBPs transferred by OrthoDB; (ii) ATtRACT Ath motif catalog (animal-dominated; Vitis absent); (iii) native UniProt GO:0003723 and Keyword RNA-binding on UP000009183; (iv) sequence RG/RGG grammar calls. The layer prioritizes RNA-related LLPS hypotheses for wet-lab triage and is documented in `data/rbp/METHOD.md`. It does **not** expand species coverage or replace experimental RBP atlases.

### Residue dual rails (methods note)

Per-residue metapredict disorder and AFDB pLDDT trajectories (v1.10.0) supply an orthogonal structure–disorder view beside PlantLLPS scores. Coverage is intentionally honest: grape disorder is complete; grape pLDDT remains partial (3 956 / 25 972 at tip; resumable); dual-rail completeness is highest for the eight-species conservative high-score candidate panel (8 819 / 9 633 with both rails). Missing AFDB models are not imputed. These tracks annotate hypotheses; they are not LLPS assays.

### Calibration and scope limits

Plant experimental positives remain sparse relative to human labels; OrthoDB joins and cross-crop hooks mitigate but do not eliminate scarcity. Cite fair-suite AUCs, not in-sample gold. Soybean, tomato and fruit-crop scores use grape-trained weights (tiers recalibrated per proteome where deposited). STRING combined_score mixes textmining and is not LLPS evidence. Stress modules are descriptive hypergeometric / short log2FC co-response patterns, not WGCNA. Logistic SHAP is exact for the linear head only. Comparator-disagreement flags are exploratory; family and structure threshold enrichments are FDR-controlled. Sequence scores are hypotheses for turbidity, microscopy and FRAP under validated MAP65-1 conditions (Sahu et al. 2023), not substitutes for experiment. RBP2GO/ATtRACT ortholog and RG/RGG flags are triage annotations, not grape RNA-binding proof. LLPS-DB “unique” flags are catalog-relative. Residue pLDDT coverage is incomplete for grape. Databases require versioning as LLPSDB and grapevine annotations advance.

### Outlook

A minted Zenodo DOI, enlarged family-held-out / multi-crop wet validation, and continuous curation would further strengthen NAR competitiveness beyond the live public host and v1.10.0 atlas / crosswalk / residue-track layer. Even without displacing Plant PS/ProX as the default plant LLPS portal, a stable, citeable grapevine condensatome fills a real specialty-crop gap for horticultural condensate biology.

## AVAILABILITY

PlantLLPS source code, trained model weights (v13), GrapeCondensatome TSV/SQLite tables, enrichment scripts, fair-benchmark deposits, explainability sidecars, crosswalk tables, residue-track sidecars and local web UIs are publicly available at:

**https://github.com/pingyguan-glitch/grape-condensatome**

Software version **v1.10.0** (MIT license; `main` tip `2818e2e319fe587b44c4803042bdcfcea05e273a`).

- **Public application URL:** https://grape-condensatome.apps.grok.bot
- **Zenodo DOI:** `[ZENODO_DOI — pending]` — *deposit PlantLLPS weights and GrapeCondensatome SQL/Parquet tables at journal submission via Zenodo–GitHub release mirroring; no DOI has been minted.*
- **Fair leaderboard:** `docs/BENCHMARK_LEADERBOARD.md`
- **LLPS DB crosswalk:** `docs/LLPS_CROSSWALK.md` · `data/llps_crosswalk/`
- **Explainability deposit:** `data/explainability/` (1 159 gold ∪ Tier A)
- **Stress modules:** `data/stress_modules/` (UI `/modules`)
- **Residue dual rails:** `data/tracks/` (UI `/protein/:id`)
- **RBP crosswalk:** `data/rbp/` (UI `/rbp`)

LLPSDB plant freeze inventory: Supplemental Table S1. Companion sequences and MAP-screen materials: Guan, companion manuscript (final journal/year/DOI pending). Key accessions: VvMAP65-1 (F6HET4), AtMAP65-1 (Q9FLP0), WVD2-like A0ACT1LGI9; soybean condensatome ortholog A0A0R0GGZ9; tomato condensatome ortholog A0A3Q7HI29.

## SUPPLEMENTARY DATA

**Supplemental Figure S1.** PlantLLPS v13 stratified cross-validation and gold-standard ROC panels (in-sample caveat marked). [PNG](./figures/supp/Figure_S1.png) · [PDF](./figures/supp/Figure_S1.pdf) · [legends](./figures/supp/LEGENDS.md).  
**Supplemental Figure S2.** Fair leave-gold-out / family hold-out / Ath→non_ath_crop transfer panels versus ProX-style, PSAP, LLPhyScore and PLAAC-approx (`docs/BENCHMARK_LEADERBOARD.md`). [PNG](./figures/supp/Figure_S2.png) · [PDF](./figures/supp/Figure_S2.pdf).  
**Supplemental Figure S3.** Feature-ablation (leave-one-group-out) and F6HET4 Ser-rich / charge-block / IDR SHAP grammar bars. [PNG](./figures/supp/Figure_S3.png) · [PDF](./figures/supp/Figure_S3.pdf).  
**Supplemental Figure S4.** Grapevine proteome-wide score distribution (*n* = 25 972) with Tier A/B and MAP anchors. [PNG](./figures/supp/Figure_S4.png) · [PDF](./figures/supp/Figure_S4.pdf).  
**Supplemental Figure S5.** Stress DE × condensatome enrichment and response-pattern modules (E-MTAB-4324 / E-GEOD-67932). [PNG](./figures/supp/Figure_S5.png) · [PDF](./figures/supp/Figure_S5.pdf).  
**Supplemental Table S1.** LLPSDB plant entries at corpus freeze (11 PS + 8 no-PS); see [`docs/supplemental-table-s1-llpsdb-plant.md`](./supplemental-table-s1-llpsdb-plant.md).  
**Supplemental Table S2.** Feature dictionary for PlantLLPS (including grammar-axis mapping for SHAP).  
**Supplemental Table S3.** Proteome-wide Tier A/B gene lists (v1.10.0 atlas; tiers unchanged since v1.0 scoring freeze).  
**Supplemental Table S4.** Stress-module enrichment summary (18 tests) and pattern-module sizes.  
**Supplemental Data S1.** Software package, scored proteome tables, fair-benchmark JSON/TSV, explainability sidecars, LLPS/RBP crosswalks, residue-track deposits and reproducibility scripts (GitHub **v1.10.0**). Regenerate figures: `python3 scripts/generate_supp_figures.py` → `docs/figures/supp/` (see `requirements-figures.txt`).

## ACKNOWLEDGEMENTS

The author thanks colleagues and the broader plant condensate community for curated resources such as LLPSDB. Acknowledgments will be updated upon submission.

## FUNDING

This work was supported by the Beijing Natural Science Foundation (6242019), National Natural Science Foundation of China (grant nos 32673384 and 32102308), the 2115 Talent Development Program of China Agricultural University, International Postdoctoral Exchange Fellowship Program (Talent-Introduction Program YJ20200335), and the China Postdoctoral Science Foundation (24329894G).

## CONFLICT OF INTEREST

The author declares no conflict of interest.

## REFERENCES

Banani SF, Lee HO, Hyman AA, Rosen MK. 2017. Biomolecular condensates: organizers of cellular biochemistry. Nat Rev Mol Cell Biol. 18:285–298.

Bolognesi B, Lorenzo Gotor N, Dhar R, Cirillo D, Baldrighi M, Tartaglia GG, Lehner B. 2016. A concentration-dependent liquid phase separation can cause toxicity upon increased protein expression. Cell Rep. 16:222–231.

Cai H, Vernon RM, Forman-Kay JD. 2022. An interpretable machine-learning algorithm to predict disordered protein phase separation based on biophysical interactions. Biomolecules. 12:1131.

Chu X, Sun T, Li Q, Xu Y, Zhang Z, Lai L, Pei J. 2022. Prediction of liquid–liquid phase separation proteins using machine learning. BMC Bioinformatics. 23:72.

Cock PJA, Antao T, Chang JT et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 25:1422–1423.

Das RK, Pappu RV. 2013. Conformations of intrinsically disordered proteins are influenced by linear sequence distributions of oppositely charged residues. Proc Natl Acad Sci USA. 110:13392–13397.

Dorone Y, Boeynaems S, Flores E et al. 2021. A prion-like protein regulator of seed germination undergoes hydration-dependent phase separation. Cell. 184:4284–4298.

Emenecker RJ, Griffith D, Holehouse AS. 2021. Metapredict: a fast, accurate, and easy-to-use predictor of consensus disorder and structure. Biophys J. 120:4312–4319.

Fang X, Wang L, Ishikawa R et al. 2019. Arabidopsis FLL2 promotes liquid–liquid phase separation of polyadenylation complexes. Nature. 569:265–269.

Guan P. Companion manuscript. Working title: A conserved Ser-rich tail marks the microtubule-bundling protein VvMAP65-1 as a phase-separation candidate. **Final journal, year, and DOI pending.**

Hardenberg M, Horvath A, Ambrus V, Fuxreiter M, Vendruscolo M. 2020. Widespread occurrence of the droplet state of proteins in the human proteome. Proc Natl Acad Sci USA. 117:33254–33262.

Jiang H, Wang S, Huang Y, He X, Cui H, Zhu X, Zheng Y. 2015. Phase transition of spindle-associated protein regulate spindle apparatus assembly. Cell. 163:108–122.

Jumper J, Evans R, Pritzel A et al. 2021. Highly accurate protein structure prediction with AlphaFold. Nature. 596:583–589.

Jung JH, Barbosa AD, Hutin S et al. 2020. A prion-like domain in ELF3 functions as a thermosensor in Arabidopsis. Nature. 585:256–260.

King MR, Petry S. 2020. Phase separation of TPX2 enhances and spatially coordinates microtubule nucleation. Nat Commun. 11:270.

Lancaster AK, Nutter-Upham A, Lindquist S, King OD. 2014. PLAAC: a web and command-line application to identify proteins with prion-like amino acid composition. Bioinformatics. 30:2501–2502.

Li Q, Peng X, Li Y et al. 2020. LLPSDB: a database of proteins undergoing liquid–liquid phase separation in vitro. Nucleic Acids Res. 48:D320–D327.

Liu X, Zhu JK, Zhao C. 2023. Liquid–liquid phase separation as a major mechanism of plant abiotic stress sensing and responses. Stress Biol. 3:56.

Mészáros B, Erdős G, Szabó B et al. 2020. PhaSePro: the database of proteins driving liquid–liquid phase separation. Nucleic Acids Res. 48:D360–D367.

van Mierlo G, Jansen JRG, Wang J, Poser I, van Heeringen SJ, Vermeulen M. 2021. Predicting protein condensate formation using machine learning. Cell Rep. 34:108705.

Ning W, Guo Y, Lin S et al. 2020. DrLLPS: a data resource of liquid–liquid phase separation in eukaryotes. Nucleic Acids Res. 48:D288–D295.

Papatheodorou I, Moreno P, Manning J et al. 2020. Expression Atlas update: from tissues to single cells. Nucleic Acids Res. 48:D77–D83.

Saar KL, Morgunov AS, Qi R et al. 2021. Learning the molecular grammar of protein condensates from sequence determinants and embeddings. Proc Natl Acad Sci USA. 118:e2019053118.

Sahu S, Chauhan P, Lumen E et al. 2023. Interplay of self-organization of microtubule asters and crosslinking protein condensates. PNAS Nexus. 2:pgad231.

Shin Y, Brangwynne GP. 2017. Liquid phase condensation in cell physiology and disease. Science. 357:eaaf4382.

Szklarczyk D, Kirsch R, Koutrouli M et al. 2023. The STRING database in 2023: protein–protein association networks and functional enrichment analyses for any sequenced genome of interest. Nucleic Acids Res. 51:D638–D646.

UniProt Consortium. 2023. UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Res. 51:D523–D531.

Varadi M, Anyango S, Deshpande M et al. 2022. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Res. 50:D439–D444.

Vernon RM, Chong PA, Tsang B et al. 2018. Pi-Pi contacts are an overlooked protein feature relevant to phase separation. eLife. 7:e31486.

Virtanen P, Gommers R, Oliphant TE et al. 2020. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods. 17:261–272.

Wang B, Zhang H, Huai J et al. 2022. Condensation of SEUSS promotes hyperosmotic stress tolerance in Arabidopsis. Nat Chem Biol. 18:1361–1369.

Wang X, Jiang B, Gu L et al. 2021. A photoregulatory mechanism of the circadian clock in Arabidopsis. Nat Plants. 7:1397–1408.

Zhu S, Gu J, Yao J et al. 2022. Liquid–liquid phase separation of RBGD2/4 is required for heat stress resistance in Arabidopsis. Dev Cell. 57:583–597.

Zhang H, et al. 2023. Plant PS proteins / ProX resource (β-isox plant condensate proteomics and portal). Mol. Plant. [complete citation to be confirmed against the published Plant PS proteins / ProX paper before submission].
