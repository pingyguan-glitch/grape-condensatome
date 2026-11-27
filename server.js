#!/usr/bin/env node
/**
 * GrapeCondensatome — Express web UI over bundled SQLite + PlantLLPS JS scorer.
 * Listens on 0.0.0.0:$PORT (default 8765).
 */

import express from "express";
import { createGunzip } from "node:zlib";
import { createReadStream, existsSync, readFileSync } from "node:fs";
import { createInterface } from "node:readline";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { DatabaseSync } from "node:sqlite";
import {
  openDb,
  getMeta,
  summarize,
  getProtein,
  searchProteins,
  listByTier,
  topProteins,
  filterProteins,
  getStringNeighborhood,
  listStringEnrichment,
  stringNetworkMeta,
  loadStringTopNeighbors,
  loadStressSidecar,
  loadStressModulesUi,
  getStressModulesUi,
  stressModulesMeta,
  listLlpsGold,
  DEFAULT_DB,
} from "./lib/db.js";
import { loadModel, scoreSequence } from "./lib/model.js";
import { explainFeatures } from "./lib/explain.js";

const __dirname = dirname(fileURLToPath(import.meta.url));
const PORT = Number(process.env.PORT || 8765);
const HOST = process.env.HOST || "0.0.0.0";
const PAGE_SIZE = 40;
const DB_PATH = process.env.GRAPE_DB || DEFAULT_DB;
const MODEL_PATH =
  process.env.GRAPE_MODEL || join(__dirname, "data", "model_formal.json");
const FASTA_PATH =
  process.env.GRAPE_FASTA || join(__dirname, "data", "UP000009183.fasta.gz");
const SOY_DB_PATH =
  process.env.SOY_DB ||
  join(__dirname, "data", "second_crop", "soybean_condensatome.sqlite");
const EXPLAIN_PATH =
  process.env.GRAPE_EXPLAIN ||
  join(__dirname, "data", "explainability", "attributions.json.gz");
const APP_VERSION = (() => {
  try {
    return readFileSync(join(__dirname, "VERSION"), "utf8").trim() || "0.0.0";
  } catch {
    return "0.0.0";
  }
})();

if (!existsSync(DB_PATH)) {
  console.error(`SQLite not found: ${DB_PATH}`);
  process.exit(1);
}

const db = openDb(DB_PATH);
const model = loadModel(MODEL_PATH);
/** @type {DatabaseSync | null} */
let soyDb = null;
if (existsSync(SOY_DB_PATH)) {
  soyDb = new DatabaseSync(SOY_DB_PATH, { readOnly: true });
  console.log(`Second-crop SQLite: ${SOY_DB_PATH}`);
} else {
  console.warn(`Soybean condensatome DB not found (${SOY_DB_PATH}); /hooks crop panel limited.`);
}
/** @type {Map<string, string>} */
const seqById = new Map();
/** @type {Map<string, object>} */
let explainById = new Map();
const app = express();
app.use(express.urlencoded({ extended: true, limit: "2mb" }));
app.use(express.static(join(__dirname, "public")));

async function loadExplainability(path = EXPLAIN_PATH) {
  if (!existsSync(path)) {
    explainById = new Map();
    console.warn(`Explainability sidecar missing: ${path}`);
    return explainById;
  }
  const gunzip = createGunzip();
  const chunks = [];
  await new Promise((resolve, reject) => {
    createReadStream(path)
      .pipe(gunzip)
      .on("data", (c) => chunks.push(c))
      .on("end", resolve)
      .on("error", reject);
  });
  const obj = JSON.parse(Buffer.concat(chunks).toString("utf8"));
  explainById = new Map(Object.entries(obj));
  console.log(`Explainability: ${explainById.size} proteins from ${path}`);
  return explainById;
}

async function loadFastaIndex(path) {
  if (!existsSync(path)) {
    console.warn(`FASTA not found (${path}); /api/download FASTA will be empty.`);
    return;
  }
  const stream = path.endsWith(".gz")
    ? createReadStream(path).pipe(createGunzip())
    : createReadStream(path);
  const rl = createInterface({ input: stream, crlfDelay: Infinity });
  let pid = null;
  let chunks = [];
  for await (const line of rl) {
    if (line.startsWith(">")) {
      if (pid) seqById.set(pid, chunks.join(""));
      const header = line.slice(1).split(/\s+/)[0];
      const parts = header.split("|");
      pid = parts.length > 1 ? parts[1] : header;
      chunks = [];
    } else {
      chunks.push(line.trim());
    }
  }
  if (pid) seqById.set(pid, chunks.join(""));
  console.log(`FASTA index: ${seqById.size} sequences from ${path}`);
}

function esc(s) {
  return String(s ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

function fmtScore(x) {
  const n = Number(x);
  return Number.isFinite(n) ? n.toFixed(4) : "—";
}

function fmtPct(x) {
  const n = Number(x);
  return Number.isFinite(n) ? n.toFixed(2) : "—";
}

function vitGeneName(geneName) {
  const g = String(geneName || "").trim();
  if (!g) return null;
  const m = g.match(/VIT_\d+s\d+g\d+/i);
  return m ? m[0] : null;
}

function genomeLinksHtml(p) {
  const vit = vitGeneName(p.gene_name);
  const links = [];
  if (vit) {
    const ens = `https://plants.ensembl.org/Vitis_vinifera/Gene/Summary?g=${encodeURIComponent(vit)}`;
    const loc = `https://plants.ensembl.org/Vitis_vinifera/Location/View?g=${encodeURIComponent(vit)}`;
    links.push(
      `<a href="${ens}" target="_blank" rel="noopener">Ensembl Plants · ${esc(vit)}</a>`
    );
    links.push(
      `<a href="${loc}" target="_blank" rel="noopener">Ensembl location / JBrowse</a>`
    );
  }
  links.push(
    `<a href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(p.protein_id)}/entry" target="_blank" rel="noopener">UniProt</a>`
  );
  if ((p.domain_pfam_ids || []).length) {
    const pf = p.domain_pfam_ids[0];
    links.push(
      `<a href="https://www.ebi.ac.uk/interpro/entry/pfam/${encodeURIComponent(pf)}/" target="_blank" rel="noopener">InterPro/Pfam ${esc(pf)}</a>`
    );
  }
  return links.join(" · ");
}

function trackHtml(length, segments, { title, tone }) {
  const L = Math.max(1, Number(length) || 1);
  const bars = (segments || [])
    .map((s) => {
      const start = Number(s.start) || 1;
      const end = Number(s.end) || start;
      const left = ((start - 1) / L) * 100;
      const width = Math.max(0.15, ((end - start + 1) / L) * 100);
      const tip = `${s.label || s.source || "seg"} ${start}–${end}`;
      return `<span class="track-seg ${tone}" style="left:${left}%;width:${width}%" title="${esc(tip)}"></span>`;
    })
    .join("");
  return `<div class="track">
    <div class="track-label">${esc(title)} <span class="muted">1–${esc(L)}</span></div>
    <div class="track-bar">${bars || `<span class="muted" style="font-size:0.8rem;padding:0.2rem">no segments</span>`}</div>
  </div>`;
}

function shapBarHtml(label, value, { maxAbs = 0.2, companion = null } = {}) {
  const v = Number(value) || 0;
  const m = Math.max(maxAbs, Math.abs(v), 1e-6);
  const pct = Math.min(50, (Math.abs(v) / m) * 50);
  const side = v >= 0 ? "pos" : "neg";
  const left = v >= 0 ? 50 : 50 - pct;
  const comp =
    companion != null && Number.isFinite(Number(companion))
      ? `<span class="muted"> · companion ${Number(companion).toFixed(3)}</span>`
      : "";
  return `<div class="shap-row">
    <div class="shap-label">${esc(label)}</div>
    <div class="shap-track"><span class="shap-zero"></span>
      <span class="shap-fill ${side}" style="left:${left}%;width:${pct}%"></span>
    </div>
    <div class="shap-val mono">${v >= 0 ? "+" : ""}${v.toFixed(4)}${comp}</div>
  </div>`;
}

function explainabilityHtml(expl, { length = 0, live = false } = {}) {
  if (!expl) {
    return `<div class="card" style="margin-top:1rem">
      <h3>Explainability (Ser-rich / charge-block / IDR)</h3>
      <p class="muted">No deposited SHAP attribution for this protein.
         Gold ∪ Tier A proteins are precomputed
         (<span class="mono">data/explainability/</span>); live paste-score shows
         feature SHAP without window occlusion.</p>
    </div>`;
  }
  const axes = expl.grammar_axes || expl.grammar?.axes || {};
  const axisOrder = ["ser_rich", "charge_block", "idr_disorder"];
  const axisLabels = {
    ser_rich: "Ser-rich / MAP65 grammar",
    charge_block: "Charge blockiness",
    idr_disorder: "IDR / disorder + LC",
  };
  const absVals = axisOrder.map((k) => {
    const a = axes[k] || {};
    const shap = a.logistic_shap_prob ?? a.hybrid ?? 0;
    return Math.abs(Number(shap) || 0);
  });
  const maxAbs = Math.max(0.08, ...absVals);
  const axisBars = axisOrder
    .map((k) => {
      const a = axes[k] || {};
      const shap = Number(a.logistic_shap_prob ?? a.hybrid ?? 0) || 0;
      return shapBarHtml(axisLabels[k] || k, shap, {
        maxAbs,
        companion: a.companion,
      });
    })
    .join("");
  const topFeats = (expl.top_features || [])
    .slice(0, 8)
    .map((f) => {
      const shap = Number(f.shap_prob_local || 0);
      return `<div><span class="k">${esc(f.feature)}</span>
        <span class="v">${shap >= 0 ? "+" : ""}${shap.toFixed(4)}
        <span class="muted">(x=${esc(
          Number.isFinite(Number(f.value)) ? Number(f.value).toFixed(4) : f.value
        )})</span></span></div>`;
    })
    .join("");
  const L = Math.max(1, Number(length) || Number(expl.length) || 1);
  const winSegs = (expl.grammar_windows || expl.occlusion?.windows || [])
    .filter((w) => Math.abs(Number(w.delta_score) || 0) > 1e-4)
    .map((w) => ({
      start: w.start,
      end: w.end,
      label: `${w.label || "win"} Δ=${Number(w.delta_score).toFixed(4)}`,
      source: w.label,
    }));
  const occTrack =
    winSegs.length > 0
      ? trackHtml(L, winSegs, { title: "Occlusion windows (Ser/charge tagged)", tone: "tone-shap" })
      : live
        ? `<p class="muted">Live score: feature SHAP only (no Gly-occlusion windows).</p>`
        : `<p class="muted">No high-|Δ| occlusion windows deposited.</p>`;
  const idrAttr = ((expl.segment_attributions || {}).idr || [])
    .slice(0, 6)
    .map(
      (s) =>
        `<div><span class="k">${esc(s.start)}–${esc(s.end)}</span>
         <span class="v">Δ̄ ${esc(Number(s.mean_window_delta).toFixed(4))}
         · ${esc(s.label || s.source || "idr")}</span></div>`
    )
    .join("");
  const method =
    expl.method_summary ||
    "Logistic: exact linear SHAP (φ=w·z). Companion terms shown beside bars. Window occlusion = Gly leave-one-window-out.";
  return `<div class="card" style="margin-top:1rem">
    <h3>Explainability (Ser-rich / charge-block / IDR)</h3>
    <p class="muted">${esc(method)}${
      live ? " · <strong>live</strong> paste-score" : " · deposited gold ∪ Tier A sidecar"
    }</p>
    <h4 style="margin-top:0.75rem">Grammar SHAP axes</h4>
    <div class="shap-panel">${axisBars}</div>
    <h4 style="margin-top:1rem">Top feature SHAP (probability-local)</h4>
    <div class="features">${topFeats || '<span class="muted">none</span>'}</div>
    ${occTrack}
    ${
      idrAttr
        ? `<h4 style="margin-top:1rem">IDR segment attributions</h4><div class="features">${idrAttr}</div>`
        : ""
    }
  </div>`;
}

function stressExprHtml(exprStress) {
  if (!exprStress || typeof exprStress !== "object") return "";
  const blocks = Object.entries(exprStress)
    .map(([acc, info]) => {
      const comps = info?.comparisons || {};
      const rows = Object.entries(comps)
        .map(([label, c]) => {
          const sig = c?.significant
            ? ' <span class="tier A" title="p&lt;0.05 and |log2FC|≥1">sig</span>'
            : "";
          return `<div><span class="k">${esc(label)}</span><span class="v">log2FC ${esc(
            c?.log2fc
          )} · p ${esc(c?.pvalue)}${sig}</span></div>`;
        })
        .join("");
      return `<h4 style="margin-top:0.75rem"><span class="mono">${esc(acc)}</span>
        · ${esc(info?.title || "")}</h4>
        <p class="muted">${esc(info?.stress || "")}
           · ${esc(info?.n_significant ?? 0)} significant contrast(s)</p>
        <div class="features">${rows || '<span class="muted">no contrasts</span>'}</div>`;
    })
    .join("");
  if (!blocks) return "";
  return `<div class="card" style="margin-top:1rem"><h3>Stress RNA-seq differentials</h3>
    <p class="muted">Expression Atlas analytics joined on <span class="mono">VIT_*</span> gene name.
       Call rule: nominal p&lt;0.05 and |log2FC|≥1.</p>
    ${blocks}</div>`;
}

function layout(title, body, query = "") {
  return `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>${esc(title)} · GrapeCondensatome</title>
<link rel="stylesheet" href="/style.css"/>
</head>
<body>
<div class="wrap">
  <header class="site-header">
    <div class="brand">
      <h1>GrapeCondensatome</h1>
      <p>Browse PlantLLPS scores for the grapevine reference proteome (UniProt UP000009183).</p>
    </div>
    <nav class="nav">
      <a href="/">Home</a>
      <a href="/tiers/A">Tier A</a>
      <a href="/tiers/B">Tier B</a>
      <a href="/tiers/C">Tier C</a>
      <a href="/tiers/D">Tier D</a>
      <a href="/compare">Compare</a>
      <a href="/gold">LLPS gold</a>
      <a href="/map">MAP screen</a>
      <a href="/hooks">Cross-crop</a>
      <a href="/network">STRING</a>
      <a href="/expression">Expression</a>
      <a href="/modules">Stress modules</a>
      <a href="/download">Download</a>
      <a href="/score">Score FASTA</a>
    </nav>
  </header>
  <div class="banner">
    PlantLLPS hybrid model (v13). Disorder: composition proxy by default; metapredict + AlphaFold pLDDT
    joined where computed. Curated LLPS gold badges + reconstructed MAP-screen ρ on /gold and /map.
    Orthologs (Ath/rice + soybean via OrthoDB), STRING Vitis neighborhoods (/network),
    stress DE × condensatome enrichment + modules (/modules), and side-by-side
    PSAP/LLPhyScore/PLAAC-approx on gold (+ proteome PLAAC).
    Second-crop soybean condensatome on /hooks. Treat ranks as hypotheses.
  </div>
  <form class="inline card" method="get" action="/search" style="margin-bottom:1.25rem">
    <label for="q" class="muted">Search UniProt ID, gene name, GO, or keyword</label>
    <input id="q" type="text" name="q" value="${esc(query)}" placeholder="e.g. F6HET4, MAP65, VIT_19s0090g01230"/>
    <button type="submit">Search</button>
  </form>
  ${body}
  <footer class="site-footer">
    PlantLLPS + GrapeCondensatome · ${esc(model.train_note?.slice(0, 160) || "")}…
  </footer>
</div>
</body>
</html>`;
}

function proteinTable(rows, empty = "No proteins found.") {
  if (!rows.length) return `<p class="muted">${esc(empty)}</p>`;
  const cells = rows
    .map((r) => {
      const pid = esc(r.protein_id);
      const tier = esc(r.tier || "?");
      const gene = esc(r.gene_name || "—");
      const pname = esc((r.protein_name || "").slice(0, 72) || "—");
      const tpm =
        r.expr_mean_tpm != null && Number.isFinite(Number(r.expr_mean_tpm))
          ? Number(r.expr_mean_tpm).toFixed(2)
          : "—";
      const ath = esc(r.ortholog_ath_gene || r.ortholog_ath || "—");
      const gma = esc(r.ortholog_gma_gene || r.ortholog_gma || "—");
      const gmaTier = r.ortholog_gma_tier
        ? `<span class="tier ${esc(r.ortholog_gma_tier)}">${esc(r.ortholog_gma_tier)}</span>`
        : "—";
      const gold =
        Number(r.llps_gold) > 0
          ? `<span class="badge gold" title="Curated LLPS evidence join">LLPS</span>`
          : "";
      return `<tr>
        <td class="mono"><a href="/protein/${pid}">${pid}</a> ${gold}</td>
        <td class="mono">${gene}</td>
        <td><span class="tier ${tier}">${tier}</span></td>
        <td>${fmtScore(r.score)}</td>
        <td>${fmtPct(r.percentile)}</td>
        <td>${fmtScore(r.map65_tail_like)}</td>
        <td>${fmtScore(r.frac_disorder_proxy)}</td>
        <td>${fmtScore(r.frac_disorder_metapredict)}</td>
        <td>${fmtScore(r.af_frac_very_low)}</td>
        <td>${fmtScore(r.plaac_nllr_approx)}</td>
        <td class="mono">${ath}</td>
        <td class="mono">${gma} ${gmaTier}</td>
        <td>${esc(tpm)}</td>
        <td class="muted">${pname}</td>
      </tr>`;
    })
    .join("");
  return `<div class="card" style="overflow-x:auto">
    <table>
      <thead><tr>
        <th>Protein</th><th>Gene</th><th>Tier</th><th>PlantLLPS</th><th>Percentile</th>
        <th>MAP65-like</th><th>Disorder proxy</th><th>metapredict</th><th>AF pLDDT&lt;50</th>
        <th>PLAAC≈</th><th>Ath ortholog</th><th>Soybean ortholog</th><th>Mean TPM</th><th>Protein name</th>
      </tr></thead>
      <tbody>${cells}</tbody>
    </table>
  </div>`;
}

function pager(base, page, total) {
  const nPages = Math.max(1, Math.ceil(total / PAGE_SIZE));
  const sep = base.includes("?") ? "&" : "?";
  const prev = page > 1 ? `${base}${sep}page=${page - 1}` : "#";
  const next = page < nPages ? `${base}${sep}page=${page + 1}` : "#";
  return `<div class="pager">
    <a class="btn secondary" href="${esc(prev)}">← Prev</a>
    <span class="muted">Page ${page} / ${nPages} · ${total} proteins</span>
    <a class="btn secondary" href="${esc(next)}">Next →</a>
  </div>`;
}

function secondCropMeta() {
  const meta = getMeta(db);
  let join = null;
  try {
    join = meta.second_crop_soybean ? JSON.parse(meta.second_crop_soybean) : null;
  } catch {
    join = null;
  }
  let soy = null;
  if (soyDb) {
    try {
      const n = soyDb.prepare("SELECT COUNT(*) AS n FROM proteins").get().n;
      const tiers = {};
      for (const r of soyDb.prepare("SELECT tier, COUNT(*) AS c FROM proteins GROUP BY tier").all()) {
        tiers[r.tier] = r.c;
      }
      soy = { path: SOY_DB_PATH, n_proteins: n, tiers };
    } catch {
      soy = { path: SOY_DB_PATH, error: "unreadable" };
    }
  }
  return { join, soy };
}

app.get("/api/health", (_req, res) => {
  const sn = stringNetworkMeta(db);
  const sm = stressModulesMeta(db);
  let explainMeta = null;
  let atlasAnalyses = null;
  try {
    const raw = getMeta(db).explainability;
    explainMeta = raw ? JSON.parse(raw) : null;
  } catch {
    explainMeta = null;
  }
  try {
    const raw = getMeta(db).atlas_analyses_stress_network_isoforms;
    atlasAnalyses = raw ? JSON.parse(raw) : null;
  } catch {
    atlasAnalyses = null;
  }
  const deep = atlasAnalyses?.analyses?.["2_stress_module_deepdive"];
  const hubs = atlasAnalyses?.analyses?.["6_string_stress_seed_hubs"];
  const iso = atlasAnalyses?.analyses?.["7_isoform_multi_protein_genes"];
  res.json({
    ok: true,
    version: APP_VERSION,
    db: DB_PATH,
    n: summarize(db).n_proteins,
    second_crop: secondCropMeta(),
    string: sn.version
      ? {
          version: sn.version,
          taxon: sn.taxon,
          min_score: sn.min_score,
          mapped: sn.mapped,
          seed_edges: sn.seed_edges,
        }
      : null,
    stress_modules: sm
      ? {
          n_modules: sm.n_modules,
          n_enrichment_tests: sm.n_enrichment_tests,
          n_enrichment_fdr_lt_0_05: sm["n_enrichment_fdr_lt_0.05"],
          n_gene_sets: sm.n_gene_sets,
          accessions: sm.accessions,
        }
      : null,
    atlas_analyses: atlasAnalyses
      ? {
          version: atlasAnalyses.version,
          deepdive_primary_fdr_terms: deep?.n_primary_term_fdr_lt_0_05 ?? deep?.["n_primary_term_fdr_lt_0.05"],
          hubs_top_n: hubs?.top_n_deposited,
          isoform_multi_vit: iso?.n_multi_accession_vit,
        }
      : null,
    explainability: {
      loaded: explainById.size,
      meta: explainMeta,
      sidecar: EXPLAIN_PATH,
    },
  });
});

app.get("/api/search", (req, res) => {
  const q = String(req.query.q || "");
  const { rows, total } = searchProteins(db, q, { limit: 20, offset: 0 });
  res.json({
    query: q,
    total,
    results: rows.map((r) => ({
      protein_id: r.protein_id,
      gene_name: r.gene_name || "",
      protein_name: r.protein_name || "",
      score: r.score,
      tier: r.tier,
      percentile: r.percentile,
      map65_tail_like: r.map65_tail_like,
      go_terms: r.go_terms || [],
      keywords: r.keywords || [],
    })),
  });
});

app.get("/", (_req, res) => {
  const summary = summarize(db);
  const meta = getMeta(db);
  const top = topProteins(db, 15);
  const tiers = summary.tiers || {};
  const stats = ["A", "B", "C", "D"]
    .map(
      (t) =>
        `<div class="card stat"><div class="n">${esc(tiers[t] || 0)}</div>
         <div class="l"><a href="/tiers/${t}">Tier ${t}</a></div></div>`
    )
    .join("");
  let exprMeta = "";
  try {
    const ea = meta.expression_atlas ? JSON.parse(meta.expression_atlas) : null;
    if (ea) {
      exprMeta = ` · Expression Atlas <span class="mono">${esc(ea.experiment || "")}</span>
        matched ${esc(ea.n_proteins_matched || "—")} / ${esc(summary.n_proteins)}`;
    }
    const sr = meta.stress_rnaseq ? JSON.parse(meta.stress_rnaseq) : null;
    if (sr?.accessions?.length) {
      exprMeta += ` · stress RNA-seq <span class="mono">${esc(sr.accessions.join(", "))}</span>
        matched ${esc(sr.n_proteins_with_any_stress || "—")} / ${esc(summary.n_proteins)}`;
    }
    const af = meta.alphafold_plddt ? JSON.parse(meta.alphafold_plddt) : null;
    if (af?.n_with_segments_field != null || af?.n_with_af_n_disorder != null) {
      exprMeta += ` · AF segments ${esc(
        af.n_with_segments_field ?? af.n_with_af_n_disorder
      )} / ${esc(af.n_with_mean || "—")} models`;
    }
  } catch {
    /* ignore */
  }
  const body = `
    <div class="grid stats">
      <div class="card stat"><div class="n">${summary.n_proteins}</div><div class="l">Proteins in DB</div></div>
      ${stats}
    </div>
    <div class="card" style="margin-bottom:1rem">
      <h2>Atlas</h2>
      <p>Status: <strong>${esc(meta.status || "unknown")}</strong>
         · ${esc(meta.proteome || "—")}
         · n_scored: ${esc(meta.n_scored || "—")}${exprMeta}</p>
      <p class="muted">${esc(meta.disclaimer || "")}</p>
    </div>
    <h2 style="font-family:var(--font-display)">Top scores</h2>
    ${proteinTable(top)}`;
  res.type("html").send(layout("Home", body));
});

app.get("/expression", (_req, res) => {
  const meta = getMeta(db);
  const summary = summarize(db);
  let ea = null;
  let sr = null;
  try {
    ea = meta.expression_atlas ? JSON.parse(meta.expression_atlas) : null;
  } catch {
    ea = null;
  }
  try {
    sr = meta.stress_rnaseq ? JSON.parse(meta.stress_rnaseq) : null;
  } catch {
    sr = null;
  }
  const baseline = ea
    ? `<div class="card"><h3>Baseline TPM · <span class="mono">${esc(
        ea.experiment || ""
      )}</span></h3>
       <p>Matched <strong>${esc(ea.n_proteins_matched)}</strong> / ${esc(
         summary.n_proteins
       )}
       · Atlas genes ${esc(ea.n_genes_in_atlas || "—")}</p>
       <p class="muted mono">${esc(ea.source || "")}</p></div>`
    : `<p class="muted">No baseline Expression Atlas join in meta.</p>`;
  const stressCards = sr?.experiments
    ? Object.entries(sr.experiments)
        .map(
          ([acc, info]) => `<div class="card" style="margin-top:1rem">
        <h3><span class="mono">${esc(acc)}</span> · ${esc(info.title || "")}</h3>
        <p>${esc(info.stress || "")}</p>
        <p>Technology: ${esc(info.technology || "—")}
           · matched <strong>${esc(info.n_proteins_matched)}</strong> / ${esc(
             sr.n_proteins
           )}
           (${esc(info.match_rate)}%)</p>
        <p class="muted">${esc(info.citation || "")}</p>
        <p class="muted">Contrasts: ${esc(
          Object.values(info.contrasts || {}).join(" · ")
        )}</p>
        <p class="muted">${esc(info.call_rule || "")}</p>
        <p class="muted mono">${esc(info.source || "")}</p>
      </div>`
        )
        .join("")
    : "";
  const body = `
    <div class="card"><h2>Expression joins</h2>
      <p class="muted">Baseline berry TPM plus multi-stress RNA-seq differentials from
      EMBL-EBI Expression Atlas (public experiment files). Gene join key:
      UniProt <span class="mono">gene_name</span> = Atlas <span class="mono">VIT_*</span>.</p>
      ${
        sr?.choice_rationale
          ? `<p class="muted">${esc(sr.choice_rationale)}</p>`
          : ""
      }
    </div>
    ${baseline}
    ${stressCards || '<p class="muted">No stress RNA-seq join yet.</p>'}
    <div class="card" style="margin-top:1rem">
      <h3>Condensatome × stress enrichment</h3>
      <p class="muted">DE ∩ Tier A∪B hypergeometric tests and response-pattern modules
         for high-tier proteins are on
         <a href="/modules">/modules</a>
         (gene sets under <span class="mono">data/stress_modules/</span>).</p>
    </div>`;
  res.type("html").send(layout("Expression", body));
});

app.get("/modules", (req, res) => {
  const sm = stressModulesMeta(db);
  const ui = getStressModulesUi();
  const tab = String(req.query.tab || "enrichment") === "modules" ? "modules" : "enrichment";
  const page = Math.max(1, Number(req.query.page) || 1);
  const params = ui?.parameters || sm?.parameters || {};
  const honesty = sm?.honesty || ui?.honesty || "";

  if (!ui && !sm) {
    const body = `
      <div class="card">
        <h2>Stress condensatome modules</h2>
        <p class="muted">Not deposited yet. Run
          <span class="mono">python3 scripts/enrich_stress_condensatome_modules.py</span>.</p>
      </div>`;
    return res.type("html").send(layout("Stress modules", body));
  }

  const enrich = ui?.enrichment || sm?.top_enrichment || [];
  const modules = ui?.modules || [];
  const rows = tab === "modules" ? modules : enrich;
  const total = rows.length;
  const pages = Math.max(1, Math.ceil(total / PAGE_SIZE));
  const slice = rows.slice((page - 1) * PAGE_SIZE, page * PAGE_SIZE);

  let tableHtml = "";
  if (tab === "enrichment") {
    const trs = slice
      .map((r) => {
        const fold =
          r.enrichment_fold != null && Number.isFinite(Number(r.enrichment_fold))
            ? Number(r.enrichment_fold).toFixed(3)
            : "—";
        const fdr =
          r.bh_fdr != null && Number.isFinite(Number(r.bh_fdr))
            ? Number(r.bh_fdr) < 0.001
              ? Number(r.bh_fdr).toExponential(2)
              : Number(r.bh_fdr).toFixed(4)
            : "—";
        return `<tr>
          <td class="mono">${esc(r.accession)}</td>
          <td class="mono">${esc(r.contrast_id)}</td>
          <td>${esc(r.direction)}</td>
          <td>${esc(r.n_de)}</td>
          <td>${esc(r.n_de_tier_ab)}</td>
          <td>${fold}</td>
          <td>${fmtP(r.hypergeom_pvalue)}</td>
          <td>${esc(fdr)}</td>
          <td class="mono">${esc(r.gene_set_id || "")}</td>
        </tr>`;
      })
      .join("");
    tableHtml = `<div class="table-wrap"><table>
      <thead><tr>
        <th>Accession</th><th>Contrast</th><th>Dir</th>
        <th>n DE</th><th>∩ A∪B</th><th>Fold</th><th>p</th><th>FDR</th><th>Gene set</th>
      </tr></thead>
      <tbody>${trs || '<tr><td colspan="9" class="muted">No enrichment rows.</td></tr>'}</tbody>
    </table></div>`;
  } else {
    const trs = slice
      .map((r) => {
        const hub = r.hub_protein_id
          ? `<a class="mono" href="/protein/${esc(r.hub_protein_id)}">${esc(
              r.hub_protein_id
            )}</a>`
          : "—";
        const tops = String(r.top_members || "")
          .split(",")
          .filter(Boolean)
          .slice(0, 5)
          .map(
            (id) =>
              `<a class="mono" href="/protein/${esc(id.trim())}">${esc(id.trim())}</a>`
          )
          .join(" ");
        return `<tr>
          <td class="mono">${esc(r.module_id)}</td>
          <td>${esc(r.n_members)}</td>
          <td>${esc(r.n_tier_a)} / ${esc(r.n_tier_b)}</td>
          <td class="mono">${esc(r.pattern || "")}</td>
          <td>${esc(r.pattern_label || "")}</td>
          <td>${hub}</td>
          <td>${tops}</td>
        </tr>`;
      })
      .join("");
    tableHtml = `<div class="table-wrap"><table>
      <thead><tr>
        <th>Module</th><th>n</th><th>A / B</th><th>Pattern</th>
        <th>Label</th><th>Hub</th><th>Top members</th>
      </tr></thead>
      <tbody>${trs || '<tr><td colspan="7" class="muted">No modules.</td></tr>'}</tbody>
    </table></div>`;
  }

  const body = `
    <div class="card">
      <h2>Stress condensatome modules</h2>
      <p>Condensatome × stress DE enrichment for
         <span class="mono">${esc((params.accessions || sm?.accessions || []).join(", "))}</span>
         and response-pattern modules among Tier A∪B proteins.</p>
      <p class="muted">Background stress-matched <strong>${esc(
        params.n_background || "—"
      )}</strong>
         · Tier A∪B rate ≈ <strong>${esc(params.background_tier_ab_rate ?? "—")}</strong>
         · enrichment tests <strong>${esc(sm?.n_enrichment_tests ?? enrich.length)}</strong>
         (FDR&lt;0.05: <strong>${esc(sm?.["n_enrichment_fdr_lt_0.05"] ?? "—")}</strong>)
         · pattern modules <strong>${esc(sm?.n_modules ?? modules.length)}</strong>
         · gene sets <strong>${esc(sm?.n_gene_sets ?? "—")}</strong>
         · seed neighborhoods <strong>${esc(sm?.n_seed_neighborhoods ?? "—")}</strong>.</p>
      <p class="muted">${esc(honesty)}</p>
      <p>Deposits:
        <span class="mono">data/stress_modules/de_enrichment.tsv</span> ·
        <span class="mono">coexpression_modules.tsv</span> ·
        <span class="mono">gene_sets.gmt.gz</span> ·
        <span class="mono">module_members.tsv.gz</span> ·
        <span class="mono">deepdive/</span> (v1.5 GO/keyword + heatmap) ·
        <span class="mono">data/string/stress_seed_hubs.tsv</span> ·
        <span class="mono">data/isoforms/</span>
      </p>
      <p class="muted">Deep-dive / hubs / isoform deltas:
         <span class="mono">python3 scripts/analyze_stress_network_isoforms.py</span>
         (primary = FDR Botrytis stage-3 down ∩ Tier A∪B; other contrasts exploratory).</p>
      <p>
        Table:
        <a href="/modules?tab=enrichment">DE enrichment</a> ·
        <a href="/modules?tab=modules">Pattern modules</a>
      </p>
      <p class="muted">${esc(total)} rows (page ${esc(page)} / ${esc(pages)}).</p>
      ${tableHtml}
      ${pagerHtml("/modules", page, pages, { tab })}
    </div>`;
  res.type("html").send(layout("Stress modules", body));
});

app.get("/search", (req, res) => {
  const q = String(req.query.q || "");
  const page = Math.max(1, Number(req.query.page || 1));
  const { rows, total } = searchProteins(db, q, {
    limit: PAGE_SIZE,
    offset: (page - 1) * PAGE_SIZE,
  });
  const base = `/search?q=${encodeURIComponent(q)}`;
  const body = `
    <div class="card"><h2>Search results</h2>
    <p class="muted">Query: <span class="mono">${esc(q)}</span> · ${total} hit(s)</p></div>
    ${proteinTable(rows)}
    ${total > PAGE_SIZE ? pager(base, page, total) : ""}`;
  res.type("html").send(layout(`Search: ${q}`, body, q));
});

app.get("/tiers/:tier", (req, res) => {
  const tier = String(req.params.tier || "").toUpperCase();
  const page = Math.max(1, Number(req.query.page || 1));
  const { rows, total } = listByTier(db, tier, {
    limit: PAGE_SIZE,
    offset: (page - 1) * PAGE_SIZE,
  });
  const body = `
    <div class="card">
      <h2>Tier <span class="tier ${esc(tier)}">${esc(tier)}</span></h2>
      <p class="muted">${total} proteins · sorted by PlantLLPS score</p>
    </div>
    ${proteinTable(rows)}
    ${pager(`/tiers/${tier}`, page, total)}`;
  res.type("html").send(layout(`Tier ${tier}`, body));
});

app.get("/protein/:id", (req, res) => {
  const p = getProtein(db, req.params.id);
  if (!p) {
    return res
      .status(404)
      .type("html")
      .send(
        layout(
          "Missing",
          `<div class="card"><p>Protein <span class="mono">${esc(req.params.id)}</span> not in DB.</p></div>`
        )
      );
  }
  const feats = p.features || {};
  const featHtml = Object.keys(feats)
    .sort()
    .map(
      (k) =>
        `<div><span class="k">${esc(k)}</span><span class="v">${esc(feats[k])}</span></div>`
    )
    .join("");
  const go = (p.go_terms || []).slice(0, 24).map((g) => esc(g)).join(", ") || "—";
  const kw = (p.keywords || []).slice(0, 24).map((g) => esc(g)).join(", ") || "—";
  const pfamList =
    (p.domain_pfam_ids || []).map((x) => esc(x)).join(", ") || "—";
  const iprList =
    (p.domain_interpro_ids || []).slice(0, 12).map((x) => esc(x)).join(", ") || "—";
  const domainTrack = trackHtml(p.length, p.domains, {
    title: "Domains",
    tone: "tone-domain",
  });
  const idrTrack = trackHtml(p.length, p.idr_segments, {
    title: "IDR",
    tone: "tone-idr",
  });
  const domainRows = (p.domains || [])
    .slice(0, 40)
    .map(
      (d) =>
        `<div><span class="k">${esc(d.start)}–${esc(d.end)}</span><span class="v">${esc(
          d.label || d.accession || d.source || ""
        )}${d.accession ? ` · <span class="mono">${esc(d.accession)}</span>` : ""} · ${esc(
          d.source || ""
        )}</span></div>`
    )
    .join("");
  const idrRows = (p.idr_segments || [])
    .slice(0, 40)
    .map(
      (d) =>
        `<div><span class="k">${esc(d.start)}–${esc(d.end)}</span><span class="v">${esc(
          d.source || ""
        )}${d.label ? ` · ${esc(d.label)}` : ""}</span></div>`
    )
    .join("");
  const body = `
    <div class="detail-hero">
      <div class="card">
        <h2 class="mono">${esc(p.protein_id)}
          ${
            Number(p.llps_gold) > 0
              ? `<span class="badge gold">LLPS evidence</span>`
              : ""
          }</h2>
        <p><strong>${esc(p.gene_name || "—")}</strong>
           · <span class="tier ${esc(p.tier)}">${esc(p.tier)}</span>
           · length ${esc(p.length)}</p>
        <p>${esc(p.protein_name || "No UniProt protein name in DB.")}</p>
        <p><span class="muted">Links:</span> ${genomeLinksHtml(p)}</p>
        <p><span class="muted">GO:</span> ${go}</p>
        <p><span class="muted">Keywords:</span> ${kw}</p>
        <p><span class="muted">Pfam:</span> ${pfamList}
           · <span class="muted">InterPro:</span> ${iprList}</p>
        <p class="muted mono">hash=${esc(p.sequence_hash)} · source=${esc(p.source)}</p>
      </div>
      <div class="card">
        <h3>Scores</h3>
        <p>PlantLLPS: <strong>${fmtScore(p.score)}</strong></p>
        <p>Companion: ${fmtScore(p.score_companion)}</p>
        <p>Logistic: ${fmtScore(p.score_logistic)}</p>
        <p>Percentile: ${fmtPct(p.percentile)}</p>
        <p>MAP65-tail-like: ${fmtScore(p.map65_tail_like)}</p>
        <p>Disorder proxy: ${fmtScore(p.frac_disorder_proxy)}</p>
        <p>metapredict frac: ${fmtScore(p.frac_disorder_metapredict)}</p>
        <p>AF mean pLDDT: ${fmtScore(p.af_mean_plddt)}
           · frac&lt;50: ${fmtScore(p.af_frac_very_low)}</p>
        <p>PLAAC-approx NLLR: ${fmtScore(p.plaac_nllr_approx)}</p>
        <p>PSAP: ${fmtScore(p.score_psap)}
           · LLPhyScore %ile: ${fmtScore(p.score_llphyscore_percentile)}
           ${p.comparator_set ? `(<span class="mono">${esc(p.comparator_set)}</span>)` : "<span class='muted'>(gold only)</span>"}</p>
        <p>Expr. mean TPM: ${
          p.expr_mean_tpm != null && Number.isFinite(Number(p.expr_mean_tpm))
            ? Number(p.expr_mean_tpm).toFixed(2)
            : "—"
        }
           ${p.expr_experiment ? `(<span class="mono">${esc(p.expr_experiment)}</span>)` : ""}</p>
        <p>Expr. max TPM: ${
          p.expr_max_tpm != null && Number.isFinite(Number(p.expr_max_tpm))
            ? Number(p.expr_max_tpm).toFixed(2)
            : "—"
        }</p>
      </div>
    </div>
    ${
      p.llps_gold_entries && p.llps_gold_entries.length
        ? `<div class="card" style="margin-top:1rem">
            <h3>Curated LLPS evidence</h3>
            <p class="muted">Joined from <span class="mono">data/curated/llps_experimental_gold.tsv</span>
               (companion gold + modest literature / LLPSDB rows). Not an exhaustive catalog.</p>
            <div class="features">${p.llps_gold_entries
              .map(
                (e) =>
                  `<div><span class="k">${esc(e.gene_symbol || e.protein_id)}</span>
                   <span class="v">${esc(e.label)} · ${esc(e.condensate_type)} · ${esc(
                     e.citation
                   )}</span></div>`
              )
              .join("")}</div></div>`
        : ""
    }
    <div class="card" style="margin-top:1rem">
      <h3>Orthologs / cross-crop</h3>
      <p>Family: <strong>${esc(p.cross_crop_family || "—")}</strong>
         · OrthoDB: <span class="mono">${esc(p.orthodb_id || "—")}</span></p>
      <p>Arabidopsis: <span class="mono">${esc(p.ortholog_ath || "—")}</span>
         ${p.ortholog_ath_gene ? `(${esc(p.ortholog_ath_gene)})` : ""}
         ${
           p.ortholog_ath
             ? `· <a href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(p.ortholog_ath)}/entry" target="_blank" rel="noopener">UniProt</a>`
             : ""
         }</p>
      <p>Rice: <span class="mono">${esc(p.ortholog_osa || "—")}</span>
         ${p.ortholog_osa_gene ? `(${esc(p.ortholog_osa_gene)})` : ""}
         ${
           p.ortholog_osa
             ? `· <a href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(p.ortholog_osa)}/entry" target="_blank" rel="noopener">UniProt</a>`
             : ""
         }</p>
      <p>Soybean (<em>Glycine max</em>): <span class="mono">${esc(p.ortholog_gma || "—")}</span>
         ${p.ortholog_gma_gene ? `(${esc(p.ortholog_gma_gene)})` : ""}
         ${
           p.ortholog_gma_tier
             ? `· <span class="tier ${esc(p.ortholog_gma_tier)}">${esc(p.ortholog_gma_tier)}</span>`
             : ""
         }
         ${
           p.ortholog_gma_score != null && Number.isFinite(Number(p.ortholog_gma_score))
             ? `· score ${fmtScore(p.ortholog_gma_score)}`
             : ""
         }
         ${
           p.ortholog_gma
             ? `· <a href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(p.ortholog_gma)}/entry" target="_blank" rel="noopener">UniProt</a>`
             : ""
         }</p>
    </div>
    <div class="card" style="margin-top:1rem">
      <h3>Domain / IDR track</h3>
      <p class="muted">Domains: UniProt FT Domain + InterPro Pfam coords where FT empty.
         IDR: AlphaFold pLDDT&lt;50 (when deposited) ∪ UniProt REGION Disordered (MobiDB-lite).
         metapredict supplies mean frac only — no per-residue IDR track.</p>
      ${domainTrack}
      ${idrTrack}
      ${
        domainRows
          ? `<h4 style="margin-top:1rem">Domain intervals</h4><div class="features">${domainRows}</div>`
          : `<p class="muted">No domain coordinates for this protein.</p>`
      }
      ${
        idrRows
          ? `<h4 style="margin-top:1rem">IDR intervals (${esc(
              p.idr_n_residues ?? ""
            )} residues · ${esc(p.idr_sources || "")})</h4><div class="features">${idrRows}</div>`
          : `<p class="muted">No IDR segments deposited for this protein.</p>`
      }
    </div>
    ${
      p.af_disorder_segments && p.af_disorder_segments.length
        ? `<div class="card" style="margin-top:1rem"><h3>AlphaFold low-pLDDT segments (&lt;50)</h3>
           <p class="muted">${esc(p.af_n_disorder_residues ?? "")} residues in
           ${p.af_disorder_segments.length} segment(s)
           ${p.af_model_id ? `· <span class="mono">${esc(p.af_model_id)}</span>` : ""}</p>
           <div class="features">${p.af_disorder_segments
             .slice(0, 40)
             .map(
               (s) =>
                 `<div><span class="k">${esc(s.start)}–${esc(s.end)}</span><span class="v">len ${esc(s.length)}</span></div>`
             )
             .join("")}</div></div>`
        : ""
    }
    ${
      p.expr_tpm
        ? `<div class="card" style="margin-top:1rem"><h3>Expression Atlas TPM by condition</h3>
           <p class="muted">Baseline <span class="mono">${esc(
             p.expr_experiment || "E-MTAB-8758"
           )}</span></p>
           <div class="features">${Object.entries(p.expr_tpm)
             .map(
               ([k, v]) =>
                 `<div><span class="k">${esc(k)}</span><span class="v">${esc(v)}</span></div>`
             )
             .join("")}</div></div>`
        : ""
    }
    ${stressExprHtml(p.expr_stress)}
    ${stringNeighborhoodHtml(p.protein_id)}
    ${explainabilityHtml(explainById.get(p.protein_id), { length: p.length })}
    <div class="card" style="margin-top:1rem">
      <h3>Features</h3>
      <div class="features">${featHtml}</div>
    </div>`;
  res.type("html").send(layout(p.protein_id, body));
});

function fmtP(x) {
  const n = Number(x);
  if (!Number.isFinite(n)) return "—";
  if (n < 1e-4) return n.toExponential(2);
  return n.toFixed(4);
}

function stringNeighborhoodHtml(proteinId) {
  const nb = getStringNeighborhood(db, proteinId);
  if (!nb) {
    return `<div class="card" style="margin-top:1rem">
      <h3>STRING neighborhood</h3>
      <p class="muted">No STRING neighborhood deposited for this protein (unmapped or no edges
         ≥ medium confidence). See <a href="/network">/network</a>.</p>
    </div>`;
  }
  const top = nb.top_neighbors || [];
  const rows = top
    .slice(0, 15)
    .map((t) => {
      const pid = t.protein_id;
      return `<tr>
        <td class="mono"><a href="/protein/${esc(pid)}">${esc(pid)}</a></td>
        <td>${esc(t.gene_name || "—")}</td>
        <td><span class="tier ${esc(t.tier || "")}">${esc(t.tier || "—")}</span></td>
        <td>${fmtScore(t.score)}</td>
        <td>${esc(t.combined_score)}</td>
      </tr>`;
    })
    .join("");
  const fold =
    nb.enrichment_fold != null && Number.isFinite(Number(nb.enrichment_fold))
      ? Number(nb.enrichment_fold).toFixed(2)
      : "—";
  return `<div class="card" style="margin-top:1rem">
    <h3>STRING neighborhood</h3>
    <p class="muted">Vitis STRING (taxid 29760) combined_score ≥ medium (see meta).
       Enrichment = Tier A∪B fraction among neighbors vs proteome background.
       Combined scores mix experiments, databases, coexpression, and <em>textmining</em> —
       not LLPS evidence. <a href="/network">Enrichment table</a></p>
    <p>STRING id: <span class="mono">${esc(nb.string_id || "—")}</span>
       · neighbors: <strong>${esc(nb.n_neighbors)}</strong>
       · Tier A∪B: ${esc(nb.n_neighbors_tier_ab)}
         (${nb.frac_tier_ab != null ? (100 * Number(nb.frac_tier_ab)).toFixed(1) : "—"}%)
       · fold: <strong>${fold}</strong>
       · hypergeom p: ${fmtP(nb.hypergeom_pvalue)}
       ${nb.is_seed ? "· <span class=\"muted\">Tier A seed</span>" : ""}</p>
    ${
      rows
        ? `<div class="table-wrap"><table>
            <thead><tr><th>Neighbor</th><th>Gene</th><th>Tier</th><th>PlantLLPS</th><th>STRING</th></tr></thead>
            <tbody>${rows}</tbody>
          </table></div>`
        : `<p class="muted">No neighbors above threshold.</p>`
    }
  </div>`;
}

app.get("/network", (req, res) => {
  const order = String(req.query.order || "pvalue") === "fold" ? "fold" : "pvalue";
  const minN = Math.max(1, Number(req.query.min_n) || 5);
  const page = Math.max(1, Number(req.query.page) || 1);
  const { rows, total, available } = listStringEnrichment(db, {
    seedsOnly: true,
    minNeighbors: minN,
    limit: PAGE_SIZE,
    offset: (page - 1) * PAGE_SIZE,
    order,
  });
  const sn = stringNetworkMeta(db);
  const pages = Math.max(1, Math.ceil(total / PAGE_SIZE));
  const tableRows = rows
    .map((r) => {
      const fold =
        r.enrichment_fold != null && Number.isFinite(Number(r.enrichment_fold))
          ? Number(r.enrichment_fold).toFixed(2)
          : "—";
      return `<tr>
        <td class="mono"><a href="/protein/${esc(r.protein_id)}">${esc(r.protein_id)}</a></td>
        <td>${esc(r.gene_name || "—")}</td>
        <td><span class="tier ${esc(r.tier || "")}">${esc(r.tier || "—")}</span></td>
        <td>${fmtScore(r.plantllps_score)}</td>
        <td>${esc(r.n_neighbors)}</td>
        <td>${esc(r.n_neighbors_tier_ab)}</td>
        <td>${fold}</td>
        <td>${fmtP(r.hypergeom_pvalue)}</td>
      </tr>`;
    })
    .join("");
  const body = `
    <div class="card">
      <h2>STRING neighborhood enrichment</h2>
      <p>1-hop neighbors of <strong>Tier A</strong> seeds in native <em>Vitis vinifera</em>
         STRING (NCBI ${esc(sn.taxon || "29760")}), combined_score ≥
         <strong>${esc(sn.min_score || "400")}</strong>. Mapped
         ${esc(sn.mapped || "—")} / 25 972 atlas IDs; deposited seed edges:
         ${esc(sn.seed_edges || "—")}.</p>
      <p class="muted">Fold = (neighbor Tier A∪B fraction) / (proteome A∪B rate ≈ 0.149).
         Hypergeometric p is descriptive only — STRING channels include textmining;
         this is <strong>not</strong> experimental condensate membership.
         Gold LLPS overlay (/gold) is separate; seeds here are score tiers.</p>
      <p>
        Sort:
        <a href="/network?order=pvalue&min_n=${minN}">p-value</a> ·
        <a href="/network?order=fold&min_n=${minN}">fold</a>
        · min neighbors:
        <a href="/network?order=${order}&min_n=3">3</a> /
        <a href="/network?order=${order}&min_n=5">5</a> /
        <a href="/network?order=${order}&min_n=10">10</a>
      </p>
      ${
        !available
          ? `<p class="muted">STRING tables not in this DB. Run
             <span class="mono">python3 scripts/enrich_string_neighborhood.py</span>.</p>`
          : total === 0
            ? `<p class="muted">No seed neighborhoods with ≥${esc(minN)} neighbors.</p>`
            : `<p class="muted">${esc(total)} seed rows (page ${esc(page)} / ${esc(pages)}).</p>
               <div class="table-wrap"><table>
                 <thead><tr>
                   <th>Protein</th><th>Gene</th><th>Tier</th><th>PlantLLPS</th>
                   <th>N</th><th>N A∪B</th><th>Fold</th><th>p</th>
                 </tr></thead>
                 <tbody>${tableRows}</tbody>
               </table></div>
               ${pagerHtml("/network", page, pages, { order, min_n: minN })}`
      }
    </div>`;
  res.type("html").send(layout("STRING network", body));
});

function pagerHtml(base, page, pages, query = {}) {
  if (pages <= 1) return "";
  const q = new URLSearchParams(
    Object.fromEntries(
      Object.entries(query).map(([k, v]) => [k, String(v)])
    )
  );
  const mk = (p) => {
    const qq = new URLSearchParams(q);
    qq.set("page", String(p));
    return `${base}?${qq.toString()}`;
  };
  const prev = page > 1 ? `<a class="btn secondary" href="${esc(mk(page - 1))}">← Prev</a>` : "";
  const next =
    page < pages ? `<a class="btn secondary" href="${esc(mk(page + 1))}">Next →</a>` : "";
  return `<div class="pager" style="margin-top:1rem">${prev} ${next}</div>`;
}

app.get("/gold", (req, res) => {
  const page = Math.max(1, Number(req.query.page || 1));
  const { rows, total } = listLlpsGold(db, {
    limit: PAGE_SIZE,
    offset: (page - 1) * PAGE_SIZE,
  });
  let curatedMeta = "";
  try {
    const m = getMeta(db);
    const g = m.llps_curated_gold ? JSON.parse(m.llps_curated_gold) : null;
    if (g) {
      curatedMeta = `Curated table rows: <strong>${esc(g.n_curated_rows)}</strong>
        · grape proteins marked: <strong>${esc(g.n_grape_proteins_marked)}</strong>`;
    }
  } catch {
    /* ignore */
  }
  const body = `
    <div class="card">
      <h2>Curated experimentally anchored LLPS</h2>
      <p class="muted">${curatedMeta || "Run scripts/enrich_curated_llps_gold.py"}</p>
      <p>Filter download: <a href="/download?evidence=llps">/download?evidence=llps</a>
         · source TSV <span class="mono">data/curated/llps_experimental_gold.tsv</span></p>
    </div>
    ${proteinTable(rows, "No curated LLPS joins in DB.")}
    ${pager("/gold", page, total)}`;
  res.type("html").send(layout("LLPS gold", body));
});

app.get("/map", (_req, res) => {
  const sidePath = join(__dirname, "data", "map_screen", "map_screen_spearman.json");
  let report = null;
  if (existsSync(sidePath)) {
    try {
      report = JSON.parse(readFileSync(sidePath, "utf8"));
    } catch {
      report = null;
    }
  }
  const ann = report?.metrics?.annotation;
  const go = report?.metrics?.go8017_only;
  const body = `
    <div class="card">
      <h2>MAP screen (reconstructed)</h2>
      <p class="muted">Companion Supplemental Data S1 (208 GO-annotated MAPs + original PSAP columns)
        was never deposited. This page summarizes the <strong>annotation reconstruction</strong>
        in <span class="mono">data/map_screen/</span> with native proteome PSAP / LLPhyScore.</p>
      ${
        ann
          ? `<p>Annotation set n=<strong>${esc(ann.n)}</strong>
               (GO:0008017 n=${esc(ann.n_go8017)}) ·
               ρ(companion vs native PSAP)≈<strong>${esc(
                 Number(ann.spearman_companion_vs_psap.rho).toFixed(3)
               )}</strong> ·
               ρ(PlantLLPS vs PSAP)≈<strong>${esc(
                 Number(ann.spearman_plantllps_vs_psap.rho).toFixed(3)
               )}</strong></p>
             <p class="muted">Companion-cited ρ=0.368 remains the manuscript transfer-gap number;
               recomputed ρ uses atlas native PSAP (see map_screen README).</p>
             <p>F6HET4 native PSAP rank ${esc(ann.F6HET4.psap_rank)} / ${esc(
                 ann.F6HET4.psap_n
               )}
               · LLPhyScore rank ${esc(ann.F6HET4.llphyscore_rank)} / ${esc(
                 ann.F6HET4.llphyscore_n
               )}
               · PlantLLPS ${fmtScore(ann.F6HET4.score_plantllps)}</p>
             <p class="muted">GO-only n=${esc(go?.n)} ·
               ρ(companion vs LLPhyScore)≈${esc(
                 Number(go?.spearman_companion_vs_llphyscore?.rho ?? NaN).toFixed(3)
               )}</p>`
          : `<p class="err">Missing data/map_screen/map_screen_spearman.json — run scripts/deposit_map_screen.py</p>`
      }
      <p><a href="/download?keyword=Microtubule">Download KW-Microtubule</a> ·
         <a href="/protein/F6HET4">VvMAP65-1 F6HET4</a></p>
    </div>`;
  res.type("html").send(layout("MAP screen", body));
});


app.get("/compare", (_req, res) => {
  // Prefer deposited gold side-by-side (includes Ath/human IDs not in grape proteome)
  const sidePath = join(__dirname, "data", "benchmarks", "comparator_side_by_side.json");
  let goldRows = [];
  if (existsSync(sidePath)) {
    try {
      goldRows = JSON.parse(readFileSync(sidePath, "utf8")).rows || [];
    } catch {
      goldRows = [];
    }
  }
  const goldCells = goldRows
    .map((r) => {
      const pid = esc(r.protein_id);
      const inDb = getProtein(db, r.protein_id);
      const link = inDb
        ? `<a href="/protein/${pid}">${pid}</a>`
        : `<span title="not in grape proteome DB">${pid}</span>`;
      return `<tr>
        <td class="mono">${link}</td>
        <td>${esc(r.label || "—")}</td>
        <td class="muted">${esc(r.note || "")}</td>
        <td>${fmtScore(r.score_plantllps)}</td>
        <td>${fmtScore(r.score_psap)}</td>
        <td>${fmtScore(r.score_llphyscore_percentile)}</td>
        <td>${fmtScore(r.plaac_nllr_approx)}</td>
        <td>${fmtScore(r.frac_disorder_metapredict)}</td>
        <td>${fmtScore(r.af_frac_very_low)}</td>
      </tr>`;
    })
    .join("");
  const meta = getMeta(db);
  let psapCov = "—";
  let llphyCov = "—";
  try {
    const psapMeta = meta.psap_proteome ? JSON.parse(meta.psap_proteome) : null;
    const llphyMeta = meta.llphyscore_proteome ? JSON.parse(meta.llphyscore_proteome) : null;
    const nProt = Number(meta.n_scored || 25972);
    if (psapMeta?.n != null) psapCov = `${psapMeta.n} / ${nProt}`;
    if (llphyMeta?.n != null) llphyCov = `${llphyMeta.n} / ${nProt}`;
  } catch {
    /* ignore */
  }
  // Proteome sample with native comparators (top PlantLLPS scores)
  const proteomeSample = db
    .prepare(
      `SELECT protein_id, gene_name, tier, score, score_psap, score_llphyscore_percentile,
              plaac_nllr_approx, frac_disorder_metapredict, af_frac_very_low FROM proteins
       ORDER BY score DESC LIMIT 40`
    )
    .all();
  const protCells = proteomeSample
    .map((r) => {
      const pid = esc(r.protein_id);
      return `<tr>
        <td class="mono"><a href="/protein/${pid}">${pid}</a></td>
        <td class="mono">${esc(r.gene_name || "—")}</td>
        <td><span class="tier ${esc(r.tier || "?")}">${esc(r.tier || "?")}</span></td>
        <td>${fmtScore(r.score)}</td>
        <td>${fmtScore(r.score_psap)}</td>
        <td>${fmtScore(r.score_llphyscore_percentile)}</td>
        <td>${fmtScore(r.plaac_nllr_approx)}</td>
        <td>${fmtScore(r.frac_disorder_metapredict)}</td>
        <td>${fmtScore(r.af_frac_very_low)}</td>
      </tr>`;
    })
    .join("");
  const body = `
    <div class="card">
      <h2>Side-by-side comparators</h2>
      <p class="muted">Gold panel: PlantLLPS vs native PSAP / LLPhyScore + PLAAC-approx.
      Proteome-wide native coverage — PSAP: <strong>${esc(psapCov)}</strong>;
      LLPhyScore: <strong>${esc(llphyCov)}</strong>.
      See <span class="mono">docs/EXTERNAL_TOOLS.md</span>.</p>
    </div>
    <h2 style="font-family:var(--font-display)">Gold standard (13+15)</h2>
    <div class="card" style="overflow-x:auto">
      <table>
        <thead><tr>
          <th>Protein</th><th>Label</th><th>Note</th><th>PlantLLPS</th>
          <th>PSAP</th><th>LLPhyScore %ile</th><th>PLAAC≈</th>
          <th>metapredict</th><th>AF pLDDT&lt;50</th>
        </tr></thead>
        <tbody>${goldCells || `<tr><td colspan="9" class="muted">Run scripts/enrich_comparator_scores.py</td></tr>`}</tbody>
      </table>
    </div>
    <h2 style="font-family:var(--font-display)">Proteome top-40 (PlantLLPS rank)</h2>
    <div class="card" style="overflow-x:auto">
      <table>
        <thead><tr>
          <th>Protein</th><th>Gene</th><th>Tier</th><th>PlantLLPS</th>
          <th>PSAP</th><th>LLPhyScore %ile</th><th>PLAAC≈</th>
          <th>metapredict</th><th>AF pLDDT&lt;50</th>
        </tr></thead>
        <tbody>${protCells}</tbody>
      </table>
    </div>`;
  res.type("html").send(layout("Compare", body));
});

const DOWNLOAD_TSV_COLS = [
  "protein_id",
  "gene_name",
  "tier",
  "score",
  "score_companion",
  "score_logistic",
  "score_psap",
  "score_llphyscore_percentile",
  "plaac_nllr_approx",
  "frac_disorder_metapredict",
  "af_mean_plddt",
  "af_frac_very_low",
  "domain_pfam_ids",
  "cross_crop_family",
  "keywords",
  "llps_gold",
];

function downloadFilters(query) {
  return {
    tier: String(query.tier || "").trim(),
    family: String(query.family || "").trim(),
    keyword: String(query.keyword || query.keywords || "").trim(),
    q: String(query.q || "").trim(),
    evidence: String(query.evidence || "").trim(),
    limit: Number(query.limit || 0) || 0,
  };
}

function rowsToTsv(rows) {
  const lines = [DOWNLOAD_TSV_COLS.join("\t")];
  for (const r of rows) {
    lines.push(
      DOWNLOAD_TSV_COLS.map((c) => {
        let v = r[c];
        if (c === "domain_pfam_ids" && Array.isArray(v)) v = v.join(";");
        if (c === "keywords" && Array.isArray(v)) v = v.join(";");
        if (v == null) return "";
        return String(v).replace(/[\t\n\r]/g, " ");
      }).join("\t")
    );
  }
  return lines.join("\n") + "\n";
}

function rowsToFasta(rows) {
  const chunks = [];
  let missing = 0;
  for (const r of rows) {
    const seq = seqById.get(r.protein_id);
    if (!seq) {
      missing += 1;
      continue;
    }
    const gene = r.gene_name || "";
    chunks.push(`>${r.protein_id}${gene ? ` ${gene}` : ""} tier=${r.tier} score=${r.score}`);
    for (let i = 0; i < seq.length; i += 60) chunks.push(seq.slice(i, i + 60));
  }
  return { text: chunks.join("\n") + (chunks.length ? "\n" : ""), missing };
}

app.get("/api/download", (req, res) => {
  const filters = downloadFilters(req.query);
  const format = String(req.query.format || "tsv").toLowerCase();
  const rows = filterProteins(db, filters);
  const stamp = [
    filters.tier && `tier${filters.tier}`,
    filters.family && `fam-${filters.family}`,
    filters.keyword && `kw`,
    filters.q && `q`,
  ]
    .filter(Boolean)
    .join("_") || "proteome";
  if (format === "fasta" || format === "fa") {
    const { text, missing } = rowsToFasta(rows);
    res.setHeader("Content-Type", "text/x-fasta; charset=utf-8");
    res.setHeader(
      "Content-Disposition",
      `attachment; filename="grape_condensatome_${stamp}.fasta"`
    );
    if (missing) res.setHeader("X-Missing-Sequences", String(missing));
    return res.send(text);
  }
  res.setHeader("Content-Type", "text/tab-separated-values; charset=utf-8");
  res.setHeader(
    "Content-Disposition",
    `attachment; filename="grape_condensatome_${stamp}.tsv"`
  );
  res.send(rowsToTsv(rows));
});

app.get("/download", (req, res) => {
  const filters = downloadFilters(req.query);
  const preview = filterProteins(db, { ...filters, limit: filters.limit || 20 });
  const count = filterProteins(db, filters).length;
  const qs = new URLSearchParams(
    Object.fromEntries(
      Object.entries(filters).filter(([k, v]) => v && k !== "limit")
    )
  ).toString();
  const body = `
    <div class="card">
      <h2>Bulk download</h2>
      <p class="muted">Filter the scored proteome by tier, cross-crop family, keyword, curated LLPS evidence, or free-text,
      then download TSV (score columns) or FASTA (sequences from UP000009183).
      API: <span class="mono">GET /api/download?tier=A&amp;evidence=llps&amp;format=tsv|fasta</span></p>
      <form method="get" action="/download" class="inline">
        <label>Tier
          <select name="tier">
            <option value="">Any</option>
            ${["A", "B", "C", "D"]
              .map(
                (t) =>
                  `<option value="${t}" ${filters.tier === t ? "selected" : ""}>${t}</option>`
              )
              .join("")}
          </select>
        </label>
        <label>Evidence
          <select name="evidence">
            <option value="">Any</option>
            <option value="llps" ${filters.evidence === "llps" ? "selected" : ""}>LLPS gold join</option>
          </select>
        </label>
        <label>Family <input name="family" value="${esc(filters.family)}" placeholder="MAP65"/></label>
        <label>Keyword <input name="keyword" value="${esc(filters.keyword)}" placeholder="Kinase"/></label>
        <label>Search <input name="q" value="${esc(filters.q)}" placeholder="VIT_ / Pfam / name"/></label>
        <button type="submit">Preview</button>
      </form>
      <p style="margin-top:0.75rem">${count} proteins match.
        <a href="/api/download?${qs}&format=tsv">Download TSV</a> ·
        <a href="/api/download?${qs}&format=fasta">Download FASTA</a>
      </p>
    </div>
    <h2 style="font-family:var(--font-display)">Preview (up to 20)</h2>
    ${proteinTable(preview)}`;
  res.type("html").send(layout("Download", body, filters.q || filters.family || filters.keyword));
});

app.get("/hooks", (req, res) => {
  const meta = getMeta(db);
  let hooks = { hooks: [] };
  try {
    hooks = meta.cross_crop_hooks ? JSON.parse(meta.cross_crop_hooks) : hooks;
  } catch {
    /* ignore */
  }
  let soyJoin = null;
  try {
    soyJoin = meta.second_crop_soybean ? JSON.parse(meta.second_crop_soybean) : null;
  } catch {
    soyJoin = null;
  }
  const page = Math.max(1, Number(req.query.page) || 1);
  const soyTier = String(req.query.soy_tier || "").trim().toUpperCase();
  const dbHooks = db
    .prepare(
      `SELECT protein_id, gene_name, tier, score, cross_crop_family, ortholog_ath, ortholog_ath_gene,
              ortholog_osa, ortholog_osa_gene, ortholog_gma, ortholog_gma_gene, ortholog_gma_score, ortholog_gma_tier
       FROM proteins WHERE cross_crop_family IS NOT NULL AND cross_crop_family != ''
       ORDER BY cross_crop_family, score DESC`
    )
    .all();
  const hookCards = (hooks.hooks || [])
    .map((h) => {
      const g = (h.grape || []).map((x) => `<a class="mono" href="/protein/${esc(x)}">${esc(x)}</a>`).join(", ") || "—";
      const a = (h.arabidopsis || []).map((x) => esc(x)).join(", ") || "—";
      const r = (h.rice || []).map((x) => esc(x)).join(", ") || "—";
      const soy = ((h.other_crops || {}).soybean || [])
        .map((x) => `<span class="mono">${esc(x)}</span>`)
        .join(", ") || "—";
      return `<div class="card" style="margin-bottom:0.75rem">
        <h3>${esc(h.family)}</h3>
        <p class="muted">${esc(h.note || "")}</p>
        <p>Grape: ${g}</p>
        <p>Arabidopsis: <span class="mono">${a}</span></p>
        <p>Rice: <span class="mono">${r}</span></p>
        <p>Soybean: ${soy}</p>
      </div>`;
    })
    .join("");

  let soyPanel = `<div class="card"><p class="muted">Soybean condensatome SQLite not loaded — run
    <span class="mono">scripts/score_second_crop_proteome.py</span>.</p></div>`;
  if (soyDb) {
    const clauses = [];
    const params = [];
    if (["A", "B", "C", "D"].includes(soyTier)) {
      clauses.push("tier = ?");
      params.push(soyTier);
    }
    const where = clauses.length ? `WHERE ${clauses.join(" AND ")}` : "";
    const total = soyDb.prepare(`SELECT COUNT(*) AS n FROM proteins ${where}`).get(...params).n;
    const offset = (page - 1) * PAGE_SIZE;
    const soyRows = soyDb
      .prepare(
        `SELECT protein_id, gene_name, protein_name, score, tier, percentile, grape_ortholog, grape_ortholog_gene
         FROM proteins ${where} ORDER BY score DESC LIMIT ? OFFSET ?`
      )
      .all(...params, PAGE_SIZE, offset);
    const tierCounts = {};
    for (const r of soyDb.prepare("SELECT tier, COUNT(*) AS c FROM proteins GROUP BY tier").all()) {
      tierCounts[r.tier] = r.c;
    }
    const soyTable = soyRows
      .map((r) => {
        const grape = r.grape_ortholog
          ? `<a class="mono" href="/protein/${esc(r.grape_ortholog)}">${esc(r.grape_ortholog)}</a>
             ${r.grape_ortholog_gene ? `(${esc(r.grape_ortholog_gene)})` : ""}`
          : "—";
        return `<tr>
          <td class="mono"><a href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(r.protein_id)}/entry" target="_blank" rel="noopener">${esc(r.protein_id)}</a></td>
          <td class="mono">${esc(r.gene_name || "—")}</td>
          <td><span class="tier ${esc(r.tier)}">${esc(r.tier)}</span></td>
          <td>${fmtScore(r.score)}</td>
          <td>${fmtPct(r.percentile)}</td>
          <td>${grape}</td>
          <td class="muted">${esc((r.protein_name || "").slice(0, 72) || "—")}</td>
        </tr>`;
      })
      .join("");
    const base = soyTier ? `/hooks?soy_tier=${encodeURIComponent(soyTier)}` : "/hooks";
    soyPanel = `
      <div class="card">
        <h2>Soybean condensatome <span class="mono">UP000008827</span></h2>
        <p class="muted">PlantLLPS v13 (composition-proxy disorder) scored on <em>Glycine max</em> reference
        proteome; tiers recalibrated on this proteome. Not a soybean-retrained model.
        OrthoDB join to grape: <strong>${esc(soyJoin?.n_matched ?? "—")}</strong> / ${esc(soyJoin?.n_grape ?? "—")}
        (${esc(soyJoin?.n_matched_with_crop_score ?? "—")} with crop scores).</p>
        <p>Tiers:
          ${["A", "B", "C", "D"]
            .map(
              (t) =>
                `<a href="/hooks?soy_tier=${t}">${t}=${esc(tierCounts[t] || 0)}</a>`
            )
            .join(" · ")}
          · <a href="/hooks">all</a>
          · n=${esc(total)}${soyTier ? ` (filter Tier ${esc(soyTier)})` : ""}
        </p>
      </div>
      <div class="card" style="overflow-x:auto">
        <table>
          <thead><tr>
            <th>Soybean ID</th><th>Gene</th><th>Tier</th><th>PlantLLPS</th><th>Percentile</th>
            <th>Grape ortholog</th><th>Protein name</th>
          </tr></thead>
          <tbody>${soyTable || `<tr><td colspan="7" class="muted">No rows</td></tr>`}</tbody>
        </table>
      </div>
      ${pager(base, page, total)}`;
  }

  const body = `
    <div class="card"><h2>Cross-crop condensatome hooks</h2>
    <p class="muted">Curated family anchors plus OrthoDB-joined Ath/rice/<strong>soybean</strong> columns on grape rows.
    Full second-crop atlas deposited under <span class="mono">data/second_crop/</span>.</p></div>
    ${hookCards || `<p class="muted">No hooks in meta — load data/cross_crop_hooks.json via enrich_orthologs.py</p>`}
    <h2 style="font-family:var(--font-display)">Grape proteome rows with family tag</h2>
    ${proteinTable(dbHooks.map((r) => getProtein(db, r.protein_id)).filter(Boolean))}
    <h2 style="font-family:var(--font-display);margin-top:1.5rem">Second crop · soybean</h2>
    ${soyPanel}`;
  res.type("html").send(layout("Cross-crop", body));
});

function scoreForm(resultHtml = "") {
  return `
    <div class="card">
      <h2>Score a sequence</h2>
      <p class="muted">Paste FASTA or raw amino-acid sequence. Scored with the bundled PlantLLPS v13 model
      (not written back to the proteome DB).</p>
      <form method="post" action="/score">
        <p><label>Optional ID<br/>
          <input type="text" name="protein_id" placeholder="QUERY1" style="min-width:12rem"/>
        </label></p>
        <p><label>Sequence<br/>
          <textarea name="sequence" placeholder=">QUERY1&#10;MKT..."></textarea>
        </label></p>
        <button type="submit">Score</button>
      </form>
    </div>
    ${resultHtml}`;
}

app.get("/score", (_req, res) => {
  res.type("html").send(layout("Score FASTA", scoreForm()));
});

app.post("/score", (req, res) => {
  const raw = String(req.body.sequence || "");
  const pid = String(req.body.protein_id || "QUERY").trim() || "QUERY";
  const lines = [];
  for (const line of raw.split(/\r?\n/)) {
    const t = line.trim();
    if (!t || t.startsWith(">")) continue;
    lines.push(t);
  }
  const seq = (lines.length ? lines.join("") : raw).replace(/[^A-Za-z]/g, "").toUpperCase();
  if (seq.length < 20) {
    return res
      .type("html")
      .send(
        layout(
          "Score FASTA",
          scoreForm(
            `<div class="card err"><p>Need a protein sequence of at least 20 amino acids.</p></div>`
          )
        )
      );
  }
  const scored = scoreSequence(seq, model);
  const liveExpl = explainFeatures(model, scored);
  const featHtml = Object.keys(scored)
    .filter(
      (k) =>
        !["score_plantllps", "score_companion", "score_logistic", "tier"].includes(k)
    )
    .sort()
    .map(
      (k) =>
        `<div><span class="k">${esc(k)}</span><span class="v">${esc(scored[k])}</span></div>`
    )
    .join("");
  const result = `
    <div class="card" style="margin-top:1rem">
      <h3>Result for <span class="mono">${esc(pid)}</span></h3>
      <p>Length ${seq.length} · Tier
         <span class="tier ${esc(scored.tier)}">${esc(scored.tier)}</span></p>
      <p>PlantLLPS: <strong>${fmtScore(scored.score_plantllps)}</strong>
         · companion ${fmtScore(scored.score_companion)}
         · logistic ${fmtScore(scored.score_logistic)}</p>
      <div class="features">${featHtml}</div>
    </div>
    ${explainabilityHtml(
      {
        grammar_axes: liveExpl.grammar.axes,
        top_features: liveExpl.top_features,
        method_summary:
          "Logistic: exact linear SHAP (φ=w·z). Companion terms shown beside bars.",
      },
      { length: seq.length, live: true }
    )}`;
  res.type("html").send(layout("Score FASTA", scoreForm(result)));
});

Promise.all([
  loadFastaIndex(FASTA_PATH).catch((err) => {
    console.error("FASTA index failed:", err);
  }),
  loadStringTopNeighbors().catch((err) => {
    console.error("STRING top neighbors failed:", err);
  }),
  loadStressSidecar().catch((err) => {
    console.error("Stress sidecar failed:", err);
  }),
  Promise.resolve().then(() => loadStressModulesUi()),
  loadExplainability().catch((err) => {
    console.error("Explainability sidecar failed:", err);
  }),
]).finally(() => {
  app.listen(PORT, HOST, () => {
    console.log(
      `GrapeCondensatome → http://${HOST === "0.0.0.0" ? "127.0.0.1" : HOST}:${PORT}/`
    );
    console.log(`Bound ${HOST}:${PORT}`);
    console.log(`SQLite: ${DB_PATH}`);
    console.log(`Model: ${MODEL_PATH}`);
    console.log(`FASTA sequences indexed: ${seqById.size}`);
    console.log(`Explainability attributions: ${explainById.size}`);
  });
});
