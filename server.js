#!/usr/bin/env node
/**
 * GrapeCondensatome — Express web UI over bundled SQLite + PlantLLPS JS scorer.
 * Listens on 0.0.0.0:$PORT (default 8765).
 */

import express from "express";
import { existsSync, readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import {
  openDb,
  getMeta,
  summarize,
  getProtein,
  searchProteins,
  listByTier,
  topProteins,
  DEFAULT_DB,
} from "./lib/db.js";
import { loadModel, scoreSequence } from "./lib/model.js";

const __dirname = dirname(fileURLToPath(import.meta.url));
const PORT = Number(process.env.PORT || 8765);
const HOST = process.env.HOST || "0.0.0.0";
const PAGE_SIZE = 40;
const DB_PATH = process.env.GRAPE_DB || DEFAULT_DB;
const MODEL_PATH =
  process.env.GRAPE_MODEL || join(__dirname, "data", "model_formal.json");

if (!existsSync(DB_PATH)) {
  console.error(`SQLite not found: ${DB_PATH}`);
  process.exit(1);
}

const db = openDb(DB_PATH);
const model = loadModel(MODEL_PATH);
const app = express();
app.use(express.urlencoded({ extended: true, limit: "2mb" }));
app.use(express.static(join(__dirname, "public")));

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
      <a href="/hooks">Cross-crop</a>
      <a href="/score">Score FASTA</a>
    </nav>
  </header>
  <div class="banner">
    PlantLLPS hybrid model (v13). Disorder: composition proxy by default; metapredict + AlphaFold pLDDT
    joined where computed. Orthologs (Ath/rice via OrthoDB) and side-by-side PSAP/LLPhyScore/PLAAC-approx
    on gold (+ proteome PLAAC). Treat ranks as hypotheses for experiment.
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
      return `<tr>
        <td class="mono"><a href="/protein/${pid}">${pid}</a></td>
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
        <th>PLAAC≈</th><th>Ath ortholog</th><th>Mean TPM</th><th>Protein name</th>
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

app.get("/api/health", (_req, res) => {
  res.json({
    ok: true,
    version: "1.0.0",
    db: DB_PATH,
    n: summarize(db).n_proteins,
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
  const uniprot = `https://www.uniprot.org/uniprotkb/${encodeURIComponent(p.protein_id)}/entry`;
  const go = (p.go_terms || []).slice(0, 24).map((g) => esc(g)).join(", ") || "—";
  const kw = (p.keywords || []).slice(0, 24).map((g) => esc(g)).join(", ") || "—";
  const body = `
    <div class="detail-hero">
      <div class="card">
        <h2 class="mono">${esc(p.protein_id)}</h2>
        <p><strong>${esc(p.gene_name || "—")}</strong>
           · <span class="tier ${esc(p.tier)}">${esc(p.tier)}</span>
           · length ${esc(p.length)}
           · <a href="${uniprot}" rel="noopener" target="_blank">UniProt</a></p>
        <p>${esc(p.protein_name || "No UniProt protein name in DB.")}</p>
        <p><span class="muted">GO:</span> ${go}</p>
        <p><span class="muted">Keywords:</span> ${kw}</p>
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
           <div class="features">${Object.entries(p.expr_tpm)
             .map(
               ([k, v]) =>
                 `<div><span class="k">${esc(k)}</span><span class="v">${esc(v)}</span></div>`
             )
             .join("")}</div></div>`
        : ""
    }
    <div class="card" style="margin-top:1rem">
      <h3>Features</h3>
      <div class="features">${featHtml}</div>
    </div>`;
  res.type("html").send(layout(p.protein_id, body));
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
  // Proteome sample with PLAAC (top scores) for side-by-side vs PlantLLPS
  const proteomeSample = db
    .prepare(
      `SELECT protein_id, gene_name, tier, score, plaac_nllr_approx, frac_disorder_metapredict,
              af_frac_very_low FROM proteins
       WHERE plaac_nllr_approx IS NOT NULL
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
        <td class="muted">—</td>
        <td class="muted">—</td>
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
      Proteome-wide PSAP/LLPhyScore are <strong>not</strong> shipped (external envs + CPU);
      PLAAC-approx and disorder tracks cover the atlas. See
      <span class="mono">docs/EXTERNAL_TOOLS.md</span>.</p>
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
    <h2 style="font-family:var(--font-display)">Proteome top-40 (PLAAC≈ / disorder)</h2>
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

app.get("/hooks", (_req, res) => {
  const meta = getMeta(db);
  let hooks = { hooks: [] };
  try {
    hooks = meta.cross_crop_hooks ? JSON.parse(meta.cross_crop_hooks) : hooks;
  } catch {
    /* ignore */
  }
  const dbHooks = db
    .prepare(
      `SELECT protein_id, gene_name, tier, score, cross_crop_family, ortholog_ath, ortholog_ath_gene,
              ortholog_osa, ortholog_osa_gene
       FROM proteins WHERE cross_crop_family IS NOT NULL AND cross_crop_family != ''
       ORDER BY cross_crop_family, score DESC`
    )
    .all();
  const hookCards = (hooks.hooks || [])
    .map((h) => {
      const g = (h.grape || []).map((x) => `<a class="mono" href="/protein/${esc(x)}">${esc(x)}</a>`).join(", ") || "—";
      const a = (h.arabidopsis || []).map((x) => esc(x)).join(", ") || "—";
      const r = (h.rice || []).map((x) => esc(x)).join(", ") || "—";
      return `<div class="card" style="margin-bottom:0.75rem">
        <h3>${esc(h.family)}</h3>
        <p class="muted">${esc(h.note || "")}</p>
        <p>Grape: ${g}</p>
        <p>Arabidopsis: <span class="mono">${a}</span></p>
        <p>Rice: <span class="mono">${r}</span></p>
      </div>`;
    })
    .join("");
  const body = `
    <div class="card"><h2>Cross-crop condensatome hooks</h2>
    <p class="muted">Curated family anchors plus OrthoDB-joined Ath/rice columns on proteome rows.</p></div>
    ${hookCards || `<p class="muted">No hooks in meta — load data/cross_crop_hooks.json via enrich_orthologs.py</p>`}
    <h2 style="font-family:var(--font-display)">Proteome rows with family tag</h2>
    ${proteinTable(dbHooks.map((r) => getProtein(db, r.protein_id)).filter(Boolean))}`;
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
    </div>`;
  res.type("html").send(layout("Score FASTA", scoreForm(result)));
});

app.listen(PORT, HOST, () => {
  console.log(`GrapeCondensatome → http://${HOST === "0.0.0.0" ? "127.0.0.1" : HOST}:${PORT}/`);
  console.log(`Bound ${HOST}:${PORT}`);
  console.log(`SQLite: ${DB_PATH}`);
  console.log(`Model: ${MODEL_PATH}`);
});
