import { DatabaseSync } from "node:sqlite";
import { createGunzip } from "node:zlib";
import { createReadStream, existsSync, readFileSync } from "node:fs";
import { createInterface } from "node:readline";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = dirname(fileURLToPath(import.meta.url));
export const DEFAULT_DB = join(__dirname, "..", "data", "grape_condensatome.sqlite");
export const DEFAULT_STRING_TOP = join(
  __dirname,
  "..",
  "data",
  "string",
  "top_neighbors.tsv.gz"
);
export const DEFAULT_STRING_STATS = join(
  __dirname,
  "..",
  "data",
  "string",
  "neighborhood_stats.tsv.gz"
);
export const DEFAULT_STRESS_SIDECAR = join(
  __dirname,
  "..",
  "data",
  "expression_atlas",
  "stress_by_protein.json.gz"
);
export const DEFAULT_STRESS_MODULES_UI = join(
  __dirname,
  "..",
  "data",
  "stress_modules",
  "ui_tables.json"
);
export const DEFAULT_FAMILY_ENRICHMENT_UI = join(
  __dirname,
  "..",
  "data",
  "family_enrichment",
  "ui_tables.json"
);
export const DEFAULT_COMPARATOR_DISAGREEMENT_UI = join(
  __dirname,
  "..",
  "data",
  "comparator_disagreement",
  "ui_tables.json"
);
export const DEFAULT_STRUCTURE_SCORE_UI = join(
  __dirname,
  "..",
  "data",
  "structure_score",
  "ui_tables.json"
);

/** @type {Map<string, Array<{protein_id:string,combined_score:number,tier?:string,score?:number,gene_name?:string}>>|null} */
let stringTopCache = null;
/** @type {Map<string, object>|null} */
let stringStatsCache = null;
/** @type {Map<string, object>|null} */
let stressByProtein = null;
/** @type {object|null} */
let stressModulesUi = null;
/** @type {object|null} */
let familyEnrichmentUi = null;
/** @type {object|null} */
let comparatorDisagreementUi = null;
/** @type {object|null} */
let structureScoreUi = null;

export async function loadStressSidecar(path = DEFAULT_STRESS_SIDECAR) {
  if (!existsSync(path)) {
    stressByProtein = new Map();
    console.warn(`Stress sidecar missing: ${path}`);
    return stressByProtein;
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
  stressByProtein = new Map(Object.entries(obj));
  console.log(`Stress sidecar: ${stressByProtein.size} proteins from ${path}`);
  return stressByProtein;
}

export function loadStressModulesUi(path = DEFAULT_STRESS_MODULES_UI) {
  if (!existsSync(path)) {
    stressModulesUi = null;
    console.warn(`Stress modules UI tables missing: ${path}`);
    return null;
  }
  try {
    stressModulesUi = JSON.parse(readFileSync(path, "utf8"));
    const nMod = stressModulesUi?.modules?.length ?? 0;
    const nEnr = stressModulesUi?.enrichment?.length ?? 0;
    console.log(`Stress modules UI: ${nEnr} enrichment rows, ${nMod} modules from ${path}`);
    return stressModulesUi;
  } catch (err) {
    stressModulesUi = null;
    console.warn(`Failed to load stress modules UI: ${err}`);
    return null;
  }
}

export function getStressModulesUi() {
  return stressModulesUi;
}

export function stressModulesMeta(db) {
  const meta = getMeta(db);
  let summary = null;
  if (meta.stress_condensatome_modules) {
    try {
      summary = JSON.parse(meta.stress_condensatome_modules);
    } catch {
      summary = null;
    }
  }
  return summary;
}

function loadJsonUi(path, label) {
  if (!existsSync(path)) {
    console.warn(`${label} UI tables missing: ${path}`);
    return null;
  }
  try {
    const obj = JSON.parse(readFileSync(path, "utf8"));
    console.log(`${label} UI loaded from ${path}`);
    return obj;
  } catch (err) {
    console.warn(`Failed to load ${label} UI: ${err}`);
    return null;
  }
}

export function loadAtlasAnalysesUi() {
  familyEnrichmentUi = loadJsonUi(DEFAULT_FAMILY_ENRICHMENT_UI, "Family enrichment");
  comparatorDisagreementUi = loadJsonUi(
    DEFAULT_COMPARATOR_DISAGREEMENT_UI,
    "Comparator disagreement"
  );
  structureScoreUi = loadJsonUi(DEFAULT_STRUCTURE_SCORE_UI, "Structure–score");
  return {
    family: familyEnrichmentUi,
    disagreement: comparatorDisagreementUi,
    structure: structureScoreUi,
  };
}

export function getFamilyEnrichmentUi() {
  return familyEnrichmentUi;
}

export function getComparatorDisagreementUi() {
  return comparatorDisagreementUi;
}

export function getStructureScoreUi() {
  return structureScoreUi;
}

function parseMetaJson(db, key) {
  const meta = getMeta(db);
  if (!meta[key]) return null;
  try {
    return JSON.parse(meta[key]);
  } catch {
    return null;
  }
}

export function familyEnrichmentMeta(db) {
  return parseMetaJson(db, "family_enrichment");
}

export function comparatorDisagreementMeta(db) {
  return parseMetaJson(db, "comparator_disagreement");
}

export function structureScoreMeta(db) {
  return parseMetaJson(db, "structure_score_joint");
}

export function openDb(path = DEFAULT_DB) {
  const db = new DatabaseSync(path, { readOnly: true });
  return db;
}

function normalizeSegments(segs) {
  if (!Array.isArray(segs)) return [];
  return segs.map((s) => {
    const start = s.start ?? s.s;
    const end = s.end ?? s.e;
    const length =
      s.length ?? (start != null && end != null ? end - start + 1 : undefined);
    return {
      ...s,
      start,
      end,
      length,
      source: s.source ?? s.src ?? "",
      label: s.label ?? s.lab,
    };
  });
}

function parseAnnotationBlob(raw) {
  if (!raw) return null;
  try {
    const obj = JSON.parse(raw);
    if (obj && typeof obj === "object") return obj;
  } catch {
    /* plain text legacy */
  }
  return null;
}

function rowToProtein(row) {
  if (!row) return null;
  let features = {};
  try {
    features = JSON.parse(row.features_json || "{}");
  } catch {
    features = {};
  }
  const annObj = parseAnnotationBlob(row.annotations);
  const gene_name = row.gene_name || annObj?.gene_name || "";
  const protein_name = row.protein_name || annObj?.protein_name || "";
  const go_terms = row.go_terms
    ? String(row.go_terms).split(";").filter(Boolean)
    : annObj?.go_terms || [];
  const keywords = row.keywords
    ? String(row.keywords).split(";").filter(Boolean)
    : annObj?.keywords || [];
  let expr_tpm = null;
  if (row.expr_tpm_json) {
    try {
      expr_tpm = JSON.parse(row.expr_tpm_json);
    } catch {
      expr_tpm = null;
    }
  }
  let expr_stress = null;
  if (row.expr_stress_json) {
    try {
      expr_stress = JSON.parse(row.expr_stress_json);
    } catch {
      expr_stress = null;
    }
  }
  if (!expr_stress && stressByProtein?.has(row.protein_id)) {
    expr_stress = stressByProtein.get(row.protein_id);
  }
  let af_disorder_segments = [];
  if (row.af_disorder_segments_json) {
    try {
      af_disorder_segments = normalizeSegments(
        JSON.parse(row.af_disorder_segments_json)
      );
    } catch {
      af_disorder_segments = [];
    }
  }
  let domains = [];
  if (row.domains_json) {
    try {
      domains = normalizeSegments(JSON.parse(row.domains_json) || []);
    } catch {
      domains = [];
    }
  }
  let idr_segments = [];
  if (row.idr_segments_json) {
    try {
      idr_segments = normalizeSegments(JSON.parse(row.idr_segments_json) || []);
    } catch {
      idr_segments = [];
    }
  }
  if (!af_disorder_segments.length && idr_segments.length) {
    af_disorder_segments = idr_segments
      .filter((s) => String(s.source || "").includes("af_plddt"))
      .map((s) => ({
        start: s.start,
        end: s.end,
        length: s.length,
        source: s.source,
      }));
  }
  const domain_pfam_ids = row.domain_pfam_ids
    ? String(row.domain_pfam_ids).split(";").filter(Boolean)
    : [];
  const domain_interpro_ids = row.domain_interpro_ids
    ? String(row.domain_interpro_ids).split(";").filter(Boolean)
    : [];
  let ortholog_obj = null;
  if (row.ortholog_json) {
    try {
      ortholog_obj = JSON.parse(row.ortholog_json);
    } catch {
      ortholog_obj = null;
    }
  }
  const metaFrac =
    row.frac_disorder_metapredict != null
      ? row.frac_disorder_metapredict
      : features.frac_disorder_metapredict;
  let llps_gold_entries = [];
  if (row.llps_gold_json) {
    try {
      llps_gold_entries = JSON.parse(row.llps_gold_json) || [];
    } catch {
      llps_gold_entries = [];
    }
  }
  let ext_llps = null;
  if (row.ext_llps_json) {
    try {
      ext_llps = JSON.parse(row.ext_llps_json);
    } catch {
      ext_llps = null;
    }
  }
  const ext_llps_dbs = row.ext_llps_dbs
    ? String(row.ext_llps_dbs)
        .split(";")
        .map((s) => s.trim())
        .filter(Boolean)
    : [];
  return {
    ...row,
    features,
    percentile: features.percentile,
    map65_tail_like: features.map65_tail_like ?? 0,
    frac_disorder_proxy: features.frac_disorder_proxy,
    frac_disorder_metapredict: metaFrac,
    gene_name,
    protein_name,
    go_terms,
    keywords,
    annotation_obj: annObj,
    expr_tpm,
    expr_stress,
    af_disorder_segments,
    domains,
    idr_segments,
    domain_pfam_ids,
    domain_interpro_ids,
    ortholog_obj,
    llps_gold: Number(row.llps_gold || 0),
    llps_gold_entries,
    ext_llps_dbs,
    ext_llps_unique: Number(row.ext_llps_unique ?? 1),
    ext_llps,
    plaac_nllr_approx:
      row.plaac_nllr_approx != null ? row.plaac_nllr_approx : features.plaac_nllr_approx,
  };
}

/** Filter proteome rows for bulk download (tier / family / keyword / q / evidence). */
export function filterProteins(
  db,
  { tier = "", family = "", keyword = "", q = "", evidence = "", limit = 0 } = {}
) {
  const clauses = [];
  const params = [];
  const t = String(tier || "").trim().toUpperCase();
  if (t && ["A", "B", "C", "D"].includes(t)) {
    clauses.push("tier = ?");
    params.push(t);
  }
  const fam = String(family || "").trim();
  if (fam) {
    clauses.push("IFNULL(cross_crop_family, '') LIKE ? COLLATE NOCASE");
    params.push(`%${fam}%`);
  }
  const kw = String(keyword || "").trim();
  if (kw) {
    clauses.push("IFNULL(keywords, '') LIKE ? COLLATE NOCASE");
    params.push(`%${kw}%`);
  }
  const ev = String(evidence || "").trim().toLowerCase();
  if (ev === "llps" || ev === "gold" || ev === "1") {
    clauses.push("IFNULL(llps_gold, 0) > 0");
  } else if (
    ev === "unique" ||
    ev === "unique_to_us" ||
    ev === "grape_unique" ||
    ev === "plantllps_unique"
  ) {
    clauses.push("IFNULL(ext_llps_unique, 1) = 1");
  } else if (ev === "shared" || ev === "catalogued" || ev === "external") {
    clauses.push("IFNULL(ext_llps_unique, 1) = 0");
  } else if (["phasepro", "llpsdb", "drllps", "cdcode"].includes(ev)) {
    clauses.push("IFNULL(ext_llps_dbs, '') LIKE ?");
    params.push(`%${ev}%`);
  }
  const query = String(q || "").trim();
  if (query) {
    clauses.push(
      `(protein_id LIKE ? COLLATE NOCASE OR IFNULL(gene_name, '') LIKE ? COLLATE NOCASE
        OR IFNULL(protein_name, '') LIKE ? COLLATE NOCASE
        OR IFNULL(keywords, '') LIKE ? COLLATE NOCASE
        OR IFNULL(cross_crop_family, '') LIKE ? COLLATE NOCASE
        OR IFNULL(domain_pfam_ids, '') LIKE ? COLLATE NOCASE
        OR IFNULL(llps_gold_json, '') LIKE ? COLLATE NOCASE
        OR IFNULL(ext_llps_dbs, '') LIKE ? COLLATE NOCASE)`
    );
    const like = `%${query}%`;
    params.push(like, like, like, like, like, like, like, like);
  }
  const where = clauses.length ? `WHERE ${clauses.join(" AND ")}` : "";
  const lim = Number(limit) > 0 ? `LIMIT ${Number(limit)}` : "";
  const sql = `SELECT * FROM proteins ${where} ORDER BY score DESC ${lim}`;
  return db.prepare(sql).all(...params).map(rowToProtein);
}

export function llpsCrosswalkMeta(db) {
  return parseMetaJson(db, "llps_crosswalk");
}

export function listUniqueCandidates(
  db,
  { tiers = "A,B,C", limit = 50, offset = 0 } = {}
) {
  const keep = String(tiers || "A,B,C")
    .split(",")
    .map((t) => t.trim().toUpperCase())
    .filter((t) => ["A", "B", "C", "D"].includes(t));
  if (!keep.length) return { rows: [], total: 0 };
  const placeholders = keep.map(() => "?").join(",");
  const total = db
    .prepare(
      `SELECT COUNT(*) AS n FROM proteins
       WHERE IFNULL(ext_llps_unique, 1) = 1 AND tier IN (${placeholders})`
    )
    .get(...keep).n;
  const rows = db
    .prepare(
      `SELECT * FROM proteins
       WHERE IFNULL(ext_llps_unique, 1) = 1 AND tier IN (${placeholders})
       ORDER BY score DESC LIMIT ? OFFSET ?`
    )
    .all(...keep, limit, offset)
    .map(rowToProtein);
  return { rows, total };
}

export function getMeta(db) {
  const rows = db.prepare("SELECT key, value FROM meta").all();
  return Object.fromEntries(rows.map((r) => [r.key, r.value]));
}

export function summarize(db) {
  const n = db.prepare("SELECT COUNT(*) AS n FROM proteins").get().n;
  const tiers = {};
  for (const r of db.prepare("SELECT tier, COUNT(*) AS c FROM proteins GROUP BY tier").all()) {
    tiers[r.tier] = r.c;
  }
  const top = db
    .prepare("SELECT protein_id, score, tier FROM proteins ORDER BY score DESC LIMIT 5")
    .all();
  return { n_proteins: n, tiers, top };
}

export function getProtein(db, proteinId) {
  const row = db
    .prepare("SELECT * FROM proteins WHERE protein_id = ? COLLATE NOCASE")
    .get(String(proteinId).trim());
  return rowToProtein(row);
}

export function searchProteins(db, query, { limit = 50, offset = 0 } = {}) {
  const q = String(query || "").trim();
  if (!q) return { rows: [], total: 0 };
  const exact = db
    .prepare("SELECT * FROM proteins WHERE protein_id = ? COLLATE NOCASE")
    .get(q);
  if (exact) return { rows: [rowToProtein(exact)], total: 1 };
  const like = `%${q}%`;
  const total = db
    .prepare(
      `SELECT COUNT(*) AS n FROM proteins
       WHERE protein_id LIKE ? COLLATE NOCASE
          OR IFNULL(annotations, '') LIKE ? COLLATE NOCASE
          OR IFNULL(gene_name, '') LIKE ? COLLATE NOCASE
          OR IFNULL(protein_name, '') LIKE ? COLLATE NOCASE
          OR IFNULL(go_terms, '') LIKE ? COLLATE NOCASE
          OR IFNULL(keywords, '') LIKE ? COLLATE NOCASE`
    )
    .get(like, like, like, like, like, like).n;
  const rows = db
    .prepare(
      `SELECT * FROM proteins
       WHERE protein_id LIKE ? COLLATE NOCASE
          OR IFNULL(annotations, '') LIKE ? COLLATE NOCASE
          OR IFNULL(gene_name, '') LIKE ? COLLATE NOCASE
          OR IFNULL(protein_name, '') LIKE ? COLLATE NOCASE
          OR IFNULL(go_terms, '') LIKE ? COLLATE NOCASE
          OR IFNULL(keywords, '') LIKE ? COLLATE NOCASE
       ORDER BY score DESC LIMIT ? OFFSET ?`
    )
    .all(like, like, like, like, like, like, limit, offset)
    .map(rowToProtein);
  return { rows, total };
}

export function listByTier(db, tier, { limit = 50, offset = 0 } = {}) {
  const t = String(tier || "").trim().toUpperCase();
  if (!["A", "B", "C", "D"].includes(t)) return { rows: [], total: 0 };
  const total = db.prepare("SELECT COUNT(*) AS n FROM proteins WHERE tier = ?").get(t).n;
  const rows = db
    .prepare(
      `SELECT * FROM proteins WHERE tier = ?
       ORDER BY score DESC LIMIT ? OFFSET ?`
    )
    .all(t, limit, offset)
    .map(rowToProtein);
  return { rows, total };
}

export function topProteins(db, limit = 20) {
  return db
    .prepare("SELECT * FROM proteins ORDER BY score DESC LIMIT ?")
    .all(limit)
    .map(rowToProtein);
}

function tableExists(db, name) {
  const row = db
    .prepare("SELECT 1 AS ok FROM sqlite_master WHERE type='table' AND name = ?")
    .get(name);
  return Boolean(row);
}

async function readTsvGz(path, onRow) {
  if (!existsSync(path)) {
    console.warn(`STRING deposit not found: ${path}`);
    return 0;
  }
  const stream = path.endsWith(".gz")
    ? createReadStream(path).pipe(createGunzip())
    : createReadStream(path);
  const rl = createInterface({ input: stream, crlfDelay: Infinity });
  let header = null;
  let n = 0;
  for await (const line of rl) {
    if (!line || line.startsWith("#")) continue;
    const parts = line.split("\t");
    if (!header) {
      header = parts;
      continue;
    }
    onRow(parts, header);
    n += 1;
  }
  return n;
}

/** Load full-proteome top-N + enrichment stats sidecars (keeps SQLite under 100 MB). */
export async function loadStringTopNeighbors(
  topPath = DEFAULT_STRING_TOP,
  statsPath = DEFAULT_STRING_STATS
) {
  stringTopCache = new Map();
  stringStatsCache = new Map();
  const nTop = await readTsvGz(topPath, (parts) => {
    const pid = parts[0];
    const neighbor = {
      protein_id: parts[1],
      combined_score: Number(parts[2]),
      tier: parts[4] || "",
      score: parts[5] != null && parts[5] !== "" ? Number(parts[5]) : null,
      gene_name: parts[6] || "",
    };
    if (!stringTopCache.has(pid)) stringTopCache.set(pid, []);
    stringTopCache.get(pid).push(neighbor);
  });
  const nStats = await readTsvGz(statsPath, (parts) => {
    const pid = parts[0];
    stringStatsCache.set(pid, {
      protein_id: pid,
      string_id: parts[1] || null,
      n_neighbors: Number(parts[2]),
      n_neighbors_tier_a: Number(parts[3]),
      n_neighbors_tier_ab: Number(parts[4]),
      frac_tier_ab: parts[5] !== "" ? Number(parts[5]) : null,
      enrichment_fold: parts[6] !== "" ? Number(parts[6]) : null,
      hypergeom_pvalue: parts[7] !== "" ? Number(parts[7]) : null,
      mean_neighbor_score: parts[8] !== "" ? Number(parts[8]) : null,
      mean_neighbor_combined: parts[9] !== "" ? Number(parts[9]) : null,
      is_seed: Number(parts[10]) || 0,
      tier: parts[11] || "",
      plantllps_score: parts[12] !== "" ? Number(parts[12]) : null,
    });
  });
  console.log(
    `STRING deposits: top-neighbors proteins=${stringTopCache.size} (rows≈${nTop}); stats=${stringStatsCache.size} (rows≈${nStats})`
  );
  return stringTopCache;
}

function topNeighborsFor(proteinId) {
  if (!stringTopCache) return [];
  const key = String(proteinId).trim();
  return stringTopCache.get(key) || stringTopCache.get(key.toUpperCase()) || [];
}

function statsFor(proteinId) {
  if (!stringStatsCache) return null;
  const key = String(proteinId).trim();
  return stringStatsCache.get(key) || stringStatsCache.get(key.toUpperCase()) || null;
}

/** STRING 1-hop neighborhood summary for one protein (SQLite seed row or TSV stats). */
export function getStringNeighborhood(db, proteinId) {
  const pid = String(proteinId).trim();
  let row = null;
  if (tableExists(db, "string_neighborhood")) {
    row = db
      .prepare("SELECT * FROM string_neighborhood WHERE protein_id = ? COLLATE NOCASE")
      .get(pid);
  }
  if (row) {
    let top = [];
    try {
      top = JSON.parse(row.top_neighbors_json || "[]") || [];
    } catch {
      top = [];
    }
    if (!top.length) top = topNeighborsFor(row.protein_id);
    return { ...row, top_neighbors: top };
  }
  const st = statsFor(pid);
  if (!st) return null;
  return { ...st, top_neighbors: topNeighborsFor(pid) };
}

/** Enrichment table for /network (default: Tier A seeds with enough neighbors). */
export function listStringEnrichment(
  db,
  { seedsOnly = true, minNeighbors = 5, limit = 50, offset = 0, order = "pvalue" } = {}
) {
  // Prefer in-memory full stats when loaded; fall back to SQLite seed table.
  if (stringStatsCache && stringStatsCache.size) {
    let rows = [...stringStatsCache.values()].filter((r) => r.n_neighbors >= minNeighbors);
    if (seedsOnly) rows = rows.filter((r) => r.is_seed);
    rows.sort((a, b) => {
      if (order === "fold") {
        return (b.enrichment_fold ?? -1) - (a.enrichment_fold ?? -1);
      }
      return (a.hypergeom_pvalue ?? 1) - (b.hypergeom_pvalue ?? 1);
    });
    const total = rows.length;
    const page = rows.slice(offset, offset + limit).map((r) => {
      const p = getProtein(db, r.protein_id);
      return {
        ...r,
        gene_name: p?.gene_name || "",
        tier: r.tier || p?.tier || "",
        plantllps_score: r.plantllps_score ?? p?.score,
      };
    });
    return { rows: page, total, available: true };
  }
  if (!tableExists(db, "string_neighborhood")) {
    return { rows: [], total: 0, available: false };
  }
  const clauses = ["n_neighbors >= ?"];
  const params = [minNeighbors];
  if (seedsOnly) clauses.push("is_seed = 1");
  const where = `WHERE ${clauses.join(" AND ")}`;
  const total = db.prepare(`SELECT COUNT(*) AS n FROM string_neighborhood ${where}`).get(...params).n;
  const orderSql =
    order === "fold"
      ? "enrichment_fold DESC, hypergeom_pvalue ASC"
      : "hypergeom_pvalue ASC, enrichment_fold DESC";
  const rows = db
    .prepare(
      `SELECT sn.*, p.gene_name, p.tier, p.score AS plantllps_score
       FROM string_neighborhood sn
       JOIN proteins p ON p.protein_id = sn.protein_id
       ${where}
       ORDER BY ${orderSql}
       LIMIT ? OFFSET ?`
    )
    .all(...params, limit, offset);
  return { rows, total, available: true };
}

export function stringNetworkMeta(db) {
  const meta = getMeta(db);
  let summary = null;
  if (meta.string_enrichment_summary) {
    try {
      summary = JSON.parse(meta.string_enrichment_summary);
    } catch {
      summary = null;
    }
  }
  return {
    version: meta.string_version || null,
    taxon: meta.string_taxon || null,
    min_score: meta.string_min_score || null,
    mapped: meta.string_mapped || null,
    seed_set: meta.string_seed_set || null,
    seed_edges: meta.string_seed_edges || null,
    neighborhood_n: meta.string_neighborhood_n || null,
    summary,
  };
}

export function listLlpsGold(db, { limit = 100, offset = 0 } = {}) {
  const total = db
    .prepare("SELECT COUNT(*) AS n FROM proteins WHERE IFNULL(llps_gold, 0) > 0")
    .get().n;
  const rows = db
    .prepare(
      `SELECT * FROM proteins WHERE IFNULL(llps_gold, 0) > 0
       ORDER BY score DESC LIMIT ? OFFSET ?`
    )
    .all(limit, offset)
    .map(rowToProtein);
  return { rows, total };
}
