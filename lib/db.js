import { DatabaseSync } from "node:sqlite";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = dirname(fileURLToPath(import.meta.url));
export const DEFAULT_DB = join(__dirname, "..", "data", "grape_condensatome.sqlite");

export function openDb(path = DEFAULT_DB) {
  const db = new DatabaseSync(path, { readOnly: true });
  return db;
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
  let af_disorder_segments = [];
  if (row.af_disorder_segments_json) {
    try {
      af_disorder_segments = JSON.parse(row.af_disorder_segments_json);
    } catch {
      af_disorder_segments = [];
    }
  }
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
    af_disorder_segments,
    ortholog_obj,
    plaac_nllr_approx:
      row.plaac_nllr_approx != null ? row.plaac_nllr_approx : features.plaac_nllr_approx,
  };
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
