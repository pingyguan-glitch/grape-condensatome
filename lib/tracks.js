/**
 * Per-residue track sidecar (metapredict disorder + AlphaFold pLDDT) + InterPro names.
 */
import { DatabaseSync } from "node:sqlite";
import { createGunzip } from "node:zlib";
import { createReadStream, existsSync, readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = dirname(fileURLToPath(import.meta.url));
export const DEFAULT_TRACKS_DB = join(
  __dirname,
  "..",
  "data",
  "tracks",
  "grape_residue_tracks.sqlite"
);
export const DEFAULT_CANDIDATE_TRACKS_DB = join(
  __dirname,
  "..",
  "data",
  "tracks",
  "crop_candidates_residue_tracks.sqlite"
);
export const DEFAULT_INTERPRO_ENTRIES = join(
  __dirname,
  "..",
  "data",
  "tracks",
  "interpro_entries.json.gz"
);

/** @type {DatabaseSync|null} */
let tracksDb = null;
/** @type {DatabaseSync|null} */
let candidateTracksDb = null;
/** @type {Map<string, {accession:string,name:string,type:string,source_database?:string}>|null} */
let interproEntries = null;
/** @type {object|null} */
let tracksCoverage = null;
/** @type {object|null} */
let candidateTracksCoverage = null;

function openOneTracksDb(path, label) {
  if (!existsSync(path)) {
    console.warn(`${label} missing: ${path}`);
    return null;
  }
  const db = new DatabaseSync(path, { readOnly: true });
  const n = db.prepare("SELECT COUNT(*) AS n FROM residue_tracks").get()?.n ?? 0;
  console.log(`${label}: ${n} proteins from ${path}`);
  return db;
}

export function openTracksDb(
  path = DEFAULT_TRACKS_DB,
  candidatePath = DEFAULT_CANDIDATE_TRACKS_DB
) {
  if (!tracksDb) {
    tracksDb = openOneTracksDb(path, "Residue tracks (grape)");
    if (tracksDb) {
      try {
        const row = tracksDb.prepare("SELECT value FROM meta WHERE key = 'coverage'").get();
        tracksCoverage = row ? JSON.parse(row.value) : null;
      } catch {
        tracksCoverage = null;
      }
    }
  }
  if (!candidateTracksDb) {
    candidateTracksDb = openOneTracksDb(candidatePath, "Residue tracks (crop candidates)");
    if (candidateTracksDb) {
      try {
        const row = candidateTracksDb
          .prepare("SELECT value FROM meta WHERE key = 'coverage'")
          .get();
        candidateTracksCoverage = row ? JSON.parse(row.value) : null;
      } catch {
        candidateTracksCoverage = null;
      }
    }
  }
  return tracksDb;
}

export function getTracksCoverage() {
  return {
    grape: tracksCoverage,
    crop_candidates: candidateTracksCoverage,
  };
}

export async function loadInterproEntries(path = DEFAULT_INTERPRO_ENTRIES) {
  if (interproEntries) return interproEntries;
  interproEntries = new Map();
  if (!existsSync(path)) {
    console.warn(`InterPro entries cache missing: ${path}`);
    return interproEntries;
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
  for (const [acc, meta] of Object.entries(obj.entries || {})) {
    interproEntries.set(acc, meta);
  }
  console.log(`InterPro entries: ${interproEntries.size} from ${path}`);
  return interproEntries;
}

export function getInterproEntry(acc) {
  return interproEntries?.get(acc) || null;
}

function decodeUnit8(blob) {
  if (!blob) return null;
  const buf = Buffer.isBuffer(blob) ? blob : Buffer.from(blob);
  const out = new Array(buf.length);
  for (let i = 0; i < buf.length; i++) out[i] = buf[i] / 255;
  return out;
}

function decodePlddt(blob) {
  if (!blob) return null;
  const buf = Buffer.isBuffer(blob) ? blob : Buffer.from(blob);
  const out = new Array(buf.length);
  for (let i = 0; i < buf.length; i++) out[i] = buf[i];
  return out;
}

/** Downsample for SVG (keep endpoints). */
export function downsample(values, maxPoints = 800) {
  if (!values || values.length <= maxPoints) return values ? [...values] : [];
  const out = [];
  const n = values.length;
  for (let i = 0; i < maxPoints; i++) {
    const idx = Math.round((i * (n - 1)) / (maxPoints - 1));
    out.push(values[idx]);
  }
  return out;
}

function rowToTracks(row) {
  if (!row) return null;
  return {
    protein_id: row.protein_id,
    length: row.length,
    disorder: decodeUnit8(row.disorder_u8),
    plddt: decodePlddt(row.plddt_u8),
    disorder_source: row.disorder_source,
    plddt_source: row.plddt_source,
    disorder_mean: row.disorder_mean,
    plddt_mean: row.plddt_mean,
  };
}

const TRACK_SQL = `SELECT protein_id, length, disorder_u8, plddt_u8,
              disorder_source, plddt_source, disorder_mean, plddt_mean
       FROM residue_tracks WHERE protein_id = ? COLLATE NOCASE`;

export function getResidueTracks(proteinId) {
  const id = String(proteinId).trim();
  if (tracksDb) {
    const row = tracksDb.prepare(TRACK_SQL).get(id);
    if (row) return rowToTracks(row);
  }
  if (candidateTracksDb) {
    const row = candidateTracksDb.prepare(TRACK_SQL).get(id);
    if (row) return rowToTracks(row);
  }
  return null;
}

export function annotateInterproIds(ids) {
  return (ids || []).map((acc) => {
    const meta = getInterproEntry(acc);
    return {
      accession: acc,
      name: meta?.name || "",
      type: meta?.type || "",
      source_database: meta?.source_database || "interpro",
    };
  });
}
