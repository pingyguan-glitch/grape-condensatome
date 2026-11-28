/** FASTA parse + browser/server scoring caps (shared with /score + /api/score). */

export const SCORE_LIMITS = Object.freeze({
  max_sequences: 25,
  min_length: 20,
  max_length: 8000,
  max_total_residues: 50_000,
  max_fasta_chars: 1_500_000,
  /** Soft client expectation; server work is synchronous per request. */
  suggested_timeout_ms: 120_000,
});

/**
 * Parse FASTA text or a raw amino-acid blob into {id, sequence} records.
 * Header IDs take the first whitespace-delimited token after `>`;
 * UniProt-style `sp|ACC|NAME` / `tr|ACC|NAME` → ACC.
 */
export function parseFasta(text, { defaultId = "QUERY" } = {}) {
  const raw = String(text || "");
  if (!raw.trim()) return [];

  const records = [];
  let curId = null;
  let chunks = [];
  let sawHeader = false;

  const flush = () => {
    if (curId == null) return;
    const sequence = chunks
      .join("")
      .replace(/[^A-Za-z]/g, "")
      .toUpperCase();
    records.push({ id: curId, sequence });
    curId = null;
    chunks = [];
  };

  for (const line of raw.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      flush();
      sawHeader = true;
      const rest = line.slice(1).trim();
      curId = normalizeFastaId(rest) || `${defaultId}${records.length + 1}`;
      chunks = [];
    } else {
      const t = line.trim();
      if (!t) continue;
      if (curId == null) {
        // Raw sequence without header — accumulate as one record.
        curId = defaultId;
      }
      chunks.push(t);
    }
  }
  flush();

  // Entire paste was one raw sequence (no `>` lines).
  if (!sawHeader && records.length === 0) {
    const sequence = raw.replace(/[^A-Za-z]/g, "").toUpperCase();
    if (sequence) records.push({ id: defaultId, sequence });
  }

  return records;
}

function normalizeFastaId(headerRest) {
  if (!headerRest) return "";
  const token = headerRest.split(/\s+/)[0] || "";
  const m = token.match(/^(?:sp|tr)\|([^|]+)\|/i);
  if (m) return m[1];
  return token.replace(/[^\w.:-]/g, "_").slice(0, 80) || "";
}

/**
 * Validate parsed records against SCORE_LIMITS.
 * @returns {{ ok: true, records } | { ok: false, error: string, code: string }}
 */
export function validateScoreBatch(records, { fastaChars = 0 } = {}) {
  const L = SCORE_LIMITS;
  if (fastaChars > L.max_fasta_chars) {
    return {
      ok: false,
      code: "fasta_too_large",
      error: `FASTA paste exceeds ${L.max_fasta_chars.toLocaleString()} characters (got ${fastaChars.toLocaleString()}).`,
    };
  }
  if (!records.length) {
    return {
      ok: false,
      code: "empty",
      error: "No sequences found. Paste FASTA (with `>` headers) or a raw amino-acid sequence.",
    };
  }
  if (records.length > L.max_sequences) {
    return {
      ok: false,
      code: "too_many_sequences",
      error: `At most ${L.max_sequences} sequences per request (got ${records.length}). Split the batch or use the CLI \`plantllps score\`.`,
    };
  }

  let total = 0;
  for (const r of records) {
    const n = r.sequence.length;
    if (n < L.min_length) {
      return {
        ok: false,
        code: "too_short",
        error: `Sequence ${r.id} is too short (${n} aa; minimum ${L.min_length}).`,
      };
    }
    if (n > L.max_length) {
      return {
        ok: false,
        code: "too_long",
        error: `Sequence ${r.id} exceeds ${L.max_length} aa (got ${n}). Use the CLI for long ORFs.`,
      };
    }
    total += n;
  }
  if (total > L.max_total_residues) {
    return {
      ok: false,
      code: "too_many_residues",
      error: `Total residues ${total.toLocaleString()} exceed cap ${L.max_total_residues.toLocaleString()}.`,
    };
  }
  return { ok: true, records };
}
