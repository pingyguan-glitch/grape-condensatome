/** Empirical percentiles vs a deposited background score distribution. */

/**
 * Load sorted PlantLLPS scores from the grape proteome SQLite.
 * @param {import("node:sqlite").DatabaseSync} db
 * @returns {number[]}
 */
export function loadSortedScores(db) {
  const rows = db.prepare("SELECT score FROM proteins WHERE score IS NOT NULL").all();
  return rows
    .map((r) => Number(r.score))
    .filter((x) => Number.isFinite(x))
    .sort((a, b) => a - b);
}

/**
 * Percentile of `score` in a sorted ascending background array.
 * Definition: 100 * (# background scores ≤ query) / N.
 * Returns null when background is empty.
 */
export function empiricalPercentile(sortedAsc, score) {
  const n = sortedAsc?.length || 0;
  if (!n || !Number.isFinite(Number(score))) return null;
  const x = Number(score);
  let lo = 0;
  let hi = n;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (sortedAsc[mid] <= x) lo = mid + 1;
    else hi = mid;
  }
  return (100.0 * lo) / n;
}
