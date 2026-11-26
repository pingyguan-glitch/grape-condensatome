/** PlantLLPS sequence features (JS port of plantllps/features.py). */

const AA = "ACDEFGHIKLMNPQRSTVWY";
const HYDROPHOBIC = new Set("AILMFVW");
const CHARGED = new Set("DEKR");
const POSITIVE = new Set("KR");
const NEGATIVE = new Set("DE");
const AROMATIC = new Set("FYW");
const DISORDER_PROMOTING = new Set("DESTPQKRNGHA");

export const FEATURE_NAMES = [
  "length",
  "frac_disorder_proxy",
  "min_entropy_12",
  "low_complexity_score",
  "ncpr",
  "kappa_charge_approx",
  "blockiness_ncpr_sd11",
  "frac_aromatic",
  "frac_ser",
  "frac_arg",
  "frac_ser_arg",
  "sp_tp_density",
  "qn_frac",
  "charge_frac",
  "pi_approx",
  "c_term_ser_frac_100",
  "c_term_ncpr_100",
  "map65_tail_like",
];

function cleanSeq(seq) {
  const s = String(seq || "")
    .toUpperCase()
    .replace(/[^A-Z]/g, "");
  return [...s].map((ch) => (AA.includes(ch) ? ch : "X")).join("");
}

function shannonEntropy(window) {
  if (!window) return 0;
  const n = window.length;
  const counts = {};
  for (const ch of window) counts[ch] = (counts[ch] || 0) + 1;
  let h = 0;
  for (const c of Object.values(counts)) {
    const p = c / n;
    h -= p * Math.log2(p);
  }
  return h;
}

function disorderProxyFraction(seq, win = 21) {
  if (!seq) return 0;
  const half = Math.floor(win / 2);
  const n = seq.length;
  let flags = 0;
  for (let i = 0; i < n; i++) {
    const a = Math.max(0, i - half);
    const b = Math.min(n, i + half + 1);
    const w = seq.slice(a, b);
    let promoting = 0;
    let buried = 0;
    for (const r of w) {
      if (DISORDER_PROMOTING.has(r)) promoting++;
      if (HYDROPHOBIC.has(r)) buried++;
    }
    const len = w.length;
    if (promoting / len >= 0.55 && buried / len <= 0.35) flags++;
  }
  return flags / n;
}

function ncpr(seq) {
  if (!seq) return 0;
  let pos = 0;
  let neg = 0;
  for (const r of seq) {
    if (POSITIVE.has(r)) pos++;
    else if (NEGATIVE.has(r)) neg++;
  }
  return (pos - neg) / seq.length;
}

function blockinessNcprSd(seq, win = 11) {
  if (seq.length < win) return 0;
  const vals = [];
  for (let i = 0; i <= seq.length - win; i++) {
    vals.push(ncpr(seq.slice(i, i + win)));
  }
  const mean = vals.reduce((a, b) => a + b, 0) / vals.length;
  const variance = vals.reduce((a, b) => a + (b - mean) ** 2, 0) / vals.length;
  return Math.sqrt(variance);
}

function kappaApprox(seq) {
  if (seq.length < 30) return 0;
  const obs = blockinessNcprSd(seq, 5);
  const charged = [...seq].filter((r) => CHARGED.has(r));
  const others = [...seq].filter((r) => !CHARGED.has(r));
  if (charged.length < 4) return 0;
  const plus = charged.filter((r) => POSITIVE.has(r));
  const minus = charged.filter((r) => NEGATIVE.has(r));
  const mixedCharged = [];
  while (plus.length || minus.length) {
    if (plus.length) mixedCharged.push(plus.pop());
    if (minus.length) mixedCharged.push(minus.pop());
  }
  const step = Math.max(1, Math.floor(seq.length / Math.max(1, mixedCharged.length)));
  const arr2 = Array(seq.length).fill("G");
  let oI = 0;
  let cI = 0;
  for (let i = 0; i < seq.length; i++) {
    if (i % step === 0 && cI < mixedCharged.length) {
      arr2[i] = mixedCharged[cI++];
    } else if (oI < others.length) {
      arr2[i] = others[oI++];
    } else if (cI < mixedCharged.length) {
      arr2[i] = mixedCharged[cI++];
    }
  }
  const mixedSd = blockinessNcprSd(arr2.join(""), 5);
  const seg = [...seq]
    .sort((a, b) => {
      const rank = (r) => (POSITIVE.has(r) ? 0 : NEGATIVE.has(r) ? 1 : 2);
      return rank(a) - rank(b);
    })
    .join("");
  const segSd = blockinessNcprSd(seg, 5);
  const denom = segSd - mixedSd;
  if (denom <= 1e-9) return 0;
  return Math.max(0, Math.min(1, (obs - mixedSd) / denom));
}

function spTpDensity(seq) {
  if (seq.length < 2) return 0;
  let count = 0;
  for (let i = 0; i < seq.length - 1; i++) {
    if ((seq[i] === "S" || seq[i] === "T") && seq[i + 1] === "P") count++;
  }
  return count / (seq.length - 1);
}

export function extractFeatures(sequence) {
  const seq = cleanSeq(sequence);
  const n = seq.length;
  if (n === 0) {
    return Object.fromEntries(FEATURE_NAMES.map((k) => [k, 0]));
  }
  const counts = {};
  for (const ch of seq) counts[ch] = (counts[ch] || 0) + 1;
  const frac = (aa) => (counts[aa] || 0) / n;

  const win = 12;
  let minEnt;
  if (n >= win) {
    let m = Infinity;
    for (let i = 0; i <= n - win; i++) {
      m = Math.min(m, shannonEntropy(seq.slice(i, i + win)));
    }
    minEnt = m;
  } else {
    minEnt = shannonEntropy(seq);
  }
  const maxEnt = Math.log2(Math.min(win, 20));
  const lcScore = 1.0 - (maxEnt > 0 ? minEnt / maxEnt : 0);

  const fracSer = frac("S");
  const fracArg = frac("R");
  const fracAro = frac("F") + frac("Y") + frac("W");
  const cterm = n >= 50 ? seq.slice(-100) : seq;
  const cSer = [...cterm].filter((r) => r === "S").length / cterm.length;
  const cNcpr = ncpr(cterm);
  const cAro = [...cterm].filter((r) => AROMATIC.has(r)).length / cterm.length;

  let map65Like = 0;
  if (cSer >= 0.08 && cNcpr > 0 && cNcpr < 0.15 && cAro < 0.05 && n >= 200) {
    map65Like = 1;
  } else if (cSer >= 0.1 && cNcpr > 0 && n >= 150) {
    map65Like = 0.5;
  }

  return {
    length: n,
    frac_disorder_proxy: disorderProxyFraction(seq),
    min_entropy_12: minEnt,
    low_complexity_score: lcScore,
    ncpr: ncpr(seq),
    kappa_charge_approx: kappaApprox(seq),
    blockiness_ncpr_sd11: blockinessNcprSd(seq, 11),
    frac_aromatic: fracAro,
    frac_ser: fracSer,
    frac_arg: fracArg,
    frac_ser_arg: fracSer + fracArg,
    sp_tp_density: spTpDensity(seq),
    qn_frac: frac("Q") + frac("N"),
    charge_frac: frac("D") + frac("E") + frac("K") + frac("R"),
    pi_approx: 7.0 + 12.0 * ncpr(seq),
    c_term_ser_frac_100: cSer,
    c_term_ncpr_100: cNcpr,
    map65_tail_like: map65Like,
  };
}

export function featuresToVector(feats, names = FEATURE_NAMES) {
  return names.map((n) => Number(feats[n] || 0));
}
