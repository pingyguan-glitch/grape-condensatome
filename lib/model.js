/** PlantLLPS scoring from bundled model_formal.json (JS port). */

import { readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { FEATURE_NAMES, extractFeatures, featuresToVector } from "./features.js";

const __dirname = dirname(fileURLToPath(import.meta.url));
const DEFAULT_MODEL = join(__dirname, "..", "data", "model_formal.json");

const COMPOSITE_WEIGHTS = {
  frac_disorder_proxy: 0.55,
  low_complexity_score: 0.25,
  blockiness_ncpr_sd11: 0.2,
};

function clip01(x) {
  return Math.max(0, Math.min(1, x));
}

function sigmoid(x) {
  if (x >= 0) {
    const z = Math.exp(-x);
    return 1 / (1 + z);
  }
  const z = Math.exp(x);
  return z / (1 + z);
}

export function loadModel(path = DEFAULT_MODEL) {
  return JSON.parse(readFileSync(path, "utf8"));
}

export function companionScore(feats) {
  const d = clip01(Number(feats.frac_disorder_proxy || 0));
  const lc = clip01(Number(feats.low_complexity_score || 0));
  const blk = clip01(Number(feats.blockiness_ncpr_sd11 || 0) / 0.2);
  let bonus = 0.08 * Number(feats.map65_tail_like || 0);
  bonus += 0.05 * clip01(Number(feats.c_term_ser_frac_100 || 0) / 0.15);
  return clip01(
    COMPOSITE_WEIGHTS.frac_disorder_proxy * d +
      COMPOSITE_WEIGHTS.low_complexity_score * lc +
      COMPOSITE_WEIGHTS.blockiness_ncpr_sd11 * blk +
      bonus
  );
}

export function logisticScore(model, feats) {
  if (!model?.trained || !model.weights?.length) return companionScore(feats);
  const names = model.feature_names || FEATURE_NAMES;
  const vec = featuresToVector(feats, names);
  const scaled = vec.map((v, i) => (v - model.means[i]) / (model.stds[i] || 1));
  let z = model.bias || 0;
  for (let i = 0; i < scaled.length; i++) z += (model.weights[i] || 0) * scaled[i];
  return sigmoid(z);
}

export function assignTier(score, feats, cutoffs = {}) {
  const map65 = Number(feats.map65_tail_like || 0);
  const disorder = Number(feats.frac_disorder_proxy || 0);
  const cSer = Number(feats.c_term_ser_frac_100 || 0);
  if (map65 >= 0.5 && disorder < 0.7) return "B";
  if (score >= 0.45 && cSer >= 0.1 && disorder < 0.55 && map65 > 0) return "B";
  const aCut = Number(cutoffs.score_A ?? 0.65);
  const cCut = Number(cutoffs.score_C ?? 0.4);
  const dA = Number(cutoffs.disorder_A ?? 0.45);
  if (score >= aCut && disorder >= dA) return "A";
  if (score >= cCut) return "C";
  return "D";
}

export function scoreSequence(sequence, model) {
  const feats = extractFeatures(sequence);
  const comp = companionScore(feats);
  const logi = logisticScore(model, feats);
  const mode = model.mode || "hybrid";
  let plant;
  if (mode === "companion_composite") plant = comp;
  else if (mode === "logistic_stub") plant = logi;
  else plant = model.trained ? 0.5 * comp + 0.5 * logi : comp;
  const tier = assignTier(plant, feats, model.tier_cutoffs || {});
  return {
    ...feats,
    score_companion: comp,
    score_logistic: logi,
    score_plantllps: plant,
    tier,
  };
}
