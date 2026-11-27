/** PlantLLPS SHAP-equivalent feature / grammar attributions (JS port). */

import { FEATURE_NAMES, featuresToVector } from "./features.js";
import { companionScore, logisticScore } from "./model.js";

export const FEATURE_GROUPS = {
  disorder: ["frac_disorder_proxy"],
  low_complexity: ["min_entropy_12", "low_complexity_score"],
  charge_blockiness: [
    "ncpr",
    "kappa_charge_approx",
    "blockiness_ncpr_sd11",
    "charge_frac",
    "pi_approx",
  ],
  aromatic: ["frac_aromatic"],
  ser_arg: ["frac_ser", "frac_arg", "frac_ser_arg"],
  sp_tp: ["sp_tp_density"],
  map65_cterm: ["c_term_ser_frac_100", "c_term_ncpr_100", "map65_tail_like"],
  qn_polar: ["qn_frac"],
  length: ["length"],
};

export const GRAMMAR_AXES = {
  ser_rich: ["ser_arg", "map65_cterm"],
  charge_block: ["charge_blockiness"],
  idr_disorder: ["disorder", "low_complexity"],
};

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

export function logisticShap(model, feats) {
  const names = model.feature_names || FEATURE_NAMES;
  const vec = featuresToVector(feats, names);
  if (!model?.trained || !model.weights?.length) {
    return {
      base_logit: 0,
      logit: 0,
      probability: companionScore(feats),
      feature_shap_logit: Object.fromEntries(names.map((n) => [n, 0])),
      feature_shap_prob_local: Object.fromEntries(names.map((n) => [n, 0])),
    };
  }
  const scaled = vec.map((v, i) => (v - model.means[i]) / (model.stds[i] || 1));
  const phi = {};
  let logit = model.bias || 0;
  for (let i = 0; i < names.length; i++) {
    const v = (model.weights[i] || 0) * scaled[i];
    phi[names[i]] = v;
    logit += v;
  }
  const p = sigmoid(logit);
  const scale = p * (1 - p);
  const phiProb = Object.fromEntries(
    Object.entries(phi).map(([k, v]) => [k, v * scale])
  );
  return {
    base_logit: model.bias || 0,
    logit,
    probability: p,
    feature_shap_logit: phi,
    feature_shap_prob_local: phiProb,
  };
}

export function companionDecomposition(feats) {
  const d = clip01(Number(feats.frac_disorder_proxy || 0));
  const lc = clip01(Number(feats.low_complexity_score || 0));
  const blk = clip01(Number(feats.blockiness_ncpr_sd11 || 0) / 0.2);
  const map65 = Number(feats.map65_tail_like || 0);
  const cSer = Number(feats.c_term_ser_frac_100 || 0);
  const terms = {
    frac_disorder_proxy: COMPOSITE_WEIGHTS.frac_disorder_proxy * d,
    low_complexity_score: COMPOSITE_WEIGHTS.low_complexity_score * lc,
    blockiness_ncpr_sd11: COMPOSITE_WEIGHTS.blockiness_ncpr_sd11 * blk,
    map65_tail_bonus: 0.08 * map65,
    c_term_ser_bonus: 0.05 * clip01(cSer / 0.15),
  };
  return { score: clip01(Object.values(terms).reduce((a, b) => a + b, 0)), terms };
}

function groupSum(phi, members) {
  return members.reduce((a, n) => a + Number(phi[n] || 0), 0);
}

export function grammarAttributions(featurePhi, companionTerms = {}) {
  const groups = Object.fromEntries(
    Object.entries(FEATURE_GROUPS).map(([g, members]) => [
      g,
      groupSum(featurePhi, members),
    ])
  );
  const companionGroup = Object.fromEntries(
    Object.keys(FEATURE_GROUPS).map((g) => [g, 0])
  );
  companionGroup.disorder = Number(companionTerms.frac_disorder_proxy || 0);
  companionGroup.low_complexity = Number(companionTerms.low_complexity_score || 0);
  companionGroup.charge_blockiness = Number(companionTerms.blockiness_ncpr_sd11 || 0);
  companionGroup.map65_cterm =
    Number(companionTerms.map65_tail_bonus || 0) +
    Number(companionTerms.c_term_ser_bonus || 0);

  const axes = {};
  for (const [axis, groupNames] of Object.entries(GRAMMAR_AXES)) {
    const shap = groupNames.reduce((a, g) => a + groups[g], 0);
    const comp = groupNames.reduce((a, g) => a + companionGroup[g], 0);
    axes[axis] = {
      logistic_shap_prob: shap,
      companion: comp,
      hybrid: shap,
      member_groups: groupNames,
    };
  }
  return { groups, companionGroup, axes };
}

/** Feature + grammar attributions for a scored feature dict + model. */
export function explainFeatures(model, feats) {
  const shap = logisticShap(model, feats);
  const comp = companionDecomposition(feats);
  const grammar = grammarAttributions(shap.feature_shap_prob_local, comp.terms);
  const ranked = (model.feature_names || FEATURE_NAMES)
    .map((name) => ({
      feature: name,
      value: Number(feats[name] || 0),
      shap_logit: shap.feature_shap_logit[name] || 0,
      shap_prob_local: shap.feature_shap_prob_local[name] || 0,
    }))
    .sort((a, b) => Math.abs(b.shap_prob_local) - Math.abs(a.shap_prob_local));
  return {
    score_logistic: logisticScore(model, feats),
    score_companion: companionScore(feats),
    logistic_shap: shap,
    companion: comp,
    grammar,
    top_features: ranked.slice(0, 12),
  };
}
