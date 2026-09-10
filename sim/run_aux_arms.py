"""
E-A auxiliary arms, at (c = 0.5, mid red noise) ONLY -- the spec's scope bound D16.

ARM 1 -- REVOCATION vs EQUIVALENT STATIC COVERAGE.
  The semantics leg of conjunct C-II.  Three conditions, matched seeds:
    static_full    c = 0.5 for the whole run (reference, higher coverage);
    revocation     c = 0.5, and at cycle R a centrality-weighted 40% of consenters
                   revoke.  Removal is RETROACTIVE: the instrument re-derives its
                   whole indicator history from the currently-consented set, so the
                   already-published pre-revocation series is rewritten;
    static_equiv   the SAME agents non-consenting from cycle 0.  Identical final
                   consent set, identical achieved coverage, no history rewrite --
                   so the only difference is the timing of the removal.
  Because the revoking set is drawn at consent assignment and applied later, the run
  carries both the as-published stream and its post-revocation re-derivation, and
  static_equiv is exactly matched rather than approximately matched.

  Two contrasts are reported:
   (i)  DETECTION TIMING -- the observer scans the as-published stream up to the
        revocation, then the re-derived stream afterwards (it cannot un-see an alarm
        it already raised, and it cannot keep a record it must erase);
   (ii) AGGREGATE INVALIDATION -- how far the published pre-revocation indicator
        series is from its own post-revocation re-derivation.  A static coverage
        reduction cannot have this property at all: there is nothing to rewrite.

  Retroactive propagation is a DESIGN CHOICE stricter than GDPR Art. 7 requires, and
  is distinct from graph unlearning (Chen et al. 2022 repairs a trained model's
  parameters; revocation holes the observed structure a third-party indicator is
  computed from, and there is no model to retrain).

  A null result on (i) costs the revocation-SEMANTICS sub-claim only (spec section 6,
  kill criterion 1); C-II's ownership/capability leg does not depend on this arm.

ARM 2 -- THE ENSIGN FIX.
  Inverse-propensity correction of the observation process with the propensity
  estimated FROM OBSERVED CONSENT STATUS ALONE: logistic regression of consent on the
  agent's degree in the OBSERVED subgraph, the only covariate a consent-bounded
  observer has.  That estimator is frozen here, so the arm is falsifiable.

  Where the propensity is load-bearing: the observer sees an observed program's
  CONSENTED members but not its non-consenting ones, so it does not know the team size.
  It must infer the team size from the estimated propensity (m_hat = k / p_hat), undo
  the partial-observation shrinkage of the recorded latency span at m_hat, and
  inverse-weight the program by the probability that a team of m_hat would have been
  observed at all.

  Reported: WHERE IT FAILS.  Under revocation the consent record is destroyed
  retroactively, so revoked agents enter the fit labelled as never-consenting; p_hat is
  then the POST-revocation consent rate applied to PRE-revocation records generated at a
  higher rate, m_hat is overestimated, and the correction over-corrects.  The estimator
  has no way to know it is wrong.

Run:  python run_aux_arms.py [--workers N] [--reps 60]
      python run_aux_arms.py --rule detection_rule_recal.json --out aux_arms_recal.json
Writes: aggregates/aux_arms.json (or --out)

FOLLOW-UP CAMPAIGN 2026-08-05: --rule re-scores the identical runs (seeds depend on the
replicate only) under the recalibrated threshold; --out keeps campaign 1's aggregate
intact so both operating points stay verifiable.
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from model import Config, run_config
from scoring import score_run, bootstrap_mean
from indicators import (auc, logistic_fit, logistic_predict, baseline_z,
                        alarm_from_z, lagging_alarm)
from runner import pmap, seed_for, write_agg, load_rule

C_AUX = 0.5
COLOUR_AUX = "red"
SEV_AUX = "mid"
CONDITIONS = ("static_full", "revocation", "static_equiv")

_RULE = None


def _rule():
    global _RULE
    if _RULE is None:
        _RULE = load_rule()
    return _RULE


def _cfg(condition: str, scen: str, seed: int) -> Config:
    base = dict(consent_c=C_AUX, consent_arm="mnar", eps_colour=COLOUR_AUX,
                eps_severity=SEV_AUX, scenario=scen, seed=seed)
    if condition == "revocation":
        return Config(revocation=True, revocation_retroactive=True, **base)
    if condition == "static_equiv":
        # same agents removed, but from cycle 0: identical final consent set,
        # identical achieved coverage, nothing to re-derive
        return Config(revocation=True, revocation_retroactive=True,
                      revocation_cycle=0, **base)
    return Config(revocation=False, **base)


def _score_two_stream(s, rule, n_cycles: int, t_rev: int) -> dict:
    """Sequential observer across a retroactive rewrite: it scans the AS-PUBLISHED
    stream while that is the record, and the RE-DERIVED stream once the revocation has
    holed the record.  An alarm already raised stands; a record already erased is
    gone."""
    L = rule["leading"]
    start = int(L["start"])
    args = (int(L["window"]), tuple(L["classes"]), start,
            float(L["min_frac"]), int(L["max_gap"]))
    za = baseline_z(s, *args)
    sb = dict(s)
    for cls in L["classes"]:
        sb[cls] = s[cls + "_rev"]
    zb = baseline_z(sb, *args)
    pre = alarm_from_z(za["z_var"], za["z_ac"], float(L["k"]), int(L["sustain"]),
                       L["combine"], start)
    t_pre = pre["alarm_cycle"]
    if t_pre is not None and t_pre < t_rev:
        t_lead = t_pre
    else:
        post = alarm_from_z(zb["z_var"], zb["z_ac"], float(L["k"]),
                            int(L["sustain"]), L["combine"], max(start, t_rev))
        t_lead = post["alarm_cycle"]
    lag = lagging_alarm(s["output_close"], rule["lagging"])
    t_lag = lag["alarm_cycle"]
    eff_lead = n_cycles if t_lead is None else t_lead
    eff_lag = n_cycles if t_lag is None else t_lag
    # aggregate invalidation over the published-then-rewritten stretch
    if L["combine"] == "var_only":
        ia = za["z_var"]
        ib = zb["z_var"]
    elif L["combine"] == "mean":
        ia = 0.5 * (za["z_var"] + za["z_ac"])
        ib = 0.5 * (zb["z_var"] + zb["z_ac"])
    else:
        ia = np.minimum(za["z_var"], za["z_ac"])
        ib = np.minimum(zb["z_var"], zb["z_ac"])
    seg = slice(start, max(start + 1, t_rev))
    da, db = ia[seg], ib[seg]
    m = np.isfinite(da) & np.isfinite(db)
    inval = float(np.mean(np.abs(da[m] - db[m]))) if m.any() else float("nan")
    k = float(L["k"])
    flip = (float(np.mean((da[m] >= k) != (db[m] >= k))) if m.any()
            else float("nan"))
    return {
        "t_lead": t_lead, "t_lag": t_lag,
        "dt_censored": eff_lag - eff_lead,
        "dt_paired": (t_lag - t_lead) if (t_lead is not None
                                         and t_lag is not None) else None,
        "lead_max_score": max(pre["max_score"], float(np.nanmax(ib[start:]))
                              if np.isfinite(ib[start:]).any() else -np.inf),
        "lead_fired": t_lead is not None,
        "lag_fired": t_lag is not None,
        "aggregate_invalidation_mean_abs_z": inval,
        "aggregate_alarm_status_flip_frac": flip,
    }


def _job(job):
    condition, scen, rep = job
    cfg = _cfg(condition, scen, seed_for(f"aux|{rep}", rep))
    out = run_config(cfg)
    s = out.series()
    if condition == "revocation":
        sc = _score_two_stream(s, _rule(), cfg.n_cycles, cfg.revocation_cycle)
    else:
        sc = score_run(s, _rule(), cfg.n_cycles)
        sc["aggregate_invalidation_mean_abs_z"] = 0.0
        sc["aggregate_alarm_status_flip_frac"] = 0.0
    sc.update({"condition": condition, "scen": scen, "rep": rep, "seed": cfg.seed,
               "nominal_c": cfg.consent_c,
               "achieved_edge_frac": float(s["achieved_edge_frac"][-1]),
               "n_revoked": int(np.sum(out.consent_revoked)),
               "n_cycles": cfg.n_cycles,
               "obs_cycles": int(np.sum(~np.isnan(s["obs_formation"])))})
    return sc


# ------------------------------------------------------------------ Ensign arm
def _shrink(k: float, m: float) -> float:
    if k >= m or m <= 1:
        return 1.0
    return (k / (k + 1.0)) / (m / (m + 1.0))


def _p_obs(m: float, p: float) -> float:
    """P(at least 2 of m team members consent) at per-agent consent probability p."""
    q = 1.0 - p
    return max(1e-3, 1.0 - q ** m - m * p * q ** (m - 1.0))


def _ipw_correct(program_log, p_hat: float, n_cycles: int) -> np.ndarray:
    """The Ensign fix as an observer with only consent status could execute it.

    The observer sees, per observed program, the number of CONSENTED members k and the
    recorded latency; it does NOT see the team's non-consenting members, so it cannot
    know the team size m.  It must infer m from the estimated consent propensity,
    m_hat = k / p_hat, undo the partial-observation shrinkage at m_hat, and weight the
    program by the inverse probability that a team of m_hat would have been observed at
    all.  This is exactly where the propensity is load-bearing -- and exactly what
    revocation corrupts."""
    num = np.zeros(n_cycles)
    den = np.zeros(n_cycles)
    p = min(0.999, max(0.05, p_hat))
    for (t, _m_true, k, v, _lat) in program_log:
        if t < 0 or t >= n_cycles:
            continue
        m_hat = max(float(k), float(k) / p)
        x = v / _shrink(float(k), m_hat)
        w = 1.0 / _p_obs(m_hat, p)
        num[t] += w * x
        den[t] += w
    out = np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)
    return out


def _ipw_job(job):
    condition, scen, rep = job
    cfg = _cfg(condition, scen, seed_for(f"aux|{rep}", rep))
    out = run_config(cfg)
    s = out.series()

    consent = np.asarray(out.consent_status, float)       # the record AS IT STANDS
    obsdeg = np.asarray(out.observed_degree, float)
    sd = obsdeg.std()
    x = (obsdeg - obsdeg.mean()) / (sd if sd > 1e-9 else 1.0)
    X = np.column_stack([np.ones_like(x), x])
    w = logistic_fit(X, consent)
    pi = np.clip(logistic_predict(X, w), 0.02, 1.0)
    p_hat = float(np.mean(pi))            # estimated marginal consent rate
    cons = consent > 0.5
    pi_pair = float(np.mean(pi[cons])) if cons.any() else 1.0

    truth = s["true_formation"]
    obs = s["obs_formation"]
    obs_ipw = _ipw_correct(out.program_log, p_hat, cfg.n_cycles)
    m = np.isfinite(truth) & np.isfinite(obs)
    m2 = np.isfinite(truth) & np.isfinite(obs_ipw)
    bias_raw = float(np.mean(obs[m] - truth[m])) if m.any() else float("nan")
    bias_ipw = float(np.mean(obs_ipw[m2] - truth[m2])) if m2.any() else float("nan")

    sc_raw = score_run(s, _rule(), cfg.n_cycles)
    s_ipw = dict(s)
    s_ipw["obs_formation"] = obs_ipw
    sc_ipw = score_run(s_ipw, _rule(), cfg.n_cycles)

    revoked = np.asarray(out.consent_revoked, float)
    # the propensity model's labels: a revoked agent is indistinguishable from an
    # agent that never consented, so its true observation propensity (it WAS observed
    # for the first R cycles) is unrecoverable from the record.
    mislabelled = int(np.sum((revoked > 0.5) & (consent < 0.5)))
    return {"condition": condition, "rep": rep, "seed": cfg.seed,
            "bias_raw": bias_raw, "bias_ipw": bias_ipw,
            "abs_bias_raw": abs(bias_raw), "abs_bias_ipw": abs(bias_ipw),
            "bias_reduction": abs(bias_raw) - abs(bias_ipw),
            "propensity_slope": float(w[1]),
            "pi_pair": pi_pair, "p_hat": p_hat,
            "n_revoked": int(revoked.sum()),
            "mislabelled_in_propensity_fit": mislabelled,
            "dt_raw": sc_raw["dt_censored"], "dt_ipw": sc_ipw["dt_censored"],
            "lead_fired_raw": sc_raw["lead_fired"],
            "lead_fired_ipw": sc_ipw["lead_fired"]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--rule", default=None,
                    help="detection-rule file in configs/ (default detection_rule.json)")
    ap.add_argument("--out", default="aux_arms.json")
    args = ap.parse_args()
    if args.rule:
        os.environ["SIM_DETECTION_RULE"] = args.rule
    reps = args.reps

    jobs = [(cond, scen, r) for cond in CONDITIONS
            for scen in ("degrade", "control") for r in range(reps)]
    rows = pmap(_job, jobs, args.workers, "arm 1: revocation vs static-equivalent")

    n_cycles = rows[0]["n_cycles"]
    cells = []
    for cond in CONDITIONS:
        d = [r for r in rows if r["condition"] == cond and r["scen"] == "degrade"]
        k = [r for r in rows if r["condition"] == cond and r["scen"] == "control"]
        dt = bootstrap_mean([r["dt_censored"] for r in d])
        cells.append({
            "condition": cond, "n_reps": len(d), "nominal_c": d[0]["nominal_c"],
            "achieved_edge_frac_mean": float(np.mean(
                [r["achieved_edge_frac"] for r in d])),
            "mean_agents_revoked": float(np.mean([r["n_revoked"] for r in d])),
            "dt_censored_mean": dt["mean"], "dt_censored_se": dt["se"],
            "lead_detect_rate": float(np.mean([r["lead_fired"] for r in d])),
            "lead_median_t": (float(np.median([r["t_lead"] for r in d
                                               if r["lead_fired"]]))
                              if any(r["lead_fired"] for r in d) else None),
            "fpr_leading": float(np.mean([r["lead_fired"] for r in k])),
            "auc_discrimination": auc([r["lead_max_score"] for r in d],
                                      [r["lead_max_score"] for r in k]),
            "aggregate_invalidation_mean_abs_z": float(np.nanmean(
                [r["aggregate_invalidation_mean_abs_z"] for r in d])),
            "aggregate_alarm_status_flip_frac": float(np.nanmean(
                [r["aggregate_alarm_status_flip_frac"] for r in d])),
            "observed_cycle_frac": float(np.mean([r["obs_cycles"] for r in d])
                                         / n_cycles),
        })
    by = {x["condition"]: x for x in cells}

    rv = {r["rep"]: r for r in rows
          if r["condition"] == "revocation" and r["scen"] == "degrade"}
    se_ = {r["rep"]: r for r in rows
           if r["condition"] == "static_equiv" and r["scen"] == "degrade"}
    shared = sorted(set(rv) & set(se_))
    pb = bootstrap_mean([rv[i]["dt_censored"] - se_[i]["dt_censored"] for i in shared])
    contrast = {
        "timing_leg": {
            "paired_dt_difference_mean": pb["mean"], "se": pb["se"],
            "lo95": pb["lo95"], "hi95": pb["hi95"], "n_pairs": pb["n"],
            "revocation_detect_rate": by["revocation"]["lead_detect_rate"],
            "static_equiv_detect_rate": by["static_equiv"]["lead_detect_rate"],
            "revocation_fpr": by["revocation"]["fpr_leading"],
            "static_equiv_fpr": by["static_equiv"]["fpr_leading"],
            "achieved_coverage_matched": [
                by["revocation"]["achieved_edge_frac_mean"],
                by["static_equiv"]["achieved_edge_frac_mean"]],
            "distinguishable": bool(pb["lo95"] > 0 or pb["hi95"] < 0),
        },
        "semantics_leg_aggregate_invalidation": {
            "revocation_mean_abs_z_shift": by["revocation"][
                "aggregate_invalidation_mean_abs_z"],
            "revocation_alarm_status_flip_frac": by["revocation"][
                "aggregate_alarm_status_flip_frac"],
            "static_equiv_mean_abs_z_shift": by["static_equiv"][
                "aggregate_invalidation_mean_abs_z"],
            "note": ("static coverage reduction cannot have a non-zero value here by "
                     "construction: there is no already-published record to re-derive. "
                     "This is the property that distinguishes revocation from an "
                     "equivalent static reduction independently of detection timing."),
            "distinguishable": bool(
                by["revocation"]["aggregate_invalidation_mean_abs_z"] > 0.25),
        },
        "kill_criterion_1_note": (
            "if the timing leg is null, spec section 6 kill criterion 1 applies to the "
            "revocation-SEMANTICS sub-claim ONLY: retain C-II's ownership/capability "
            "leg and lead on C-III."),
    }

    ijobs = [(cond, "degrade", r) for cond in ("static_full", "revocation")
             for r in range(reps)]
    irows = pmap(_ipw_job, ijobs, args.workers, "arm 2: Ensign fix (IPW)")
    ipw = {}
    for cond in ("static_full", "revocation"):
        rs = [r for r in irows if r["condition"] == cond]
        br = bootstrap_mean([r["bias_reduction"] for r in rs])
        ipw[cond] = {
            "n_reps": len(rs),
            "abs_bias_raw": float(np.nanmean([r["abs_bias_raw"] for r in rs])),
            "abs_bias_ipw": float(np.nanmean([r["abs_bias_ipw"] for r in rs])),
            "bias_reduction_mean": br["mean"], "bias_reduction_se": br["se"],
            "bias_reduction_lo95": br["lo95"], "bias_reduction_hi95": br["hi95"],
            "frac_runs_bias_reduced": float(np.mean(
                [r["bias_reduction"] > 0 for r in rs])),
            "propensity_slope_mean": float(np.mean(
                [r["propensity_slope"] for r in rs])),
            "pi_pair_mean": float(np.mean([r["pi_pair"] for r in rs])),
            "p_hat_mean": float(np.mean([r["p_hat"] for r in rs])),
            "mean_agents_revoked": float(np.mean([r["n_revoked"] for r in rs])),
            "mislabelled_in_propensity_fit_mean": float(np.mean(
                [r["mislabelled_in_propensity_fit"] for r in rs])),
            "dt_raw_mean": float(np.mean([r["dt_raw"] for r in rs])),
            "dt_ipw_mean": float(np.mean([r["dt_ipw"] for r in rs])),
            "lead_detect_raw": float(np.mean([r["lead_fired_raw"] for r in rs])),
            "lead_detect_ipw": float(np.mean([r["lead_fired_ipw"] for r in rs])),
        }
    ipw["_where_it_fails"] = {
        "static_bias_reduction": ipw["static_full"]["bias_reduction_mean"],
        "revocation_bias_reduction": ipw["revocation"]["bias_reduction_mean"],
        "ratio_revocation_over_static": (
            ipw["revocation"]["bias_reduction_mean"]
            / ipw["static_full"]["bias_reduction_mean"]
            if abs(ipw["static_full"]["bias_reduction_mean"]) > 1e-9 else None),
        "correction_survives_revocation": bool(
            ipw["revocation"]["bias_reduction_mean"]
            > 0.5 * ipw["static_full"]["bias_reduction_mean"]),
        "mean_agents_mislabelled_under_revocation": ipw["revocation"][
            "mislabelled_in_propensity_fit_mean"],
        "mechanism": ("revocation erases the consent record retroactively, so revoked "
                      "agents enter the propensity fit labelled as never-consenting; "
                      "the estimator is fitted on corrupted labels and the propensity "
                      "it recovers is not the propensity that generated the observed "
                      "data.  The correction has no way to know it is wrong."),
    }

    write_agg(args.out, {
        "_experiment": "E-A auxiliary arms at (c=0.5, mid red noise) -- spec D16 bound",
        "design": {"c": C_AUX, "colour": COLOUR_AUX, "severity": SEV_AUX,
                   "conditions": list(CONDITIONS), "reps": reps,
                   "rule_file": os.environ.get("SIM_DETECTION_RULE",
                                               "detection_rule.json"),
                   "revocation_cycle": Config().revocation_cycle,
                   "revocation_frac": Config().revocation_frac,
                   "n_cycles": n_cycles,
                   "matched_seeds": "seed depends on replicate only; all three "
                                    "conditions share it",
                   "ensign_estimator": "logistic(consent ~ 1 + observed-subgraph "
                                       "degree), frozen"},
        "detection_rule": _rule(),
        "revocation_arm": {"cells": cells, "contrast": contrast},
        "ensign_arm": ipw,
        "per_run_rows": rows,
        "ensign_rows": irows,
    })

    print("\n  revocation arm:")
    for x in cells:
        print(f"    {x['condition']:13s} achieved={x['achieved_edge_frac_mean']:.3f} "
              f"revoked={x['mean_agents_revoked']:5.1f} "
              f"dt={x['dt_censored_mean']:+7.1f}+-{x['dt_censored_se']:4.1f} "
              f"DR={x['lead_detect_rate']:.2f} FPR={x['fpr_leading']:.2f} "
              f"AUC={x['auc_discrimination']:.2f} "
              f"inval={x['aggregate_invalidation_mean_abs_z']:.2f} "
              f"flip={x['aggregate_alarm_status_flip_frac']:.2f}")
    t = contrast["timing_leg"]
    print(f"    paired dt(revocation - static_equiv) = {t['paired_dt_difference_mean']:+.1f} "
          f"[{t['lo95']:+.1f},{t['hi95']:+.1f}] distinguishable={t['distinguishable']}")
    s_ = contrast["semantics_leg_aggregate_invalidation"]
    print(f"    aggregate invalidation |dz|={s_['revocation_mean_abs_z_shift']:.2f} "
          f"alarm-status flips={s_['revocation_alarm_status_flip_frac']:.1%} "
          f"distinguishable={s_['distinguishable']}")
    print("\n  Ensign fix:")
    for cond in ("static_full", "revocation"):
        v = ipw[cond]
        print(f"    {cond:13s} |bias| raw={v['abs_bias_raw']:.3f} "
              f"ipw={v['abs_bias_ipw']:.3f} reduction={v['bias_reduction_mean']:+.3f} "
              f"[{v['bias_reduction_lo95']:+.3f},{v['bias_reduction_hi95']:+.3f}] "
              f"improved in {v['frac_runs_bias_reduced']:.0%} of runs, "
              f"mislabelled={v['mislabelled_in_propensity_fit_mean']:.1f}")
    print(f"    survives revocation: "
          f"{ipw['_where_it_fails']['correction_survives_revocation']}")


if __name__ == "__main__":
    main()
