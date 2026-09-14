"""
Cross-campaign analysis of the recalibration (follow-up campaign 2026-08-05).

No new runs: everything here is derived from two committed aggregates,
`ea_headline.json` (campaign 1, threshold k = 5.0, lagging 3.5 sigma sustained 5) and
`ea_headline_recal.json` (recalibrated, k = 6.5, lagging 3.5 sigma sustained 3), on the
SAME seeds -- a seed depends on (c, arm, replicate) only, so the two campaigns scored the
identical runs and every contrast below is paired.

Three things it computes, all of which the follow-up report needs and none of which is a
single-aggregate lookup:

  1 DELTA-T DECOMPOSITION.  Two thresholds moved at once, and they push Delta t the same
    way, so the drop must be attributed:
        Delta t             = eff_lag - eff_lead        (per run, censored at run length)
        Delta t_new - old   = mean(eff_lag_new - eff_lag_old)      <- comparator effect
                            - mean(eff_lead_new - eff_lead_old)   <- leading-rule effect
    The comparator effect is the cost of honouring matched false-alarm calibration (both detectors calibrated to the
    same tightened target, comparator gets the earliest admissible rule); the leading
    effect is the cost of the stricter alarm threshold itself.  The identity is asserted.

  2 POOLED FALSE-ALARM RATES with exact binomial and Wilson 95% intervals, so "above the
    3.6% bar" is a statistical statement rather than a point comparison.

  3 A SEED-NOISE GAUGE ON FPR.  At c = 1.0 the MNAR and MCAR arms are the SAME physical
    configuration -- consent is universal, so the assignment rule is irrelevant -- run on
    disjoint seed sets.  Their FPR difference is therefore pure seed noise, which is how
    much a 30-run per-cell FPR can be trusted.  (RESULTS.md section 3.3 uses the same
    device for Delta t.)

Run:  python recal_analysis.py
Writes: aggregates/recal_analysis.json
"""

from __future__ import annotations

import json
import math
import os

import numpy as np

from runner import write_agg, AGG

OLD = "ea_headline.json"
NEW = "ea_headline_recal.json"
EB_N500 = "eb_joint_sweep_n500.json"
EXTERNAL_BAR = 0.036
Z = 1.959963985


def _load(name):
    with open(os.path.join(AGG, name)) as f:
        return json.load(f)


def wilson(k: int, n: int) -> list[float]:
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    c = (p + Z * Z / (2 * n)) / (1 + Z * Z / n)
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / (1 + Z * Z / n)
    return [c - h, c + h]


def _eff(r, key: str, n_cycles: int) -> float:
    v = r[key]
    return float(n_cycles if v is None else v)


def _fpr_block(d, tag: str) -> dict:
    rows = [r for r in d["per_run_rows"] if r["scen"] == "control"]
    n = len(rows)
    k = int(sum(1 for r in rows if r["lead_fired"]))
    kl = int(sum(1 for r in rows if r["lag_fired"]))
    lo, hi = wilson(k, n)
    return {
        "campaign": tag,
        "rule_k": d["detection_rule"]["leading"]["k"],
        "lagging_k_sd": d["detection_rule"]["lagging"]["k_sd"],
        "lagging_sustain": d["detection_rule"]["lagging"]["sustain"],
        "n_control_runs": n, "n_false_alarms": k,
        "fpr_leading": k / n, "fpr_leading_wilson95": [lo, hi],
        "fpr_lagging": kl / n,
        "declared_target": d["detection_rule"].get("nominal_fpr_target"),
        "external_bar": EXTERNAL_BAR,
        "above_external_bar_point": bool(k / n > EXTERNAL_BAR),
        "above_external_bar_ci_lower": bool(lo > EXTERNAL_BAR),
    }


def compute() -> dict:
    old, new = _load(OLD), _load(NEW)
    n_cycles = old["design"]["n_cycles"]
    assert new["design"]["n_cycles"] == n_cycles

    def index(d):
        out = {}
        for r in d["per_run_rows"]:
            out[(r["c"], r["arm"], r["colour"], r["sev"], r["scen"], r["rep"])] = r
        return out

    io, inew = index(old), index(new)
    shared = sorted(set(io) & set(inew))
    # the two campaigns must have scored the identical runs wherever they overlap
    mismatched = [k for k in shared if io[k]["seed"] != inew[k]["seed"]]

    cells = {}
    for key in shared:
        c, arm, col, sev, scen, rep = key
        if scen != "degrade":
            continue
        cells.setdefault((c, arm, col, sev), []).append(key)

    decomp = []
    for (c, arm, col, sev), keys in sorted(cells.items()):
        lo_ = [_eff(io[k], "t_lead", n_cycles) for k in keys]
        ln_ = [_eff(inew[k], "t_lead", n_cycles) for k in keys]
        go_ = [_eff(io[k], "t_lag", n_cycles) for k in keys]
        gn_ = [_eff(inew[k], "t_lag", n_cycles) for k in keys]
        dt_old = float(np.mean(np.asarray(go_) - np.asarray(lo_)))
        dt_new = float(np.mean(np.asarray(gn_) - np.asarray(ln_)))
        lag_eff = float(np.mean(np.asarray(gn_) - np.asarray(go_)))
        lead_eff = -float(np.mean(np.asarray(ln_) - np.asarray(lo_)))
        decomp.append({
            "c": c, "arm": arm, "colour": col, "severity": sev,
            "n_shared_reps": len(keys),
            "dt_campaign1": dt_old, "dt_recalibrated": dt_new,
            "dt_change": dt_new - dt_old,
            "attributable_to_comparator": lag_eff,
            "attributable_to_leading_threshold": lead_eff,
            "identity_residual": (dt_new - dt_old) - (lag_eff + lead_eff),
            "mean_lead_alarm_cycle_campaign1": float(np.mean(lo_)),
            "mean_lead_alarm_cycle_recalibrated": float(np.mean(ln_)),
            "mean_lag_alarm_cycle_campaign1": float(np.mean(go_)),
            "mean_lag_alarm_cycle_recalibrated": float(np.mean(gn_)),
        })

    hs = new["design"]["headline_severity"]
    head = [x for x in decomp if x["arm"] == "mnar" and x["severity"] == hs]
    pooled_shift = {
        "n_cells": len(decomp),
        "mean_dt_change": float(np.mean([x["dt_change"] for x in decomp])),
        "mean_attributable_to_comparator": float(np.mean(
            [x["attributable_to_comparator"] for x in decomp])),
        "mean_attributable_to_leading_threshold": float(np.mean(
            [x["attributable_to_leading_threshold"] for x in decomp])),
        "max_abs_identity_residual": float(np.max(
            [abs(x["identity_residual"]) for x in decomp])),
        "comparator_share_of_change": (
            float(np.mean([x["attributable_to_comparator"] for x in decomp])
                  / np.mean([x["dt_change"] for x in decomp]))),
    }

    # seed-noise gauge: at c = 1.0 the two arms are the same physical configuration
    gauge = {}
    for camp, d in (("campaign1", old), ("recalibrated", new)):
        per = {}
        for col in new["design"]["colours"]:
            vals = {}
            for arm in ("mnar", "mcar"):
                ctl = [r for r in d["per_run_rows"]
                       if r["scen"] == "control" and r["c"] == 1.0
                       and r["arm"] == arm and r["colour"] == col and r["sev"] == hs]
                vals[arm] = (float(np.mean([r["lead_fired"] for r in ctl]))
                             if ctl else float("nan"))
            per[col] = {"mnar_fpr": vals["mnar"], "mcar_fpr": vals["mcar"],
                        "abs_difference": abs(vals["mnar"] - vals["mcar"])}
        gauge[camp] = {
            "per_colour": per,
            "max_abs_difference": float(max(v["abs_difference"]
                                            for v in per.values())),
            "note": ("at c = 1.0 consent is universal, so MNAR and MCAR are the same "
                     "physical configuration on disjoint seed sets: their FPR "
                     "difference is pure seed noise on 30 control runs"),
        }

    # per-coverage false-alarm rates WITH intervals: the pooled excess in section 8.3 is
    # only interpretable per c if the per-c sample size can resolve a 3.5% rate, and at
    # 30 reps/cell it cannot.  This block says exactly which levels are resolvable.
    by_c = {}
    mid_ctl = [r for r in new["per_run_rows"]
               if r["scen"] == "control" and r["sev"] == hs]
    for c in sorted({r["c"] for r in mid_ctl}):
        sel = [r for r in mid_ctl if r["c"] == c]
        n = len(sel)
        k = int(sum(1 for r in sel if r["lead_fired"]))
        lo, hi = wilson(k, n)
        by_c[str(c)] = {"n_control_runs": n, "n_false_alarms": k, "fpr": k / n,
                        "wilson95": [lo, hi],
                        "resolvably_above_external_bar": bool(lo > EXTERNAL_BAR)}

    # coverage-and-coupling campaign, derived: the reported cell statistic is the mean of the SIGNED correlation,
    # while the paired coupling-loss statistic is a difference of |corr| per run.  The two
    # differ because a fraction of runs ANTI-TRACK (positive correlation: the published
    # indicator moves the wrong way against the latent state), and that fraction is itself
    # a reportable cost of restricted coverage.  Both are recomputed here from the
    # committed per-run rows of the >= 500-pair aggregate.
    eb_anti = {}
    try:
        eb = _load(EB_N500)
    except FileNotFoundError:
        eb = None
    if eb is not None:
        for cell in eb["cells"]:
            sel = [r["corr_indicator_latent"] for r in eb["per_run_rows"]
                   if r["c"] == cell["c"] and r["g_name"] == cell["g_name"]
                   and r["mode"] == cell["mode"]]
            a = np.asarray([x for x in sel
                            if x is not None and np.isfinite(x)], float)
            eb_anti[f"{cell['c']}|{cell['g_name']}|{cell['mode']}"] = {
                "n": int(a.size),
                "mean_signed_corr": float(a.mean()),
                "mean_abs_corr": float(np.abs(a).mean()),
                "frac_anti_tracking": float(np.mean(a > 0)),
            }

    # coverage-and-coupling campaign, derived: a BOUND on the (c x g) interaction.  The interaction is the change in
    # the paired coupling loss between the lowest and the highest coverage; its standard
    # error is the root-sum-square of the two per-cell bootstrap s.e.s, so the smallest
    # interaction these runs could have detected at 80% power is 2.8016 x that.  Reporting
    # this converts "not resolved" into "bounded below X", which is a result.
    eb_bound = {}
    if eb is not None:
        for g, v in eb["interaction"].items():
            cs = sorted(float(k) for k in v["per_c"])
            lo_c, hi_c = str(cs[0]), str(cs[-1])
            a, b = v["per_c"][lo_c], v["per_c"][hi_c]
            rng = a["tracking_loss_mean"] - b["tracking_loss_mean"]
            se = math.sqrt(a["se"] ** 2 + b["se"] ** 2)
            mde = (Z + 0.841621234) * se
            eb_bound[g] = {
                "c_low": cs[0], "c_high": cs[-1],
                "loss_at_low_c": a["tracking_loss_mean"],
                "loss_at_high_c": b["tracking_loss_mean"],
                "interaction_range": rng,
                "range_se": se,
                "range_ci95": [rng - Z * se, rng + Z * se],
                "mde_80pct_power_for_range": mde,
                "range_below_mde": bool(abs(rng) < mde),
            }

    iso = {}
    for key, e in new["iso_line_dt_zero"].items():
        o = old["iso_line_dt_zero"].get(key, {})
        iso[key] = {"campaign1_c_at_dt_zero": o.get("c_at_dt_zero"),
                    "recalibrated_c_at_dt_zero": e.get("c_at_dt_zero"),
                    "recalibrated_band": e.get("c_at_dt_zero_band"),
                    "recalibrated_n_reps": e.get("n_reps")}

    return {
        "_experiment": ("cross-campaign analysis of the threshold recalibration "
                        "(no new runs; derived from two committed aggregates)"),
        "sources": {"campaign1": OLD, "recalibrated": NEW,
                    "eb_n500": EB_N500},
        "paired_on": ("identical seeds: a seed depends on (c, arm, replicate) only, so "
                      "both campaigns scored the same runs"),
        "n_shared_runs": len(shared),
        "n_seed_mismatches": len(mismatched),
        "delta_t_decomposition": decomp,
        "delta_t_decomposition_pooled": pooled_shift,
        "headline_mnar_mid": [
            {"colour": x["colour"], "c": x["c"], "dt_campaign1": x["dt_campaign1"],
             "dt_recalibrated": x["dt_recalibrated"],
             "attributable_to_comparator": x["attributable_to_comparator"],
             "attributable_to_leading_threshold":
                 x["attributable_to_leading_threshold"]}
            for x in sorted(head, key=lambda y: (y["colour"], -y["c"]))],
        "pooled_false_alarm_rates": [_fpr_block(old, "campaign1"),
                                     _fpr_block(new, "recalibrated")],
        "fpr_seed_noise_gauge": gauge,
        "fpr_by_coverage_recalibrated_mid_severity": by_c,
        "eb_signed_vs_abs_and_antitracking": eb_anti,
        "eb_interaction_bound": eb_bound,
        "iso_line_shift": iso,
    }


def main() -> None:
    payload = compute()
    write_agg("recal_analysis.json", payload)
    p = payload["delta_t_decomposition_pooled"]
    print(f"  cells decomposed: {p['n_cells']}, identity residual max "
          f"{p['max_abs_identity_residual']:.2e}")
    print(f"  mean Delta t change {p['mean_dt_change']:+.2f} cycles = "
          f"{p['mean_attributable_to_comparator']:+.2f} (comparator) "
          f"{p['mean_attributable_to_leading_threshold']:+.2f} (leading threshold)")
    for b in payload["pooled_false_alarm_rates"]:
        print(f"  {b['campaign']:13s} FPR {b['fpr_leading']:.4f} "
              f"[{b['fpr_leading_wilson95'][0]:.4f},"
              f"{b['fpr_leading_wilson95'][1]:.4f}] on {b['n_control_runs']} controls "
              f"-- above the 3.6% bar (CI lower): {b['above_external_bar_ci_lower']}")
    print("  FPR seed-noise gauge (c=1.0, MNAR vs MCAR, same configuration): "
          f"campaign1 max {payload['fpr_seed_noise_gauge']['campaign1']['max_abs_difference']:.3f}, "
          f"recalibrated max {payload['fpr_seed_noise_gauge']['recalibrated']['max_abs_difference']:.3f}")


if __name__ == "__main__":
    main()
