"""
coverage-and-noise campaign -- lead time under consent-bounded observation (frozen spec section 3).

Grid: consent coverage c x proxy-error colour x severity x consent arm, with a
PAIRED no-transition control run on MATCHED SEEDS in every cell (Arnold & Ferrari
2025 template; the warrant is Boettiger & Hastings 2012b -- conditioning on systems
known to have transitioned inflates apparent EWS performance).

Matched-seed construction: the seed depends on (c, consent arm, replicate) ONLY.
The degradation and control scenarios, and all nine proxy-error settings, therefore
share one substrate seed per replicate, so every contrast in this file is paired.

Reported per cell (never lead time alone):
  * Delta t, the lead-time advantage of the leading indicator over the lagging
    output metric, with a bootstrap standard error;
  * FPR at the frozen rule, measured on the paired control arm;
  * a discrimination statistic (AUC of the leading rule's score, transition vs
    matched control);
  * nominal c against the ACHIEVED observed-edge fraction (Cencetti et al.'s
    input/output split) -- their divergence under MNAR is itself a result.

The detection rule is loaded from configs/detection_rule.json and applied UNCHANGED
in every cell.  The c sweep is the EXPERIMENTER's cross-scenario comparison over
counterfactual consent regimes, never an operation available to the in-model
instrument within a run.

Run:  python run_ea_headline.py [--workers N] [--reps 30] [--quick]
      python run_ea_headline.py --rule detection_rule_recal.json \
                               --out ea_headline_recal.json \
                               --c-extra 0.333 --iso-c 0.25,0.333,0.5 --iso-reps 100
Writes: aggregates/ea_headline.json (or --out)

FOLLOW-UP CAMPAIGN 2026-08-05 (spec section 9, amendments large-null recalibration / coverage-grid refinement).  Three flags were
added; none changes the design:
  --rule      score under a different pre-declared rule file (the recalibrated
              threshold).  The runs themselves are identical: a seed depends on
              (c, arm, replicate) only, so re-running reproduces campaign 1's runs
              bit-for-bit and only the scoring changes.
  --iso-c / --iso-reps   raise replication on the cells adjacent to the Delta t = 0
              iso-line, at the headline severity only.  Adding replicates to an
              existing cell is not a design change.
  --c-extra   add a grid POINT in c (amendment coverage-grid refinement, coverage-and-noise campaign only, reason: iso-line
              localisation).  The added point is run at the headline severity only,
              since the iso-line is defined there.
"""

from __future__ import annotations

import argparse
import itertools
import os

import numpy as np

from model import Config, run_config
from scoring import score_run, bootstrap_mean
from indicators import auc
from runner import pmap, seed_for, write_agg, load_rule

C_GRID = (0.25, 0.5, 0.75, 1.0)
COLOURS = ("white", "red", "systematic")
SEVERITIES = ("low", "mid", "high")
ARMS = ("mnar", "mcar")           # mnar = primary (centrality-correlated), mcar = ref
SCENARIOS = ("degrade", "control")
HEADLINE_SEVERITY = "mid"

_RULE = None


def _rule():
    global _RULE
    if _RULE is None:
        _RULE = load_rule()
    return _RULE


def _job(job):
    c, arm, colour, sev, scen, rep = job
    cfg = Config(consent_c=c, consent_arm=arm, eps_colour=colour, eps_severity=sev,
                 scenario=scen, seed=seed_for(f"ea|{c}|{arm}", rep))
    out = run_config(cfg)
    s = out.series()
    sc = score_run(s, _rule(), cfg.n_cycles)
    sc.update({"c": c, "arm": arm, "colour": colour, "sev": sev, "scen": scen,
               "rep": rep, "seed": cfg.seed,
               "achieved_edge_frac": float(s["achieved_edge_frac"][-1]),
               "obs_cycles_with_data": int(np.sum(~np.isnan(s["obs_formation"]))),
               "n_cycles": cfg.n_cycles})
    return sc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--rule", default=None,
                    help="detection-rule file in configs/ (default detection_rule.json)")
    ap.add_argument("--out", default="ea_headline.json")
    ap.add_argument("--c-extra", default="",
                    help="comma-separated extra c grid points, headline severity only")
    ap.add_argument("--iso-c", default="",
                    help="comma-separated c values to raise replication on")
    ap.add_argument("--iso-reps", type=int, default=0)
    args = ap.parse_args()
    if args.rule:
        os.environ["SIM_DETECTION_RULE"] = args.rule
    reps = 4 if args.quick else args.reps
    colours = ("white", "red") if args.quick else COLOURS
    sevs = (HEADLINE_SEVERITY,) if args.quick else SEVERITIES

    extra = tuple(float(x) for x in args.c_extra.split(",") if x.strip())
    c_grid = tuple(sorted(set(C_GRID) | set(extra)))
    iso_c = tuple(float(x) for x in args.iso_c.split(",") if x.strip())
    iso_reps = args.iso_reps

    def reps_for(c: float, sev: str) -> int:
        """Replicate count per cell.  Extra grid points exist at the headline severity
        only; iso-line-adjacent cells are raised there too."""
        if c in extra and sev != HEADLINE_SEVERITY:
            return 0
        if iso_reps and c in iso_c and sev == HEADLINE_SEVERITY:
            return max(iso_reps, reps)
        return reps

    jobs = [(c, arm, col, sev, scen, r)
            for c, arm, col, sev, scen in itertools.product(
                c_grid, ARMS, colours, sevs, SCENARIOS)
            for r in range(reps_for(c, sev))]
    rows = pmap(_job, jobs, args.workers, "coverage-and-noise campaign grid")

    n_cycles = rows[0]["n_cycles"]
    cells = []
    for c, arm, col, sev in itertools.product(c_grid, ARMS, colours, sevs):
        d = [r for r in rows if r["c"] == c and r["arm"] == arm
             and r["colour"] == col and r["sev"] == sev and r["scen"] == "degrade"]
        k = [r for r in rows if r["c"] == c and r["arm"] == arm
             and r["colour"] == col and r["sev"] == sev and r["scen"] == "control"]
        if not d or not k:
            continue
        dt = bootstrap_mean([r["dt_censored"] for r in d])
        dtp = bootstrap_mean([r["dt_paired"] for r in d])
        fired_lead = [r for r in d if r["lead_fired"]]
        cells.append({
            "c": c, "arm": arm, "colour": col, "severity": sev,
            "n_reps": len(d),
            "dt_censored_mean": dt["mean"], "dt_censored_se": dt["se"],
            "dt_censored_lo95": dt["lo95"], "dt_censored_hi95": dt["hi95"],
            "dt_paired_mean": dtp["mean"], "dt_paired_se": dtp["se"],
            "dt_paired_n": dtp["n"],
            "lead_detect_rate": float(np.mean([r["lead_fired"] for r in d])),
            "lag_detect_rate": float(np.mean([r["lag_fired"] for r in d])),
            "lead_median_t": (float(np.median([r["t_lead"] for r in fired_lead]))
                              if fired_lead else None),
            "lag_median_t": (float(np.median([r["t_lag"] for r in d
                                              if r["lag_fired"]]))
                             if any(r["lag_fired"] for r in d) else None),
            "fpr_leading": float(np.mean([r["lead_fired"] for r in k])),
            "fpr_lagging": float(np.mean([r["lag_fired"] for r in k])),
            "auc_discrimination": auc([r["lead_max_score"] for r in d],
                                      [r["lead_max_score"] for r in k]),
            "nominal_c": c,
            "achieved_edge_frac_mean": float(np.mean(
                [r["achieved_edge_frac"] for r in d])),
            "achieved_over_nominal": float(np.mean(
                [r["achieved_edge_frac"] for r in d]) / c),
            "observed_cycle_frac": float(np.mean(
                [r["obs_cycles_with_data"] for r in d]) / n_cycles),
        })

    # ---- Delta t = 0 iso-line in c, per (arm, colour) at the headline severity --
    iso = {}
    for arm, col in itertools.product(ARMS, colours):
        cs = [x for x in cells if x["arm"] == arm and x["colour"] == col
              and x["severity"] == HEADLINE_SEVERITY]
        cs.sort(key=lambda x: x["c"])
        if len(cs) < 2:
            continue
        entry = {"c": [x["c"] for x in cs],
                 "dt": [x["dt_censored_mean"] for x in cs],
                 "dt_se": [x["dt_censored_se"] for x in cs],
                 "n_reps": [x["n_reps"] for x in cs]}
        cross = None
        for a, b in zip(cs[:-1], cs[1:]):
            if (a["dt_censored_mean"] <= 0.0) != (b["dt_censored_mean"] <= 0.0):
                y0, y1 = a["dt_censored_mean"], b["dt_censored_mean"]
                cross = a["c"] + (b["c"] - a["c"]) * (0.0 - y0) / (y1 - y0)
                break
        entry["c_at_dt_zero"] = cross
        if cross is None:
            entry["note"] = ("no sign change of Delta t across the swept c range: "
                             + ("Delta t > 0 at every c (coverage floor below 0.25)"
                                if cs[0]["dt_censored_mean"] > 0 else
                                "Delta t <= 0 at every c"))
        # band from seed variation: iso-line location under +-1 se on Delta t
        band = []
        for sgn in (-1.0, +1.0):
            vals = [x["dt_censored_mean"] + sgn * x["dt_censored_se"] for x in cs]
            b_cross = None
            for i in range(len(cs) - 1):
                if (vals[i] <= 0.0) != (vals[i + 1] <= 0.0):
                    b_cross = (cs[i]["c"] + (cs[i + 1]["c"] - cs[i]["c"])
                               * (0.0 - vals[i]) / (vals[i + 1] - vals[i]))
                    break
            band.append(b_cross)
        entry["c_at_dt_zero_band"] = band
        iso[f"{arm}|{col}"] = entry

    # ---- pooled false-alarm accounting against both bars ----------------------
    ctl_all = [r for r in rows if r["scen"] == "control"]
    mid_ctl = [r for r in ctl_all if r["sev"] == HEADLINE_SEVERITY]
    pooled = {
        "n_control_runs": len(ctl_all),
        "fpr_leading_pooled": float(np.mean([r["lead_fired"] for r in ctl_all])),
        "fpr_lagging_pooled": float(np.mean([r["lag_fired"] for r in ctl_all])),
        "declared_target": _rule().get("nominal_fpr_target"),
        "external_bar": 0.036,
        "fpr_leading_by_colour_mid_severity": {
            col: float(np.mean([r["lead_fired"] for r in mid_ctl
                                if r["colour"] == col])) for col in colours},
        "fpr_leading_by_c_mid_severity": {
            str(c): float(np.mean([r["lead_fired"] for r in mid_ctl
                                   if r["c"] == c])) for c in c_grid
            if any(r["c"] == c for r in mid_ctl)},
        "per_cell_fpr_min": float(min(x["fpr_leading"] for x in cells)),
        "per_cell_fpr_max": float(max(x["fpr_leading"] for x in cells)),
        "per_cell_fpr_median": float(np.median([x["fpr_leading"] for x in cells])),
        "meets_declared_target": bool(
            float(np.mean([r["lead_fired"] for r in ctl_all]))
            <= (_rule().get("nominal_fpr_target") or 0.05)),
        "meets_external_bar": bool(
            float(np.mean([r["lead_fired"] for r in ctl_all])) <= 0.036),
    }

    # ---- confusion time-series material (figure hysteresis analysis) --------------------------
    confusion = {}
    for arm, col in itertools.product(ARMS, colours):
        for c in c_grid:
            d = [r for r in rows if r["c"] == c and r["arm"] == arm
                 and r["colour"] == col and r["sev"] == HEADLINE_SEVERITY
                 and r["scen"] == "degrade"]
            k = [r for r in rows if r["c"] == c and r["arm"] == arm
                 and r["colour"] == col and r["sev"] == HEADLINE_SEVERITY
                 and r["scen"] == "control"]
            if not d:
                continue
            confusion[f"{arm}|{col}|{c}"] = {
                "alarm_cycles_transition": [r["t_lead"] for r in d],
                "alarm_cycles_control": [r["t_lead"] for r in k],
                "lag_alarm_cycles_transition": [r["t_lag"] for r in d],
            }

    write_agg(args.out, {
        "_experiment": "coverage-and-noise campaign lead time under consent-bounded observation",
        "design": {"c_grid": list(c_grid), "colours": list(colours),
                   "severities": list(sevs), "arms": list(ARMS),
                   "headline_severity": HEADLINE_SEVERITY, "reps": reps,
                   "n_cycles": n_cycles,
                   "paired_control": "matched seeds; seed depends on (c, arm, rep) only",
                   "external_fpr_bar": 0.036,
                   "c_grid_frozen": list(C_GRID),
                   "c_extra_points": list(extra),
                   "iso_line_c_values_raised": list(iso_c),
                   "iso_line_reps": iso_reps or reps,
                   "rule_file": os.environ.get("SIM_DETECTION_RULE",
                                               "detection_rule.json"),
                   "n_runs": len(rows)},
        "detection_rule": _rule(),
        "cells": cells,
        "pooled_false_alarms": pooled,
        "iso_line_dt_zero": iso,
        "confusion_series": confusion,
        "per_run_rows": rows,
    })

    print("\n  c    arm   colour      dt     se   leadDR  FPR   AUC  achieved/c")
    for x in sorted(cells, key=lambda y: (y["colour"], y["arm"], -y["c"])):
        if x["severity"] != HEADLINE_SEVERITY:
            continue
        print(f"  {x['c']:.2f} {x['arm']:5s} {x['colour']:10s} "
              f"{x['dt_censored_mean']:+7.1f} {x['dt_censored_se']:5.1f} "
              f"{x['lead_detect_rate']:5.2f} {x['fpr_leading']:5.2f} "
              f"{x['auc_discrimination']:5.2f} {x['achieved_over_nominal']:5.2f}")
    print("\n  iso-line c at Delta t = 0:")
    for k, v in iso.items():
        print(f"    {k}: {v['c_at_dt_zero']}  band={v['c_at_dt_zero_band']}")


if __name__ == "__main__":
    main()
