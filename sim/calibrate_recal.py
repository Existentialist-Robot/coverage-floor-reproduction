"""
Threshold RE-calibration on a LARGE dedicated null (follow-up campaign 2026-08-05).

Why this exists (RESULTS.md section 7 item 2 -- the honest failure of campaign 1):
the frozen threshold k = 5.0 was the smallest admissible k on only 40 calibration
control runs, and that overfits the calibration sample -- the pooled false-alarm rate
measured on the 2 160 EVALUATION control runs came out at 8.98%, above its own 5%
nominal and 2.5x Falmagne, Stephenson & Levin (2026)'s 3.6% external bar.

What is recalibrated, and what is NOT (spec section 9, amendment A-8):
  * NOT the rule FORM.  Statistics, window, signal class, class mixing, combination
    and sustain length stay exactly as frozen in configs/detection_rule.json --
    rolling variance + rolling lag-1 autocorrelation of the in-window linearly
    detrended observed formation-latency series, window 30, standardised against the
    run's own reference period, combined as the mean of the two z-scores, sustained
    5 cycles.  Only the THRESHOLD is re-picked.  Re-picking the form would need the
    degradation arm and would no longer be the same procedure.
  * The SELECTION PROCEDURE is the frozen one with a tightened target, declared
    before any number below was looked at:
        leading  -- the SMALLEST k on the k grid whose calibration null FPR <= 3.5%;
        lagging  -- among (k_sd, sustain) whose calibration null FPR <= 3.5% and whose
                    calibration detection rate >= 0.5, the EARLIEST-detecting one,
                    i.e. the comparator keeps every advantage it had in campaign 1.
    The 3.5% target sits just under the 3.6% external bar so that the operating point
    is admissible against the bar by construction rather than by luck.
  * The NULL is dedicated, large and fresh: N_NULL >= 600 control runs on seed block
    950000+, disjoint from the 900000-900039 block used by campaign 1's calibration
    and from every scored seed (all scored seeds are < 900000 by construction --
    runner.SCORED_SEED_CEILING).  Calibration stays at c = 1.0, white / mid, MNAR, the
    same configuration campaign 1 calibrated on: the instrument may not re-tune itself
    per consent regime.

Exactness note: FPR curves are computed from a per-run sufficient statistic rather
than by re-running the alarm loop per candidate threshold.  For a sustain-s level
rule the alarm exists iff  max_i min(score[i..i+s-1]) >= k  over admissible windows,
so the per-run scalar S_s = max_i min(...) gives the whole FPR(k) curve exactly, and
lets the verifier recompute every curve from the committed per-run rows.

Run:  python calibrate_recal.py [--workers N] [--n-null 600]
Writes: configs/detection_rule_recal.json, aggregates/calibration_null_large.json
"""

from __future__ import annotations

import argparse
import itertools
import json
import os

import numpy as np

from model import Config, run_config
from indicators import baseline_z, alarm_from_z, lagging_alarm, fill_gaps
from runner import pmap, write_agg, load_rule, HERE, SCORED_SEED_CEILING

NULL_SEED_BASE = 950_000          # fresh block, disjoint from campaign 1's 900000+
CAL_DEGRADE_BASE = 900_000        # campaign 1's calibration degradation seeds, reused
N_CAL_DEGRADE = 40
TARGET_FPR = 0.035                # tightened target, just under the 3.6% external bar
EXTERNAL_BAR = 0.036
OLD_NOMINAL = 0.05
K_GRID = [round(2.0 + 0.25 * i, 3) for i in range(0, 57)]        # 2.00 .. 16.00
KSD_GRID = [round(1.0 + 0.25 * i, 3) for i in range(0, 33)]      # 1.00 .. 9.00
SUSTAIN_GRID = (3, 5)


def _cfg(scen: str, seed: int) -> Config:
    return Config(scenario=scen, seed=seed, consent_c=1.0, consent_arm="mnar",
                  eps_colour="white", eps_severity="mid")


def max_sustained_min(x: np.ndarray, sustain: int, start: int) -> float:
    """max over windows [i, i+sustain) with i >= start of min(x) inside the window,
    NaN treated as -inf.  A sustain-s level rule fires iff this statistic >= k."""
    y = np.asarray(x, float).copy()
    y[~np.isfinite(y)] = -np.inf
    n = y.size
    best = -np.inf
    for i in range(max(0, start), n - sustain + 1):
        m = y[i:i + sustain].min()
        if m > best:
            best = m
    return float(best)


def _job(job):
    """One calibration run reduced to the per-run sufficient statistics."""
    scen, seed = job
    cfg = _cfg(scen, seed)
    s = run_config(cfg).series()
    frozen = load_rule("detection_rule.json")
    L = frozen["leading"]
    z = baseline_z(s, int(L["window"]), tuple(L["classes"]), int(L["start"]),
                   float(L["min_frac"]), int(L["max_gap"]))
    score = 0.5 * (z["z_var"] + z["z_ac"])          # combine = "mean", as frozen
    lead_stat = max_sustained_min(score, int(L["sustain"]), int(L["start"]))
    # lagging: standardised drop below the baseline mean, in baseline sigmas
    G = frozen["lagging"]
    base_end = int(G["baseline_end"])
    y = fill_gaps(np.asarray(s["output_close"], float), int(G["max_gap"]))
    base = y[:base_end]
    base = base[~np.isnan(base)]
    mu, sd = float(base.mean()), float(base.std(ddof=1))
    drop = (mu - y) / (sd if sd > 1e-12 else 1e-12)
    row = {"scen": scen, "seed": int(seed),
           "lead_sustained_score": lead_stat,
           "lead_max_score": float(np.nanmax(score[int(L["start"]):])),
           "n_cycles": cfg.n_cycles}
    for su in SUSTAIN_GRID:
        row[f"lag_sustained_drop_s{su}"] = max_sustained_min(drop, su, base_end)
    if scen == "degrade":
        # detection cycle of the leading rule at each candidate k is needed only for
        # the lagging comparator's median-detection tie-break, so record the lagging
        # alarm cycle per (k_sd, sustain) instead: cheap, exact, and recomputable.
        det = {}
        for k_sd, su in itertools.product(KSD_GRID, SUSTAIN_GRID):
            r = lagging_alarm(s["output_close"],
                              {**G, "k_sd": k_sd, "sustain": su})
            det[f"{k_sd}|{su}"] = r["alarm_cycle"]
        row["lag_alarm_cycle_by_rule"] = det
        row["lead_alarm_cycle_frozen_k"] = alarm_from_z(
            z["z_var"], z["z_ac"], float(L["k"]), int(L["sustain"]),
            L["combine"], int(L["start"]))["alarm_cycle"]
    return row


def fpr_curve(stats: list[float], grid: list[float]) -> list[dict]:
    a = np.asarray(stats, float)
    return [{"k": k, "fpr": float(np.mean(a >= k)),
             "n_false_alarms": int(np.sum(a >= k))} for k in grid]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--n-null", type=int, default=600)
    args = ap.parse_args()
    n_null = args.n_null
    if n_null < 600:
        raise SystemExit("the follow-up campaign declares a null of at least 600 runs")

    frozen = load_rule("detection_rule.json")
    L = dict(frozen["leading"])
    G = dict(frozen["lagging"])

    jobs = [("control", NULL_SEED_BASE + r) for r in range(n_null)]
    jobs += [("degrade", CAL_DEGRADE_BASE + r) for r in range(N_CAL_DEGRADE)]
    rows = pmap(_job, jobs, args.workers, "calibration null (control) + cal degrade")
    null = [r for r in rows if r["scen"] == "control"]
    deg = [r for r in rows if r["scen"] == "degrade"]
    assert len(null) == n_null and len(deg) == N_CAL_DEGRADE

    # ------------------------------------------------------------ leading threshold
    lead_curve = fpr_curve([r["lead_sustained_score"] for r in null], K_GRID)
    adm = [x for x in lead_curve if x["fpr"] <= TARGET_FPR]
    if not adm:
        raise SystemExit("no k on the grid reaches the 3.5% target -- stop and report")
    k_new = min(x["k"] for x in adm)
    k_row = next(x for x in lead_curve if x["k"] == k_new)
    old_on_new_null = next((x for x in lead_curve
                            if abs(x["k"] - float(L["k"])) < 1e-9), None)

    # ------------------------------------------------------------ lagging threshold
    lag_rows = []
    for k_sd, su in itertools.product(KSD_GRID, SUSTAIN_GRID):
        fpr = float(np.mean([r[f"lag_sustained_drop_s{su}"] >= k_sd for r in null]))
        hits = [r["lag_alarm_cycle_by_rule"][f"{k_sd}|{su}"] for r in deg]
        fired = [h for h in hits if h is not None]
        lag_rows.append({"k_sd": k_sd, "sustain": su, "fpr": fpr,
                         "n_false_alarms": int(round(fpr * n_null)),
                         "detect_rate": len(fired) / len(deg),
                         "median_t": float(np.median(fired)) if fired else None})
    lag_adm = [r for r in lag_rows if r["fpr"] <= TARGET_FPR
               and r["detect_rate"] >= 0.5 and r["median_t"] is not None]
    if not lag_adm:
        raise SystemExit("no admissible lagging comparator at the 3.5% target")
    lag_adm.sort(key=lambda r: (r["median_t"], r["k_sd"]))
    lag_new = lag_adm[0]

    # ------------------------------------------------------------ the new rule file
    rule = {
        "_declared": (
            "RECALIBRATED 2026-08-05 on a dedicated null of "
            f"{n_null} control runs, seed block {NULL_SEED_BASE}+ (fresh, disjoint "
            "from campaign 1's calibration block 900000-900039 and from every scored "
            "seed, all of which are < 900000). Configuration c=1.0, white/mid, MNAR -- "
            "identical to campaign 1's calibration configuration; only the sample size "
            "and the FPR target changed. Rule FORM is unchanged from "
            "configs/detection_rule.json: only the leading threshold k and the lagging "
            "threshold k_sd were re-picked. Selection, declared in advance: leading = "
            "smallest k with calibration null FPR <= 3.5%; lagging = among rules with "
            "calibration null FPR <= 3.5% and calibration detection rate >= 0.5, the "
            "earliest-detecting. Applied UNCHANGED at every coverage level, noise "
            "colour, severity and consent arm. See spec section 9 amendment A-8."),
        "supersedes": "configs/detection_rule.json",
        "nominal_fpr_target": TARGET_FPR,
        "external_fpr_bar": EXTERNAL_BAR,
        "n_calibration_seeds": n_null,
        "calibration_seed_base": NULL_SEED_BASE,
        "scored_seed_ceiling": SCORED_SEED_CEILING,
        "leading": {**L, "k": k_new},
        "lagging": {**G, "k_sd": lag_new["k_sd"], "sustain": lag_new["sustain"]},
        "calibration_result": {
            "leading": {"k": k_new, "calibration_fpr": k_row["fpr"],
                        "n_false_alarms": k_row["n_false_alarms"],
                        "n_null": n_null,
                        "previous_k": float(L["k"]),
                        "previous_k_fpr_on_this_null": (
                            old_on_new_null["fpr"] if old_on_new_null else None)},
            "lagging": lag_new},
    }
    with open(os.path.join(HERE, "configs", "detection_rule_recal.json"), "w") as f:
        json.dump(rule, f, indent=1)

    write_agg("calibration_null_large.json", {
        "_experiment": ("threshold recalibration on a large dedicated null "
                        "(follow-up campaign 2026-08-05)"),
        "design": {
            "n_null_control_runs": n_null,
            "null_seed_base": NULL_SEED_BASE,
            "n_calibration_degrade_runs": N_CAL_DEGRADE,
            "calibration_degrade_seed_base": CAL_DEGRADE_BASE,
            "config": {"c": 1.0, "arm": "mnar", "colour": "white",
                       "severity": "mid", "scenario": "control"},
            "target_fpr": TARGET_FPR, "external_bar": EXTERNAL_BAR,
            "previous_nominal_fpr": OLD_NOMINAL,
            "k_grid": K_GRID, "k_sd_grid": KSD_GRID,
            "sustain_grid": list(SUSTAIN_GRID),
            "rule_form_unchanged": True,
            "what_was_repicked": ["leading.k", "lagging.k_sd", "lagging.sustain"],
            "selection_rule": ("leading: smallest k with null FPR <= 3.5%. lagging: "
                               "among rules with null FPR <= 3.5% and calibration "
                               "detection rate >= 0.5, the earliest-detecting."),
        },
        "frozen_rule_campaign1": frozen,
        "recalibrated_rule": rule,
        "leading_fpr_curve": lead_curve,
        "lagging_grid": lag_rows,
        "null_score_quantiles": {
            q: float(np.quantile([r["lead_sustained_score"] for r in null], q))
            for q in (0.5, 0.9, 0.95, 0.96, 0.965, 0.99, 1.0)},
        "chosen": {"leading_k": k_new, "leading_calibration_fpr": k_row["fpr"],
                   "lagging_k_sd": lag_new["k_sd"],
                   "lagging_sustain": lag_new["sustain"],
                   "lagging_calibration_fpr": lag_new["fpr"]},
        "campaign1_threshold_on_this_null": {
            "k": float(L["k"]),
            "fpr": old_on_new_null["fpr"] if old_on_new_null else None,
            "note": ("the campaign-1 threshold's false-alarm rate measured on 600 "
                     "fresh nulls instead of the 40 it was picked on"),
        },
        # committed in full so the verifier can recompute both selections from the rows:
        # the null rows carry the sufficient statistics for every FPR curve, and the
        # calibration degradation rows carry the lagging alarm cycle per candidate rule.
        "per_run_rows": rows,
    })

    print(f"\n  leading : k {L['k']} -> {k_new}  "
          f"(null FPR {old_on_new_null['fpr'] if old_on_new_null else float('nan'):.4f}"
          f" -> {k_row['fpr']:.4f} on {n_null} runs)")
    print(f"  lagging : k_sd {G['k_sd']}/s{G['sustain']} -> "
          f"{lag_new['k_sd']}/s{lag_new['sustain']}  "
          f"(null FPR {lag_new['fpr']:.4f}, cal median_t {lag_new['median_t']})")
    print("  wrote configs/detection_rule_recal.json")


if __name__ == "__main__":
    main()
