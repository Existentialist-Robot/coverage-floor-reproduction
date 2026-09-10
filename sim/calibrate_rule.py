"""
Detection-rule calibration on a DECLARED-DISJOINT seed block.

Pre-registration discipline (frozen specification, section 3):
  * The detection rule is chosen HERE, on CALIBRATION seeds >= 900000, and written
    to configs/detection_rule.json.  Every scored run in the E-A / auxiliary / E-B
    grids uses seeds < 900000, so no scored seed informs the rule.
  * Calibration happens at c = 1.0 only.  The rule is then applied UNCHANGED at
    every coverage level: the instrument cannot re-tune itself per consent regime.
  * Both detectors are calibrated against the SAME nominal false-alarm rate on the
    paired no-transition control arm, so lead time is compared at equal FPR
    (Boettiger & Hastings 2012b; Arnold & Ferrari 2025).
  * SELECTION CRITERION, fixed before inspecting any number: among variants whose
    calibration control FPR is <= NOMINAL_FPR and whose calibration detection rate
    is >= 0.5, take the one with the largest median lead time; ties on AUC.

  * AMENDMENT 1 (see spec section 9): the frozen spec said "pre-declared
    significance rule on the indicator trend".  A sequential Kendall-tau trend rule
    was implemented and tested first and CANNOT control the false-alarm rate on this
    substrate -- detection rate tracked FPR at every threshold (grid preserved in
    aggregates/calibration_rule_grids.json, key "rejected_kendall_tau_grid"),
    because overlapping rolling windows give tau a very wide null distribution and
    the first-crossing-over-the-run statistic is a maximum over many correlated
    draws.  The frozen rule is therefore a BASELINE-REFERENCED LEVEL rule on the
    same two variance-family statistics.  The rejected grid is a reportable result:
    it is the in-paper false-alarm critique measured rather than asserted.

Run:  python calibrate_rule.py
Writes: configs/detection_rule.json, aggregates/calibration_rule_grids.json
"""

from __future__ import annotations

import itertools
import json
import os

import numpy as np

from model import Config, run_config
from indicators import (baseline_z, alarm_from_z, class_trends,
                        alarm_from_trends, lagging_alarm, auc)

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_calcache.npz")
CAL_SEED_BASE = 900_000
N_CAL = 40
NOMINAL_FPR = 0.05
CLASSES = ("obs_formation", "obs_escalation")
SERIES_KEYS = ("obs_formation", "obs_escalation", "output_close",
               "true_formation", "true_escalation", "a_t", "giant_frac",
               "mean_degree")


def load_or_run() -> dict[str, list[dict[str, np.ndarray]]]:
    if os.path.exists(CACHE):
        d = np.load(CACHE)
        return {scen: [{k: d[f"{scen}/{r}/{k}"] for k in SERIES_KEYS}
                       for r in range(N_CAL)]
                for scen in ("degrade", "control")}
    store: dict[str, list] = {}
    flat: dict[str, np.ndarray] = {}
    for scen in ("degrade", "control"):
        rows = []
        for r in range(N_CAL):
            s = run_config(Config(scenario=scen, seed=CAL_SEED_BASE + r,
                                  consent_c=1.0, consent_arm="mnar",
                                  eps_colour="white", eps_severity="mid")).series()
            rows.append(s)
            for k in SERIES_KEYS:
                flat[f"{scen}/{r}/{k}"] = s[k]
        store[scen] = rows
        print(f"  calibration runs: {scen} done", flush=True)
    np.savez_compressed(CACHE, **flat)
    return store


# ---------------------------------------------------------------- lagging rule
def calibrate_lagging(deg, ctl, baseline_end: int) -> dict:
    rows = []
    for k, sustain in itertools.product((1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0),
                                        (3, 5)):
        rule = {"baseline_end": baseline_end, "k_sd": k, "sustain": sustain,
                "max_gap": 3}
        fp = sum(1 for s in ctl
                 if lagging_alarm(s["output_close"], rule)["alarm_cycle"] is not None)
        det = [h for h in (lagging_alarm(s["output_close"], rule)["alarm_cycle"]
                           for s in deg) if h is not None]
        rows.append({"k_sd": k, "sustain": sustain, "fpr": fp / len(ctl),
                     "detect_rate": len(det) / len(deg),
                     "median_t": float(np.median(det)) if det else None})
    ok = [r for r in rows if r["fpr"] <= NOMINAL_FPR and r["detect_rate"] >= 0.5]
    # the most generous admissible comparator: earliest median detection
    ok.sort(key=lambda r: r["median_t"])
    return {"grid": rows, "chosen": ok[0] if ok else None}


# ---------------------------------------------------------------- leading rule
def calibrate_leading(deg, ctl, baseline_end: int, lag_median: float) -> dict:
    rows = []
    cls_opts = (CLASSES, ("obs_formation",))
    for window, classes in itertools.product((30, 40, 50), cls_opts):
        dz = [baseline_z(s, window, classes, baseline_end) for s in deg]
        cz = [baseline_z(s, window, classes, baseline_end) for s in ctl]
        for sustain, combine in itertools.product((3, 5), ("and", "mean", "var_only")):
            for k in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0):
                cr = [alarm_from_z(t["z_var"], t["z_ac"], k, sustain, combine,
                                   baseline_end) for t in cz]
                fpr = sum(1 for r in cr if r["alarm_cycle"] is not None) / len(cr)
                if fpr > NOMINAL_FPR:
                    continue
                dr = [alarm_from_z(t["z_var"], t["z_ac"], k, sustain, combine,
                                   baseline_end) for t in dz]
                det = [r["alarm_cycle"] for r in dr if r["alarm_cycle"] is not None]
                rows.append({
                    "window": window, "classes": list(classes),
                    "sustain": sustain, "combine": combine, "k": k,
                    "fpr": fpr, "detect_rate": len(det) / len(dr),
                    "median_t": float(np.median(det)) if det else None,
                    "auc": auc([r["max_score"] for r in dr],
                               [r["max_score"] for r in cr])})
                break        # lowest admissible k for this variant
    adm = [r for r in rows if r["detect_rate"] >= 0.5 and r["median_t"] is not None]
    for r in adm:
        r["median_lead"] = lag_median - r["median_t"]
    adm.sort(key=lambda r: (-r["median_lead"], -r["auc"]))
    return {"grid": rows, "chosen": adm[0] if adm else None}


# ------------------------------------------- the rejected Kendall-tau rule
def kendall_tau_grid(deg, ctl) -> list[dict]:
    rows = []
    for window, span in ((30, 25), (40, 25), (40, 40), (50, 40)):
        dtr = [class_trends(s, window, span, CLASSES) for s in deg]
        ctr = [class_trends(s, window, span, CLASSES) for s in ctl]
        for sustain, combine in itertools.product((3, 5), ("and", "var_only")):
            for thr in (0.3, 0.4, 0.5, 0.6, 0.7):
                cr = [alarm_from_trends(t["tau_var"], t["tau_ac"], thr, sustain,
                                        combine) for t in ctr]
                dr = [alarm_from_trends(t["tau_var"], t["tau_ac"], thr, sustain,
                                        combine) for t in dtr]
                det = [r["alarm_cycle"] for r in dr if r["alarm_cycle"] is not None]
                rows.append({
                    "window": window, "tau_span": span, "sustain": sustain,
                    "combine": combine, "tau_threshold": thr,
                    "fpr": sum(1 for r in cr
                               if r["alarm_cycle"] is not None) / len(cr),
                    "detect_rate": len(det) / len(dr),
                    "median_t": float(np.median(det)) if det else None})
    return rows


def main() -> None:
    print(f"calibration: {2 * N_CAL} runs at c=1.0, white/mid", flush=True)
    runs = load_or_run()
    deg, ctl = runs["degrade"], runs["control"]
    base = Config()
    baseline_end = base.degrade_start + 40      # declared reference period

    lag = calibrate_lagging(deg, ctl, baseline_end)
    if lag["chosen"] is None:
        raise SystemExit("no admissible lagging comparator at nominal FPR")
    print("  lagging  :", lag["chosen"], flush=True)

    lead = calibrate_leading(deg, ctl, baseline_end, lag["chosen"]["median_t"])
    if lead["chosen"] is None:
        raise SystemExit("no admissible leading rule at nominal FPR -- "
                         "kill-criterion-2 territory, stop and report")
    ch = lead["chosen"]
    print("  leading  :", ch, flush=True)

    tau_rows = kendall_tau_grid(deg, ctl)
    print("  rejected Kendall-tau rule: min FPR at detect_rate>=0.5 = "
          f"{min((r['fpr'] for r in tau_rows if r['detect_rate'] >= 0.5), default=None)}",
          flush=True)

    rule = {
        "_declared": ("Frozen on calibration seeds >= 900000 (n=40 per arm) at "
                      "c=1.0, mid white noise. Applied UNCHANGED at every coverage "
                      "level, noise colour, severity and consent arm. See spec "
                      "section 9 amendment 1 for the trend-rule -> level-rule change."),
        "nominal_fpr_target": NOMINAL_FPR,
        "n_calibration_seeds": N_CAL,
        "scored_seed_ceiling": CAL_SEED_BASE,
        "leading": {
            "statistics": ["rolling_variance", "rolling_lag1_autocorrelation"],
            "detrend": "in-window linear",
            "window": ch["window"],
            "classes": ch["classes"],
            "class_mixing": "mean of per-class standardised statistic "
                            "(variance-family only; no covariance across node pairs)",
            "standardisation": "z against the same run's baseline period "
                               f"[window-1, {baseline_end})",
            "combine": ch["combine"],
            "k": ch["k"],
            "sustain": ch["sustain"],
            "start": baseline_end,
            "min_frac": 0.5,
            "max_gap": 3,
        },
        "lagging": {
            "metric": "program-close output, trailing mean over "
                      f"{base.output_report_window} cycles",
            "baseline_end": baseline_end,
            "k_sd": lag["chosen"]["k_sd"],
            "sustain": lag["chosen"]["sustain"],
            "max_gap": 3,
        },
        "calibration_result": {"leading": ch, "lagging": lag["chosen"]},
    }
    os.makedirs(os.path.join(HERE, "configs"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "aggregates"), exist_ok=True)
    with open(os.path.join(HERE, "configs", "detection_rule.json"), "w") as f:
        json.dump(rule, f, indent=1)
    with open(os.path.join(HERE, "aggregates", "calibration_rule_grids.json"),
              "w") as f:
        json.dump({"nominal_fpr": NOMINAL_FPR, "n_cal_seeds": N_CAL,
                   "baseline_end": baseline_end,
                   "leading_level_grid": lead["grid"],
                   "lagging_grid": lag["grid"],
                   "rejected_kendall_tau_grid": tau_rows}, f, indent=1)
    print("wrote configs/detection_rule.json")


if __name__ == "__main__":
    main()
