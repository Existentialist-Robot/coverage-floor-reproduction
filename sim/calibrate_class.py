"""
Threshold calibration for the BOTH-CLASS leading rule (remediation campaign 2026-08-05).

Why this exists: the manuscript
says the leading indicators are "per-class rising variance and lag-1 autocorrelation on
formation AND escalation latency", but the frozen rule's `classes` list is
`["obs_formation"]` only.  Escalation latency is simulated and never scored, so the rule's
advertised `class_mixing` is vacuous -- there is one class.  Two ways to fix a
misstatement: change the prose, or run the rule the prose describes.  This script does
the second, so the choice between them is made on evidence.

What is calibrated, and what is NOT:
  * The only change from configs/detection_rule_recal.json is
    `leading.classes = ["obs_formation", "obs_escalation"]`.  Statistics, window,
    detrending, standardisation, combination, sustain, start, min_frac and max_gap are
    byte-identical, and the LAGGING comparator is inherited verbatim -- the comparator
    must not move when the leading channel changes, or the lead-time contrast is not a
    contrast.
  * The leading threshold is re-picked because adding a signal class changes the null
    distribution of the composite score, so re-using k = 6.5 would compare two rules at
    different false-alarm rates.  Selection is the frozen A-2/A-8 procedure with the
    same target: the SMALLEST k on the same 0.25-step grid whose calibration null FPR is
    <= 3.5%.
  * The null is the SAME 600 control runs on seed block 950000+ that A-8 used (same
    seeds, same configuration c = 1.0 / white / mid / MNAR), so the two rules are
    calibrated on identical data and the only difference between them is the class set.

Run:  python calibrate_class.py [--workers N] [--n-null 600]
Writes: configs/detection_rule_bothclass.json,
        aggregates/calibration_null_bothclass.json
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

from calibrate_recal import (K_GRID, NULL_SEED_BASE, TARGET_FPR, EXTERNAL_BAR,
                            max_sustained_min, _cfg, fpr_curve)
from indicators import baseline_z
from model import run_config
from runner import pmap, write_agg, load_rule, HERE, SCORED_SEED_CEILING

BOTH_CLASSES = ("obs_formation", "obs_escalation")


def _job(job):
    """One calibration run reduced to the per-run sufficient statistic of the
    BOTH-CLASS composite score (the same statistic A-8 uses: for a sustain-s level
    rule the alarm exists iff max_i min(score[i..i+s-1]) >= k)."""
    scen, seed = job
    cfg = _cfg(scen, seed)
    s = run_config(cfg).series()
    L = load_rule("detection_rule_recal.json")["leading"]
    z = baseline_z(s, int(L["window"]), BOTH_CLASSES, int(L["start"]),
                   float(L["min_frac"]), int(L["max_gap"]))
    score = 0.5 * (z["z_var"] + z["z_ac"])          # combine = "mean", as frozen
    return {"scen": scen, "seed": int(seed),
            "lead_sustained_score": max_sustained_min(score, int(L["sustain"]),
                                                      int(L["start"])),
            "lead_max_score": float(np.nanmax(score[int(L["start"]):])),
            "n_cycles": cfg.n_cycles}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--n-null", type=int, default=600)
    args = ap.parse_args()
    if args.n_null < 600:
        raise SystemExit("the null must match A-8's 600 runs to be the same null")

    recal = load_rule("detection_rule_recal.json")
    jobs = [("control", NULL_SEED_BASE + r) for r in range(args.n_null)]
    rows = pmap(_job, jobs, args.workers, "both-class calibration null")

    curve = fpr_curve([r["lead_sustained_score"] for r in rows], K_GRID)
    adm = [x for x in curve if x["fpr"] <= TARGET_FPR]
    if not adm:
        raise SystemExit("no k on the grid reaches the 3.5% target for the both-class "
                         "rule -- stop and report")
    k_new = min(x["k"] for x in adm)
    k_row = next(x for x in curve if x["k"] == k_new)
    at_recal_k = next((x for x in curve
                       if abs(x["k"] - float(recal["leading"]["k"])) < 1e-9), None)

    rule = {
        "_declared": (
            "BOTH-CLASS variant of configs/detection_rule_recal.json, calibrated "
            f"2026-08-05 on the SAME {args.n_null}-run null (seed block "
            f"{NULL_SEED_BASE}+, configuration c=1.0 / white / mid / MNAR) that "
            "amendment A-8 used. The ONLY difference from the recalibrated rule is "
            "leading.classes = [obs_formation, obs_escalation]; the lagging comparator "
            "is inherited verbatim. The leading threshold is re-picked by the frozen "
            "A-2/A-8 procedure -- smallest k on the 0.25-step grid with calibration "
            "null FPR <= 3.5% -- because adding a class changes the null distribution "
            "of the composite score. See spec section 9 amendment A-11."),
        "supersedes": None,
        "variant_of": "configs/detection_rule_recal.json",
        "nominal_fpr_target": TARGET_FPR,
        "external_fpr_bar": EXTERNAL_BAR,
        "n_calibration_seeds": args.n_null,
        "calibration_seed_base": NULL_SEED_BASE,
        "scored_seed_ceiling": SCORED_SEED_CEILING,
        "leading": {**recal["leading"], "classes": list(BOTH_CLASSES), "k": k_new},
        "lagging": dict(recal["lagging"]),
        "calibration_result": {
            "leading": {"k": k_new, "calibration_fpr": k_row["fpr"],
                        "n_false_alarms": k_row["n_false_alarms"],
                        "n_null": args.n_null,
                        "formation_only_k": float(recal["leading"]["k"]),
                        "formation_only_k_fpr_on_bothclass_score": (
                            at_recal_k["fpr"] if at_recal_k else None)},
            "lagging": "inherited unchanged from detection_rule_recal.json"},
    }
    with open(os.path.join(HERE, "configs", "detection_rule_bothclass.json"), "w") as f:
        json.dump(rule, f, indent=1)

    write_agg("calibration_null_bothclass.json", {
        "_experiment": ("threshold calibration for the both-class leading rule "
                        "(remediation campaign 2026-08-05, review finding 11)"),
        "design": {
            "n_null_control_runs": args.n_null,
            "null_seed_base": NULL_SEED_BASE,
            "config": {"c": 1.0, "arm": "mnar", "colour": "white",
                       "severity": "mid", "scenario": "control"},
            "target_fpr": TARGET_FPR, "external_bar": EXTERNAL_BAR,
            "k_grid": K_GRID,
            "classes": list(BOTH_CLASSES),
            "same_null_as": "calibration_null_large.json",
            "what_was_repicked": ["leading.k"],
            "what_changed_in_the_rule": ["leading.classes"],
            "selection_rule": ("smallest k with calibration null FPR <= 3.5%, the "
                               "frozen A-2/A-8 procedure at the same target"),
        },
        "formation_only_rule": recal,
        "bothclass_rule": rule,
        "leading_fpr_curve": curve,
        "chosen": {"leading_k": k_new, "leading_calibration_fpr": k_row["fpr"]},
        "formation_only_threshold_on_bothclass_score": {
            "k": float(recal["leading"]["k"]),
            "fpr": at_recal_k["fpr"] if at_recal_k else None,
            "note": ("what the formation-only threshold would have done on the "
                     "both-class score: the reason the threshold must be re-picked "
                     "rather than re-used"),
        },
        "per_run_rows": rows,
    })

    print(f"\n  both-class leading k = {k_new} (null FPR {k_row['fpr']:.4f} on "
          f"{args.n_null} runs); formation-only k = "
          f"{recal['leading']['k']} would run at "
          f"{at_recal_k['fpr'] if at_recal_k else float('nan'):.4f} on this score")
    print("  wrote configs/detection_rule_bothclass.json")


if __name__ == "__main__":
    main()
