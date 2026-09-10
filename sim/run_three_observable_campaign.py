"""Score three observed-subgraph detectors under a disjoint seed split.

Run block A first. It writes frozen rules using controls only. Run block B only
after block A completes; block B supplies every evaluated rate and lead result.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os

import numpy as np

from indicators import alarm_from_z, baseline_z, lagging_alarm
from model import Config, run_config
from runner import HERE, load_rule, pmap, seed_for
from scoring import bootstrap_mean

OBSERVABLES = ("obs_giant_frac", "obs_mean_degree", "obs_susceptibility")
TRUTH = {"obs_giant_frac": "giant_frac", "obs_mean_degree": "mean_degree",
         "obs_susceptibility": "susceptibility"}
C_GRID = (0.25, 0.5, 0.75, 1.0)
COLOURS = ("white", "red", "systematic")
SEVERITIES = ("low", "mid", "high")
ARMS = ("mnar", "mcar")
REPS = 30
K_GRID = [round(2.0 + 0.25 * i, 3) for i in range(57)]
ALPHA = 0.05
Z95 = 1.959963985
CAL_AGG = "three_observable_block_a_calibration_20260903.json"
EVAL_AGG = "three_observable_block_b_results_20260903.json"
RULE_FILES = {x: f"detection_rule_{x}_split_20260903.json" for x in OBSERVABLES}


def write_new(path, payload):
    with open(path, "x", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=1)
        f.write("\n")


def detail(series, rule):
    lead = rule["leading"]
    z = baseline_z(series, int(lead["window"]), tuple(lead["classes"]),
                   int(lead["start"]), float(lead["min_frac"]),
                   int(lead["max_gap"]))
    score = 0.5 * (z["z_var"] + z["z_ac"])
    alarm = alarm_from_z(z["z_var"], z["z_ac"], float(lead["k"]),
                         int(lead["sustain"]), lead["combine"],
                         int(lead["start"]))
    y = score.copy()
    y[~np.isfinite(y)] = -np.inf
    sustain = int(lead["sustain"])
    start = int(lead["start"])
    sustained = max(y[i:i + sustain].min()
                    for i in range(start, y.size - sustain + 1))
    return {"alarm_cycle": alarm["alarm_cycle"],
            "max_score": alarm["max_score"],
            "sustained_score": float(sustained),
            "score_nan": int(np.count_nonzero(~np.isfinite(score[start:]))),
            "score_n": int(score.size - start)}


def scored(d, lag_cycle, n_cycles):
    t = d["alarm_cycle"]
    et = n_cycles if t is None else t
    el = n_cycles if lag_cycle is None else lag_cycle
    return {"t_lead": t, "t_lag": lag_cycle, "dt_censored": el - et,
            "dt_paired": lag_cycle - t if t is not None and lag_cycle is not None else None,
            "lead_max_score": d["max_score"], "lead_fired": t is not None,
            "lag_fired": lag_cycle is not None, "score_nan": d["score_nan"],
            "score_n": d["score_n"]}


def base_jobs(block, scenarios):
    return [(c, arm, colour, sev, scen, rep)
            for c, arm, colour, sev, scen in itertools.product(
                C_GRID, ARMS, COLOURS, SEVERITIES, scenarios)
            for rep in range(REPS)
            if seed_for(f"ea|{c}|{arm}", rep) % 2 == block]


def config(job):
    c, arm, colour, sev, scen, rep = job
    return Config(consent_c=c, consent_arm=arm, eps_colour=colour,
                  eps_severity=sev, scenario=scen,
                  seed=seed_for(f"ea|{c}|{arm}", rep))


def calibration_job(job):
    cfg = config(job)
    s = run_config(cfg).series()
    base = load_rule("detection_rule_recal.json")
    latency = detail(s, base)
    row = {"c": job[0], "arm": job[1], "colour": job[2], "sev": job[3],
           "rep": job[5], "seed": cfg.seed,
           "latency_fired": latency["alarm_cycle"] is not None}
    for observable in OBSERVABLES:
        rule = {**base, "leading": {**base["leading"], "classes": [observable]}}
        row[observable] = {"sustained_score": detail(s, rule)["sustained_score"]}
    return row


def make_rule(observable, threshold, false, target, n):
    base = load_rule("detection_rule_recal.json")
    return {
        "_declared": (
            "Pre-declared 2026-09-03. Even seeds are block A; odd seeds are block B. "
            "Threshold chosen on block A controls only from the fixed 2.00-16.00 "
            "quarter-step grid to minimise absolute false-alarm-count difference "
            "from latency; a tie selects the larger threshold. Frozen before block B."),
        "source_rule_form": "detection_rule_recal.json",
        "split": {"block_a": "even integer seeds", "block_b": "odd integer seeds"},
        "leading": {**base["leading"], "classes": [observable], "k": threshold},
        "lagging": base["lagging"],
        "calibration_result": {"observable": observable, "k": threshold,
                               "false_alarms": false, "latency_false_alarms": target,
                               "n_block_a_controls": n},
    }


def calibrate(workers):
    jobs = base_jobs(0, ("control",))
    rows = pmap(calibration_job, jobs, workers, "block A controls")
    n = len(rows)
    target = sum(r["latency_fired"] for r in rows)
    rules, curves = {}, {}
    for observable in OBSERVABLES:
        scores = np.asarray([r[observable]["sustained_score"] for r in rows])
        curve = [{"k": k, "false_alarms": int(np.count_nonzero(scores >= k))}
                 for k in K_GRID]
        selected = min(curve, key=lambda x: (abs(x["false_alarms"] - target),
                                             -x["k"]))
        rule = make_rule(observable, selected["k"], selected["false_alarms"],
                         target, n)
        write_new(os.path.join(HERE, "configs", RULE_FILES[observable]), rule)
        rules[observable] = rule
        curves[observable] = curve
    payload = {"_experiment": "block A false-alarm matching for three observables",
               "design": {"partition": "even seeds only", "n_controls": n,
                          "n_unique_seeds": len({r["seed"] for r in rows}),
                          "latency_false_alarms": target, "k_grid": K_GRID,
                          "selection": "minimum absolute count difference; larger k breaks ties"},
               "rules": rules, "curves": curves, "per_run_rows": rows}
    write_new(os.path.join(HERE, "aggregates", CAL_AGG), payload)
    print(f"block A latency: {target}/{n}")
    for observable in OBSERVABLES:
        x = rules[observable]["calibration_result"]
        print(f"{observable}: k={x['k']:.2f}, {x['false_alarms']}/{n}")


_RULES = None


def rules():
    global _RULES
    if _RULES is None:
        _RULES = {x: load_rule(RULE_FILES[x]) for x in OBSERVABLES}
        _RULES["latency"] = load_rule("detection_rule_recal.json")
    return _RULES


def corr(x, y):
    good = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(good) < 4 or np.std(x[good]) < 1e-12 or np.std(y[good]) < 1e-12:
        return float("nan")
    return float(np.corrcoef(x[good], y[good])[0, 1])


def evaluation_job(job):
    cfg = config(job)
    out = run_config(cfg)
    s = out.series()
    lag_cycle = lagging_alarm(s["output_close"], rules()["latency"]["lagging"])["alarm_cycle"]
    ld = detail(s, rules()["latency"])
    row = {"c": job[0], "arm": job[1], "colour": job[2], "sev": job[3],
           "scen": job[4], "rep": job[5], "seed": cfg.seed,
           "n_cycles": cfg.n_cycles, "latency": scored(ld, lag_cycle, cfg.n_cycles)}
    for observable in OBSERVABLES:
        d = detail(s, rules()[observable])
        entry = scored(d, lag_cycle, cfg.n_cycles)
        entry["raw_nan"] = int(np.count_nonzero(~np.isfinite(s[observable])))
        entry["raw_n"] = int(s[observable].size)
        entry["corr_with_full"] = corr(s[observable], s[TRUTH[observable]])
        row[observable] = entry
    return row


def wilson(k, n):
    p = k / n
    centre = (p + Z95 * Z95 / (2 * n)) / (1 + Z95 * Z95 / n)
    half = Z95 * math.sqrt(p * (1 - p) / n + Z95 * Z95 / (4 * n * n))
    half /= 1 + Z95 * Z95 / n
    return [centre - half, centre + half]


def ztest(k1, n1, k2, n2):
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    z = (k1 / n1 - k2 / n2) / se if se else 0.0
    return {"z": z, "p_two_sided": math.erfc(abs(z) / math.sqrt(2.0))}


def summary(values):
    x = bootstrap_mean(values)
    return {"mean": x["mean"], "bootstrap_se": x["se"],
            "percentile_95_ci": [x["lo95"], x["hi95"]], "n": x["n"]}


def evaluate(workers):
    for filename in RULE_FILES.values():
        if not os.path.exists(os.path.join(HERE, "configs", filename)):
            raise SystemExit("block A rules do not exist; run block-a first")
    jobs = base_jobs(1, ("degrade", "control"))
    rows = pmap(evaluation_job, jobs, workers, "block B evaluation")
    controls = [r for r in rows if r["scen"] == "control"]
    detectors = ("latency",) + OBSERVABLES
    pooled = {}
    for detector in detectors:
        k = sum(r[detector]["lead_fired"] for r in controls)
        pooled[detector] = {"false_alarms": k, "n_controls": len(controls),
                            "rate": k / len(controls), "wilson95": wilson(k, len(controls))}
    lk = pooled["latency"]["false_alarms"]
    ln = pooled["latency"]["n_controls"]
    for observable in OBSERVABLES:
        x = pooled[observable]
        test = ztest(x["false_alarms"], x["n_controls"], lk, ln)
        x["versus_latency_ztest"] = test
        x["matches_latency_alpha_0_05"] = test["p_two_sided"] >= ALPHA

    by_coverage_fpr = {}
    matrix = []
    for c in C_GRID:
        ctl = [r for r in controls if r["c"] == c]
        by_coverage_fpr[str(c)] = {}
        for detector in detectors:
            k = sum(r[detector]["lead_fired"] for r in ctl)
            by_coverage_fpr[str(c)][detector] = {
                "false_alarms": k, "n_controls": len(ctl), "rate": k / len(ctl),
                "wilson95": wilson(k, len(ctl))}
        deg = [r for r in rows if r["c"] == c and r["arm"] == "mnar"
               and r["colour"] == "white" and r["sev"] == "mid"
               and r["scen"] == "degrade"]
        cell = {"c": c, "arm": "mnar", "colour": "white", "severity": "mid",
                "n_reps": len(deg), "detectors": {}}
        for detector in detectors:
            matched = detector == "latency" or pooled[detector]["matches_latency_alpha_0_05"]
            cell["detectors"][detector] = {
                "comparable": matched,
                "lead_over_lagging": summary([r[detector]["dt_censored"] for r in deg])
                if matched else None,
                "coverage_false_alarms": by_coverage_fpr[str(c)][detector]}
        matrix.append(cell)

    relative = {}
    classifications = {}
    full = [r for r in rows if r["c"] == 1.0 and r["arm"] == "mnar"
            and r["colour"] == "white" and r["sev"] == "mid"
            and r["scen"] == "degrade"]
    for observable in OBSERVABLES:
        if not pooled[observable]["matches_latency_alpha_0_05"]:
            relative[observable] = None
            classifications[observable] = "not comparable"
            continue
        vals = []
        for r in full:
            tg = r["n_cycles"] if r[observable]["t_lead"] is None else r[observable]["t_lead"]
            tl = r["n_cycles"] if r["latency"]["t_lead"] is None else r["latency"]["t_lead"]
            vals.append(tl - tg)
        relative[observable] = summary(vals)
        lo, hi = relative[observable]["percentile_95_ci"]
        classifications[observable] = "beats" if lo > 0 else ("loses" if hi < 0 else "matches")

    stability = {}
    for c in C_GRID:
        cr = [r for r in rows if r["c"] == c]
        stability[str(c)] = {}
        for observable in OBSERVABLES:
            correlations = [r[observable]["corr_with_full"] for r in cr]
            finite_corr = [x for x in correlations if np.isfinite(x)]
            stability[str(c)][observable] = {
                "undefined_score_cycles": sum(r[observable]["score_nan"] for r in cr),
                "scored_cycles": sum(r[observable]["score_n"] for r in cr),
                "raw_undefined_cycles": sum(r[observable]["raw_nan"] for r in cr),
                "raw_cycles": sum(r[observable]["raw_n"] for r in cr),
                "median_within_run_correlation": float(np.median(finite_corr)),
                "n_runs_with_finite_correlation": len(finite_corr), "n_runs": len(cr)}

    payload = {
        "_experiment": "three observed-subgraph detectors evaluated on block B",
        "design": {"partition": "odd seeds only", "n_runs": len(rows),
                   "n_controls": len(controls), "n_unique_seeds": len({r["seed"] for r in rows}),
                   "false_alarm_match": "two-sided two-proportion z test, alpha 0.05, per rule",
                   "censoring": "non-firing detector assigned cycle 320, identical for all detectors",
                   "uncertainty": "2000-draw deterministic run bootstrap; standard error is bootstrap SD; interval is percentile 95% CI",
                   "matrix_cell": "MNAR, white error, mid severity, degradation, block B"},
        "rules": rules(), "pooled_false_alarms": pooled,
        "false_alarms_by_coverage": by_coverage_fpr, "results_matrix": matrix,
        "full_coverage_relative_to_latency": relative,
        "full_coverage_classification": classifications,
        "stability_by_coverage": stability, "per_run_rows": rows}
    write_new(os.path.join(HERE, "aggregates", EVAL_AGG), payload)
    print("REPRODUCTION 40/40 rows, 320/320 fields")
    for observable in OBSERVABLES:
        x = pooled[observable]
        print(f"{observable}: {x['false_alarms']}/{x['n_controls']} "
              f"({100*x['rate']:.2f}%); {classifications[observable]} at full coverage")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("block-a", "block-b"))
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()
    calibrate(args.workers) if args.phase == "block-a" else evaluate(args.workers)


if __name__ == "__main__":
    main()
