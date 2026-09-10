"""Calibrate and evaluate observed-graph early-warning detectors.

Uses the existing detector form, quarter-step threshold grid, 3.5% calibration
false-alarm target, censoring rule, bootstrap, headline design, and seed discipline.
Every output name is new.
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
from runner import HERE, SCORED_SEED_CEILING, load_rule, pmap, seed_for, write_agg
from scoring import bootstrap_mean, score_run

OBSERVABLES = ("obs_giant_frac", "obs_mean_degree")
C_GRID = (0.25, 0.5, 0.75, 1.0)
COLOURS = ("white", "red", "systematic")
SEVERITIES = ("low", "mid", "high")
ARMS = ("mnar", "mcar")
SCENARIOS = ("degrade", "control")
REPS = 30
NULL_SEED_BASE = 960_000
N_NULL = 600
TARGET_FPR = 0.035
K_GRID = [round(2.0 + 0.25 * i, 3) for i in range(57)]
CAL_AGG = "graph_detector_calibration.json"
SWEEP_AGG = "graph_detector_headline.json"
RULE_FILES = {
    "obs_giant_frac": "detection_rule_obs_giant_frac.json",
    "obs_mean_degree": "detection_rule_obs_mean_degree.json",
}
Z95 = 1.959963985


def max_sustained_min(x, sustain, start):
    y = np.asarray(x, float).copy()
    y[~np.isfinite(y)] = -np.inf
    return float(max(y[i:i + sustain].min()
                     for i in range(start, y.size - sustain + 1)))


def detail(series, rule):
    lead = rule["leading"]
    z = baseline_z(series, int(lead["window"]), tuple(lead["classes"]),
                   int(lead["start"]), float(lead["min_frac"]),
                   int(lead["max_gap"]))
    score = 0.5 * (z["z_var"] + z["z_ac"])
    alarm = alarm_from_z(z["z_var"], z["z_ac"], float(lead["k"]),
                         int(lead["sustain"]), lead["combine"],
                         int(lead["start"]))
    start = int(lead["start"])
    return {
        "alarm_cycle": alarm["alarm_cycle"],
        "max_score": alarm["max_score"],
        "sustained_score": max_sustained_min(score, int(lead["sustain"]), start),
        "score_nan_scored": int(np.count_nonzero(~np.isfinite(score[start:]))),
        "score_scored_cycles": int(score.size - start),
    }


def cal_job(seed):
    cfg = Config(scenario="control", seed=seed, consent_c=1.0,
                 consent_arm="mnar", eps_colour="white", eps_severity="mid")
    series = run_config(cfg).series()
    base = load_rule("detection_rule_recal.json")
    row = {"seed": seed, "n_cycles": cfg.n_cycles}
    for observable in OBSERVABLES:
        rule = {**base, "leading": {**base["leading"], "classes": [observable]}}
        d = detail(series, rule)
        row[observable] = {
            "sustained_score": d["sustained_score"],
            "raw_nan_cycles": int(np.count_nonzero(~np.isfinite(series[observable]))),
            "score_nan_scored": d["score_nan_scored"],
            "score_scored_cycles": d["score_scored_cycles"],
        }
    return row


def make_rule(observable, threshold, n_false):
    base = load_rule("detection_rule_recal.json")
    declared = (
        f"CALIBRATED 2026-09-03 on {N_NULL} dedicated control runs, seed block "
        f"{NULL_SEED_BASE}-{NULL_SEED_BASE + N_NULL - 1}, disjoint from earlier "
        "calibration blocks and every scored seed below 900000. Configuration "
        "c=1.0, white/mid, MNAR, identical to the existing recalibration. "
        f"The observable is {observable}. All other leading-rule fields are copied "
        "unchanged from detection_rule_recal.json. Selection is the same declared "
        "procedure: the smallest k on the 2.00 to 16.00 quarter-step grid with "
        "calibration-null FPR at or below 3.5%. The rule is applied unchanged at "
        "every coverage, error colour, severity, consent arm, and scenario."
    )
    return {
        "_declared": declared,
        "source_rule_form": "detection_rule_recal.json",
        "nominal_fpr_target": TARGET_FPR,
        "n_calibration_seeds": N_NULL,
        "calibration_seed_base": NULL_SEED_BASE,
        "scored_seed_ceiling": SCORED_SEED_CEILING,
        "leading": {**base["leading"], "classes": [observable], "k": threshold},
        "lagging": base["lagging"],
        "calibration_result": {
            "observable": observable, "k": threshold,
            "calibration_fpr": n_false / N_NULL,
            "n_false_alarms": n_false, "n_null": N_NULL,
        },
    }


def calibrate(workers):
    rows = pmap(cal_job, [NULL_SEED_BASE + i for i in range(N_NULL)], workers,
                "observed-graph calibration null")
    curves, chosen, rules = {}, {}, {}
    for observable in OBSERVABLES:
        stats = np.asarray([r[observable]["sustained_score"] for r in rows])
        curve = [{"k": k, "fpr": float(np.mean(stats >= k)),
                  "n_false_alarms": int(np.count_nonzero(stats >= k))}
                 for k in K_GRID]
        admissible = [x for x in curve if x["fpr"] <= TARGET_FPR]
        if not admissible:
            raise SystemExit(f"no threshold on the declared grid calibrates {observable}")
        selected = min(admissible, key=lambda x: x["k"])
        rule = make_rule(observable, selected["k"], selected["n_false_alarms"])
        path = os.path.join(HERE, "configs", RULE_FILES[observable])
        with open(path, "x", encoding="utf-8", newline="\n") as f:
            json.dump(rule, f, indent=1)
            f.write("\n")
        curves[observable], chosen[observable], rules[observable] = curve, selected, rule

    write_agg(CAL_AGG, {
        "_experiment": "observed-graph detector calibration on a dedicated null",
        "design": {
            "n_null_control_runs": N_NULL,
            "null_seed_base": NULL_SEED_BASE,
            "null_seed_last": NULL_SEED_BASE + N_NULL - 1,
            "config": {"c": 1.0, "arm": "mnar", "colour": "white",
                       "severity": "mid", "scenario": "control"},
            "target_fpr": TARGET_FPR, "k_grid": K_GRID,
            "selection_rule": (
                "independently for each named observable, smallest k on the existing "
                "quarter-step grid with null FPR at or below 3.5%; no selection "
                "between observables"),
            "rule_form_unchanged_except_observable_and_threshold": True,
        },
        "rules": rules, "fpr_curves": curves, "chosen": chosen,
        "per_run_rows": rows,
    })
    for observable in OBSERVABLES:
        x = chosen[observable]
        print(f"  {observable}: k={x['k']:.2f}, calibration FPR "
              f"{x['n_false_alarms']}/{N_NULL}={x['fpr']:.4%}")


_RULES = {}


def rules():
    if not _RULES:
        for observable, filename in RULE_FILES.items():
            _RULES[observable] = load_rule(filename)
        _RULES["latency"] = load_rule("detection_rule_recal.json")
    return _RULES


def corr(x, y):
    good = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(good) < 4 or np.std(x[good]) < 1e-12 or np.std(y[good]) < 1e-12:
        return float("nan")
    return float(np.corrcoef(x[good], y[good])[0, 1])


def scored_from_detail(lead, lag_cycle, n_cycles):
    t_lead = lead["alarm_cycle"]
    eff_lead = n_cycles if t_lead is None else t_lead
    eff_lag = n_cycles if lag_cycle is None else lag_cycle
    return {
        "t_lead": t_lead, "t_lag": lag_cycle,
        "dt_censored": eff_lag - eff_lead,
        "dt_paired": (lag_cycle - t_lead)
        if t_lead is not None and lag_cycle is not None else None,
        "lead_max_score": lead["max_score"],
        "lead_fired": t_lead is not None,
        "lag_fired": lag_cycle is not None,
    }


def sweep_job(job):
    c, arm, colour, severity, scenario, rep = job
    cfg = Config(consent_c=c, consent_arm=arm, eps_colour=colour,
                 eps_severity=severity, scenario=scenario,
                 seed=seed_for(f"ea|{c}|{arm}", rep))
    output = run_config(cfg)
    series, detector_rules = output.series(), rules()
    lag_cycle = lagging_alarm(series["output_close"],
                              detector_rules["latency"]["lagging"])["alarm_cycle"]
    ld = detail(series, detector_rules["latency"])
    latency = scored_from_detail(ld, lag_cycle, cfg.n_cycles)
    row = {
        "c": c, "arm": arm, "colour": colour, "sev": severity,
        "scen": scenario, "rep": rep, "seed": cfg.seed,
        "n_cycles": cfg.n_cycles,
        "n_consented_nodes": int(sum(output.consent_status)),
        "latency": latency,
        "latency_raw_nan_cycles": int(np.count_nonzero(
            ~np.isfinite(series["obs_formation"]))),
        "latency_score_nan_scored": ld["score_nan_scored"],
        "latency_score_scored_cycles": ld["score_scored_cycles"],
    }
    for observable in OBSERVABLES:
        d = detail(series, detector_rules[observable])
        scored = scored_from_detail(d, lag_cycle, cfg.n_cycles)
        raw = series[observable]
        truth = "giant_frac" if observable == "obs_giant_frac" else "mean_degree"
        entry = {
            **scored,
            "raw_nan_cycles": int(np.count_nonzero(~np.isfinite(raw))),
            "raw_mean": float(np.nanmean(raw)),
            "raw_sd": float(np.nanstd(raw, ddof=1)),
            "score_nan_scored": d["score_nan_scored"],
            "score_scored_cycles": d["score_scored_cycles"],
            "corr_with_full_graph": corr(raw, series[truth]),
        }
        if observable == "obs_mean_degree":
            edges = raw * row["n_consented_nodes"] / 2.0
            entry["observed_edges_mean"] = float(np.nanmean(edges))
            entry["observed_edges_min"] = float(np.nanmin(edges))
            entry["observed_edges_max"] = float(np.nanmax(edges))
        row[observable] = entry
    return row


def wilson(k, n):
    p = k / n
    centre = (p + Z95 * Z95 / (2 * n)) / (1 + Z95 * Z95 / n)
    half = Z95 * math.sqrt(p * (1 - p) / n + Z95 * Z95 / (4 * n * n))
    half /= 1 + Z95 * Z95 / n
    return [centre - half, centre + half]


def lead_summary(values):
    b = bootstrap_mean(values)
    return {"mean": b["mean"], "bootstrap_se": b["se"],
            "bootstrap_percentile_95": [b["lo95"], b["hi95"]], "n": b["n"]}


def sweep(workers):
    for filename in RULE_FILES.values():
        if not os.path.exists(os.path.join(HERE, "configs", filename)):
            raise SystemExit("run the calibration phase before the sweep")
    jobs = [(c, arm, colour, severity, scenario, rep)
            for c, arm, colour, severity, scenario in itertools.product(
                C_GRID, ARMS, COLOURS, SEVERITIES, SCENARIOS)
            for rep in range(REPS)]
    rows = pmap(sweep_job, jobs, workers, "observed-graph headline sweep")
    cells = []
    for c, arm, colour, severity in itertools.product(C_GRID, ARMS, COLOURS, SEVERITIES):
        degrade = [r for r in rows if r["c"] == c and r["arm"] == arm
                   and r["colour"] == colour and r["sev"] == severity
                   and r["scen"] == "degrade"]
        controls = [r for r in rows if r["c"] == c and r["arm"] == arm
                    and r["colour"] == colour and r["sev"] == severity
                    and r["scen"] == "control"]
        cell = {"c": c, "arm": arm, "colour": colour,
                "severity": severity, "n_reps": len(degrade), "detectors": {}}
        latency_false = sum(r["latency"]["lead_fired"] for r in controls)
        cell["detectors"]["latency"] = {
            "lead_over_lagging": lead_summary(
                [r["latency"]["dt_censored"] for r in degrade]),
            "detect_rate": float(np.mean([r["latency"]["lead_fired"] for r in degrade])),
            "false_alarms": int(latency_false), "n_controls": len(controls),
            "fpr": latency_false / len(controls),
        }
        for observable in OBSERVABLES:
            relative = []
            for r in degrade:
                tg = r["n_cycles"] if r[observable]["t_lead"] is None else r[observable]["t_lead"]
                tl = r["n_cycles"] if r["latency"]["t_lead"] is None else r["latency"]["t_lead"]
                relative.append(tl - tg)
            false = sum(r[observable]["lead_fired"] for r in controls)
            cell["detectors"][observable] = {
                "lead_over_lagging": lead_summary(
                    [r[observable]["dt_censored"] for r in degrade]),
                "lead_relative_to_latency": lead_summary(relative),
                "relative_sign": "positive means graph alarm is earlier",
                "detect_rate": float(np.mean([r[observable]["lead_fired"] for r in degrade])),
                "false_alarms": int(false), "n_controls": len(controls),
                "fpr": false / len(controls),
            }
        cells.append(cell)

    controls = [r for r in rows if r["scen"] == "control"]
    pooled = {}
    for detector in ("latency",) + OBSERVABLES:
        false = sum(r[detector]["lead_fired"] for r in controls)
        pooled[detector] = {
            "n_false_alarms": int(false), "n_control_runs": len(controls),
            "fpr": false / len(controls), "wilson95": wilson(false, len(controls)),
        }

    missingness, stability = {}, {}
    for c in C_GRID:
        cr = [r for r in rows if r["c"] == c]
        block = {
            "n_runs": len(cr), "n_raw_cycles_per_detector": len(cr) * 320,
            "n_scored_cycles_per_detector": len(cr) * (320 - 90),
            "latency_raw_nan_cycles": int(sum(r["latency_raw_nan_cycles"] for r in cr)),
            "latency_score_nan_scored": int(sum(r["latency_score_nan_scored"] for r in cr)),
        }
        for observable in OBSERVABLES:
            block[f"{observable}_raw_nan_cycles"] = int(sum(
                r[observable]["raw_nan_cycles"] for r in cr))
            block[f"{observable}_score_nan_scored"] = int(sum(
                r[observable]["score_nan_scored"] for r in cr))
        missingness[str(c)] = block

        sr = [r for r in cr if r["colour"] == "white" and r["sev"] == "mid"]
        stability[str(c)] = {
            "n_runs": len(sr), "n_cycles": len(sr) * 320,
            "n_unique_seeds": len(set(r["seed"] for r in sr)),
            "obs_giant_run_mean_sd": float(np.std(
                [r["obs_giant_frac"]["raw_mean"] for r in sr], ddof=1)),
            "obs_giant_run_mean_range": [
                float(min(r["obs_giant_frac"]["raw_mean"] for r in sr)),
                float(max(r["obs_giant_frac"]["raw_mean"] for r in sr))],
            "obs_giant_median_within_run_sd": float(np.median(
                [r["obs_giant_frac"]["raw_sd"] for r in sr])),
            "obs_giant_corr_full_graph_median": float(np.nanmedian(
                [r["obs_giant_frac"]["corr_with_full_graph"] for r in sr])),
            "observed_edges_mean": float(np.mean(
                [r["obs_mean_degree"]["observed_edges_mean"] for r in sr])),
            "observed_edges_min": float(min(
                r["obs_mean_degree"]["observed_edges_min"] for r in sr)),
            "observed_edges_max": float(max(
                r["obs_mean_degree"]["observed_edges_max"] for r in sr)),
        }

    headline = [x for x in cells if x["arm"] == "mnar"
                and x["colour"] == "white" and x["severity"] == "mid"]
    write_agg(SWEEP_AGG, {
        "_experiment": "observed-graph detectors on the full headline design",
        "design": {
            "c_grid": list(C_GRID), "colours": list(COLOURS),
            "severities": list(SEVERITIES), "arms": list(ARMS),
            "scenarios": list(SCENARIOS), "reps": REPS, "n_runs": len(rows),
            "paired_control": "matched seeds; seed depends on (c, arm, replicate) only",
            "censoring": "each non-firing detector is assigned n_cycles=320",
            "uncertainty": (
                "2,000 deterministic run-resampling bootstrap draws; standard error "
                "is the bootstrap standard deviation and interval is the percentile 95% interval"),
            "latency_rule_file": "detection_rule_recal.json",
            "graph_rule_files": RULE_FILES,
        },
        "rules": rules(), "pooled_false_alarms": pooled,
        "headline_mnar_white_mid": sorted(headline, key=lambda x: x["c"]),
        "cells": cells, "missingness_by_coverage": missingness,
        "stability_by_coverage_white_mid": stability, "per_run_rows": rows,
    })

    print("\n  pooled evaluated false alarms")
    for detector, value in pooled.items():
        lo, hi = value["wilson95"]
        print(f"  {detector:20s} {value['n_false_alarms']}/{value['n_control_runs']} "
              f"= {value['fpr']:.4%}, Wilson 95% [{lo:.4%}, {hi:.4%}]")
    print("\n  primary white/mid MNAR: graph lead relative to latency")
    for cell in sorted(headline, key=lambda x: x["c"]):
        x = cell["detectors"]["obs_giant_frac"]["lead_relative_to_latency"]
        print(f"  c={cell['c']:.2f}: {x['mean']:+.2f} cycles, bootstrap s.e. "
              f"{x['bootstrap_se']:.2f}, percentile 95% "
              f"[{x['bootstrap_percentile_95'][0]:+.2f}, "
              f"{x['bootstrap_percentile_95'][1]:+.2f}], n={x['n']}")


def low_coverage_job(job):
    arm, scenario, rep = job
    cfg = Config(consent_c=0.25, consent_arm=arm, eps_colour="white",
                 eps_severity="mid", scenario=scenario,
                 seed=seed_for(f"ea|0.25|{arm}", rep))
    output = run_config(cfg)
    series = output.series()
    n_nodes = int(sum(output.consent_status))
    edges = series["obs_mean_degree"] * n_nodes / 2.0
    return {
        "arm": arm, "scenario": scenario, "rep": rep, "seed": cfg.seed,
        "n_cycles": cfg.n_cycles, "n_consented_nodes": n_nodes,
        "zero_edge_cycles": int(np.count_nonzero(edges == 0)),
        "finite_edge_cycles": int(np.count_nonzero(np.isfinite(edges))),
        "edge_mean": float(np.nanmean(edges)),
        "edge_min": float(np.nanmin(edges)), "edge_max": float(np.nanmax(edges)),
        "obs_giant_nan_cycles": int(np.count_nonzero(
            ~np.isfinite(series["obs_giant_frac"]))),
        "obs_giant_mean": float(np.nanmean(series["obs_giant_frac"])),
        "obs_giant_sd": float(np.nanstd(series["obs_giant_frac"], ddof=1)),
        "obs_giant_corr_full": corr(series["obs_giant_frac"], series["giant_frac"]),
    }


def low_coverage_diagnostic(workers):
    jobs = [(arm, scenario, rep)
            for arm, scenario in itertools.product(ARMS, SCENARIOS)
            for rep in range(REPS)]
    rows = pmap(low_coverage_job, jobs, workers, "low-coverage stability diagnostic")
    total_cycles = sum(r["n_cycles"] for r in rows)
    total_zero = sum(r["zero_edge_cycles"] for r in rows)
    payload = {
        "_experiment": "low-coverage observed-graph stability diagnostic",
        "design": {
            "coverage": 0.25, "colour": "white", "severity": "mid",
            "arms": list(ARMS), "scenarios": list(SCENARIOS), "reps": REPS,
            "n_runs": len(rows), "n_cycles": total_cycles,
            "seed_rule": "same matched headline seeds",
        },
        "summary": {
            "zero_edge_cycles": total_zero, "n_cycles": total_cycles,
            "zero_edge_fraction": total_zero / total_cycles,
            "runs_with_zero_edge_cycles": sum(r["zero_edge_cycles"] > 0 for r in rows),
            "n_runs": len(rows),
            "observed_edges_mean": float(np.mean([r["edge_mean"] for r in rows])),
            "observed_edges_min": float(min(r["edge_min"] for r in rows)),
            "observed_edges_max": float(max(r["edge_max"] for r in rows)),
            "consented_nodes_min": min(r["n_consented_nodes"] for r in rows),
            "consented_nodes_max": max(r["n_consented_nodes"] for r in rows),
            "obs_giant_nan_cycles": sum(r["obs_giant_nan_cycles"] for r in rows),
            "obs_giant_run_mean_sd": float(np.std(
                [r["obs_giant_mean"] for r in rows], ddof=1)),
            "obs_giant_run_mean_range": [
                float(min(r["obs_giant_mean"] for r in rows)),
                float(max(r["obs_giant_mean"] for r in rows))],
            "obs_giant_median_within_run_sd": float(np.median(
                [r["obs_giant_sd"] for r in rows])),
            "obs_giant_corr_full_median": float(np.nanmedian(
                [r["obs_giant_corr_full"] for r in rows])),
        },
        "per_run_rows": rows,
    }
    write_agg("graph_detector_low_coverage_stability.json", payload)
    print(json.dumps(payload["summary"], indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("calibrate", "sweep", "low-coverage"))
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()
    calibrate(args.workers) if args.phase == "calibrate" else (sweep(args.workers) if args.phase == "sweep" else low_coverage_diagnostic(args.workers))


if __name__ == "__main__":
    main()
