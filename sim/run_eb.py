"""
coverage-and-coupling campaign -- the joint (c x g) sweep (frozen spec section 4, coupling-response and severed-control design).

The one genuinely unprecedented measurement in the design: does restricted observation
PROTECT a leading indicator from gaming by hiding the target, or COMPOUND its loss?
The univariate g sweep is pre-empted by Bauch et al. (2016), who already report that
behavioural coupling mutes an EWS against an uncoupled control; the joint sweep is
what is new, and it can produce a non-monotone optimum in c that no cited work
predicts.

Grid: c in {0.25, 0.5, 1.0} x g in {0, mid, high}, degradation scenario.
Output: corr(observed leading indicator, latent A(t)).
  SIGN CONVENTION -- the indicator rises as A(t) falls, so a well-tracking instrument
  gives a NEGATIVE correlation; corr -> 0 means the indicator has stopped tracking the
  latent state.  |corr| is the tracking quality.

Coupling-severed control (Bauch et al.'s coupled/uncoupled design): the
identical substrate and seed, with the response DRAWN FROM THE SAME RNG STREAM and then
discarded, so severing the coupling silences that dial without re-phasing the substrate
streams.

Agents respond to the RECORDED PROXY under partial coverage, never to the true state
(spec section 2, dial 4): they preferentially recruit partners whose ties are already
in the record, because those are the ties the published indicator is computed from.
The response includes a withdrawal/attrition channel mirroring Bol et al. (2018)'s
measured participation effect: agents in the bottom quantile of the recorded score
withdraw at a rate proportional to g.

Run:  python run_eb.py [--workers N] [--reps 30]
      python run_eb.py --rule detection_rule_recal.json --out eb_joint_sweep_n500.json \
                       --reps 500
Writes: aggregates/eb_joint_sweep.json (or --out)

FOLLOW-UP CAMPAIGN 2026-08-05 -- resolving the C-III evidence gap (RESULTS.md section 7
item 1).  Nothing about the design changes: the grid, the substrate (N = 400), the
statistic, the coupling implementation and the severed control are all as frozen.  Only
the replicate count rises, from 120 to >= 500 matched pairs per cell, and the scoring
uses the recalibrated threshold.  A DECISION RULE was pre-declared before these runs:

  * both paired coupling-loss CIs (g = mid, g = high) exclude zero and agree in sign
      -> the interaction is RESOLVED; report the direction;
  * both exclude zero but disagree in sign across g
      -> that is itself the result: a coupling-strength-dependent regime change,
         reported with its c-dependence;
  * either CI still straddles zero
      -> UNRESOLVED at this substrate, reported with the power analysis and the minimum
         detectable effect at the achieved number of pairs.

The rule is applied mechanically in `decision` inside the aggregate, so the verdict is
not a post-hoc reading of the table.
"""

from __future__ import annotations

import argparse
import itertools
import math
import os

import numpy as np

from model import Config, run_config
from scoring import indicator_vs_latent, score_run, bootstrap_mean
from runner import pmap, seed_for, write_agg, load_rule

C_GRID = (0.25, 0.5, 1.0)
G_LEVELS = {"0": 0.0, "mid": 0.5, "high": 1.0}
COUPLING = ("coupled", "severed")

_RULE = None


def _rule():
    global _RULE
    if _RULE is None:
        _RULE = load_rule()
    return _RULE


def _job(job):
    c, gname, mode, rep = job
    g = G_LEVELS[gname]
    cfg = Config(consent_c=c, consent_arm="mnar", eps_colour="white",
                 eps_severity="mid", scenario="degrade",
                 coupling_g=g, coupling_severed=(mode == "severed"),
                 seed=seed_for(f"eb|{c}", rep))
    out = run_config(cfg)
    s = out.series()
    iv = indicator_vs_latent(s, _rule())
    sc = score_run(s, _rule(), cfg.n_cycles)
    n_act = np.asarray(out.n_active, float)
    return {"c": c, "g_name": gname, "g": g, "mode": mode, "rep": rep,
            "seed": cfg.seed,
            "corr_indicator_latent": iv["corr_indicator_latent"],
            "corr_var_leg": iv["corr_var_leg"],
            "t_lead": sc["t_lead"], "lead_fired": sc["lead_fired"],
            "dt_censored": sc["dt_censored"],
            "achieved_edge_frac": float(s["achieved_edge_frac"][-1]),
            "active_end_frac": float(n_act[-1] / cfg.n_agents),
            "attrition_frac": float(1.0 - n_act[-1] / cfg.n_agents),
            "withdrawn_total": int(out.withdrawn_total),
            "n_cycles": cfg.n_cycles}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--rule", default=None,
                    help="detection-rule file in configs/ (default detection_rule.json)")
    ap.add_argument("--out", default="eb_joint_sweep.json")
    args = ap.parse_args()
    if args.rule:
        os.environ["SIM_DETECTION_RULE"] = args.rule
    reps = args.reps

    jobs = [(c, g, m, r) for c, g, m, r in itertools.product(
        C_GRID, G_LEVELS, COUPLING, range(reps))]
    rows = pmap(_job, jobs, args.workers, "coverage-and-coupling campaign (c x g) joint sweep")

    cells = []
    for c, g, m in itertools.product(C_GRID, G_LEVELS, COUPLING):
        rs = [r for r in rows if r["c"] == c and r["g_name"] == g
              and r["mode"] == m]
        cr = bootstrap_mean([r["corr_indicator_latent"] for r in rs])
        cells.append({
            "c": c, "g_name": g, "g": G_LEVELS[g], "mode": m, "n_reps": len(rs),
            "corr_mean": cr["mean"], "corr_se": cr["se"],
            "corr_lo95": cr["lo95"], "corr_hi95": cr["hi95"],
            "abs_corr_mean": abs(cr["mean"]),
            "corr_var_leg_mean": float(np.nanmean(
                [r["corr_var_leg"] for r in rs])),
            "lead_detect_rate": float(np.mean([r["lead_fired"] for r in rs])),
            "dt_censored_mean": float(np.mean([r["dt_censored"] for r in rs])),
            "achieved_edge_frac_mean": float(np.mean(
                [r["achieved_edge_frac"] for r in rs])),
            "attrition_frac_mean": float(np.mean([r["attrition_frac"] for r in rs])),
            "withdrawn_total_mean": float(np.mean([r["withdrawn_total"] for r in rs])),
        })
    look = {(x["c"], x["g_name"], x["mode"]): x for x in cells}

    # ---- the interaction: does low coverage protect or compound? --------------
    # tracking loss caused by coupling, at matched seeds and matched coverage:
    #   loss = |corr| severed  -  |corr| coupled     (loss > 0 => coupling degraded it)
    interaction = {}
    for gname in G_LEVELS:
        if G_LEVELS[gname] == 0.0:
            continue
        per_c = {}
        for c in C_GRID:
            cp = {r["rep"]: r for r in rows if r["c"] == c
                  and r["g_name"] == gname and r["mode"] == "coupled"}
            sv = {r["rep"]: r for r in rows if r["c"] == c
                  and r["g_name"] == gname and r["mode"] == "severed"}
            shared = sorted(set(cp) & set(sv))
            paired = [abs(sv[i]["corr_indicator_latent"])
                      - abs(cp[i]["corr_indicator_latent"]) for i in shared]
            b = bootstrap_mean(paired)
            per_c[c] = {"tracking_loss_mean": b["mean"], "se": b["se"],
                        "lo95": b["lo95"], "hi95": b["hi95"], "n_pairs": b["n"],
                        "abs_corr_coupled": look[(c, gname, "coupled")]["abs_corr_mean"],
                        "abs_corr_severed": look[(c, gname, "severed")]["abs_corr_mean"],
                        "significant": bool(b["lo95"] > 0 or b["hi95"] < 0)}
        losses = [per_c[c]["tracking_loss_mean"] for c in C_GRID]   # c ascending
        lo_c, hi_c = losses[0], losses[-1]
        if hi_c - lo_c > 0.02:
            direction = ("PROTECT -- coupling costs less tracking at low coverage "
                         "than at full coverage")
        elif lo_c - hi_c > 0.02:
            direction = ("COMPOUND -- coupling costs more tracking at low coverage "
                         "than at full coverage")
        else:
            direction = "NEITHER -- the coupling cost does not vary with coverage"
        # non-monotone optimum in c for tracking quality under coupling
        aq = [look[(c, gname, "coupled")]["abs_corr_mean"] for c in C_GRID]
        best = int(np.argmax(aq))
        interaction[gname] = {
            "per_c": per_c,
            "loss_at_lowest_c": lo_c, "loss_at_full_c": hi_c,
            "direction": direction,
            "abs_corr_coupled_by_c": dict(zip([str(c) for c in C_GRID], aq)),
            "argmax_abs_corr_c": C_GRID[best],
            "non_monotone_optimum_in_c": bool(0 < best < len(C_GRID) - 1),
        }

    # ---- the pre-declared decision rule (follow-up campaign 2026-08-05) -------
    # "the paired coupling-loss CI at g" is the bootstrap 95% CI of the paired
    # (severed - coupled) |corr| difference POOLED OVER THE THREE COVERAGE LEVELS at
    # that g -- the main effect of coupling at matched seeds.  The c-dependence of the
    # same quantity is reported alongside it and is what carries the interaction.
    # Power: minimum detectable effect at 80% power, two-sided alpha = 0.05, from the
    # observed paired s.d.:  MDE = (z_0.975 + z_0.8) * sd / sqrt(n) = 2.8016 sd/sqrt(n).
    Z80 = 1.959963985 + 0.841621234
    decision = {"rule": ("pre-declared: both CIs exclude zero and agree in sign -> "
                         "RESOLVED (report direction); both exclude zero and disagree "
                         "in sign -> COUPLING-STRENGTH-DEPENDENT REGIME CHANGE (report "
                         "with c-dependence); either straddles zero -> UNRESOLVED "
                         "(report power analysis + MDE)"),
                "ci_definition": ("bootstrap 95% CI of the paired (severed - coupled) "
                                  "|corr| difference, pooled over c at fixed g"),
                "per_g": {}}
    for gname in G_LEVELS:
        if G_LEVELS[gname] == 0.0:
            continue
        pooled = []
        for c in C_GRID:
            cp = {r["rep"]: r for r in rows if r["c"] == c
                  and r["g_name"] == gname and r["mode"] == "coupled"}
            sv = {r["rep"]: r for r in rows if r["c"] == c
                  and r["g_name"] == gname and r["mode"] == "severed"}
            for i in sorted(set(cp) & set(sv)):
                pooled.append(abs(sv[i]["corr_indicator_latent"])
                              - abs(cp[i]["corr_indicator_latent"]))
        b = bootstrap_mean(pooled)
        arr = np.asarray([v for v in pooled if np.isfinite(v)], float)
        sd = float(arr.std(ddof=1)) if arr.size > 1 else float("nan")
        n = int(arr.size)
        decision["per_g"][gname] = {
            "n_pairs_pooled": n,
            "tracking_loss_mean": b["mean"], "se": b["se"],
            "lo95": b["lo95"], "hi95": b["hi95"],
            "paired_sd": sd,
            "excludes_zero": bool(b["lo95"] > 0 or b["hi95"] < 0),
            "sign": int(np.sign(b["mean"])),
            "mde_80pct_power": (Z80 * sd / math.sqrt(n) if n > 1 else float("nan")),
            "n_pairs_needed_for_observed_effect": (
                int(math.ceil((Z80 * sd / abs(b["mean"])) ** 2))
                if abs(b["mean"]) > 1e-12 else None),
            "per_c_tracking_loss": {str(c): interaction[gname]["per_c"][c]
                                    for c in C_GRID},
        }
    pg = decision["per_g"]
    both_excl = all(v["excludes_zero"] for v in pg.values())
    signs = {v["sign"] for v in pg.values()}
    if both_excl and len(signs) == 1:
        s = signs.pop()
        decision["verdict"] = "RESOLVED"
        decision["direction"] = (
            "coupling COSTS tracking quality at every coupling strength (severed "
            "tracks the latent state better than coupled)" if s > 0 else
            "coupling IMPROVES tracking quality at every coupling strength (coupled "
            "tracks the latent state better than severed)")
    elif both_excl:
        decision["verdict"] = "COUPLING-STRENGTH-DEPENDENT REGIME CHANGE"
        decision["direction"] = (
            "the sign of the coupling effect changes with coupling strength: "
            + ", ".join(f"g={g}: {v['tracking_loss_mean']:+.4f} "
                        f"[{v['lo95']:+.4f},{v['hi95']:+.4f}]"
                        for g, v in pg.items()))
    else:
        decision["verdict"] = "UNRESOLVED"
        decision["direction"] = (
            "at least one paired coupling-loss CI still straddles zero at "
            + ", ".join(f"n={v['n_pairs_pooled']} pairs (g={g})"
                        for g, v in pg.items())
            + "; the minimum detectable effect at 80% power is "
            + ", ".join(f"{v['mde_80pct_power']:.4f} (g={g})"
                        for g, v in pg.items()))
    decision["interaction_c_dependence"] = {
        g: {"loss_at_lowest_c": interaction[g]["loss_at_lowest_c"],
            "loss_at_full_c": interaction[g]["loss_at_full_c"],
            "direction": interaction[g]["direction"],
            "per_c_significant": {str(c): interaction[g]["per_c"][c]["significant"]
                                  for c in C_GRID}}
        for g in pg}

    write_agg(args.out, {
        "_experiment": "coverage-and-coupling campaign joint (c x g) sweep -- protection vs compounding",
        "design": {"c_grid": list(C_GRID), "g_levels": G_LEVELS,
                   "modes": list(COUPLING), "reps": reps,
                   "noise": "white / mid", "consent_arm": "mnar",
                   "n_cycles": rows[0]["n_cycles"],
                   "n_agents": 400,
                   "rule_file": os.environ.get("SIM_DETECTION_RULE",
                                               "detection_rule.json"),
                   "n_runs": len(rows),
                   "sign_convention": ("the indicator rises as A(t) falls, so a "
                                       "well-tracking instrument gives a NEGATIVE "
                                       "corr; |corr| is tracking quality"),
                   "severed_control": ("identical substrate and seed; the coupling "
                                       "response is drawn from the same RNG stream "
                                       "and discarded (Bauch et al. 2016 design)"),
                   "note_g0": ("at g=0 the coupled and severed conditions are the "
                               "same run by construction -- the coupling code is "
                               "never entered; both are reported for symmetry")},
        "detection_rule": _rule(),
        "cells": cells,
        "interaction": interaction,
        "decision": decision,
        "per_run_rows": rows,
    })

    print("\n   c     g     mode      corr      se   |corr|  attrition  detectDR")
    for x in sorted(cells, key=lambda y: (y["c"], y["g"], y["mode"])):
        print(f"  {x['c']:.2f} {x['g_name']:5s} {x['mode']:8s} "
              f"{x['corr_mean']:+7.3f} {x['corr_se']:6.3f} {x['abs_corr_mean']:6.3f} "
              f"{x['withdrawn_total_mean']:9.0f} {x['lead_detect_rate']:8.2f}")
    print("\n  interaction (tracking loss from coupling, by coverage):")
    for g, v in interaction.items():
        print(f"    g={g}: loss(c=0.25)={v['loss_at_lowest_c']:+.3f} "
              f"loss(c=1.0)={v['loss_at_full_c']:+.3f} -> {v['direction']}")
        print(f"      |corr| by c: {v['abs_corr_coupled_by_c']} "
              f"argmax c={v['argmax_abs_corr_c']} "
              f"non_monotone={v['non_monotone_optimum_in_c']}")
    print(f"\n  pre-declared decision rule -> {decision['verdict']}")
    for g, v in decision["per_g"].items():
        print(f"    g={g:5s} n_pairs={v['n_pairs_pooled']:4d} "
              f"loss={v['tracking_loss_mean']:+.4f} "
              f"[{v['lo95']:+.4f},{v['hi95']:+.4f}] "
              f"excludes_zero={v['excludes_zero']} "
              f"sd={v['paired_sd']:.4f} MDE80={v['mde_80pct_power']:.4f}")
    print(f"    {decision['direction']}")


if __name__ == "__main__":
    main()
