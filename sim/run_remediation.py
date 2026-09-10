"""
Remediation campaign 2026-08-05 -- the new runs the review requires.

Everything here is a SENSITIVITY sweep around the frozen headline cell (white proxy
noise, mid severity, MNAR consent, N = 400).  No claim, statistic, detection-rule form,
grid semantics or kill criterion changes; the substrate parameters that are varied are
varied *in order to report the headline's dependence on them*, which is the opposite of
re-tuning.  Spec section 9 amendments A-11 ... A-14.

Seeds are the E-A seeds -- `seed_for("ea|{c}|{arm}", rep)` -- so every variant run is
matched, seed for seed, to the headline run it is being compared against, and the base
variant (`knee0.30_sig0.100`) reproduces `ea_headline_recal.json` bit-for-bit wherever
the two overlap.  The verifier asserts that.

FOUR LEGS
  gain     (finding 2, the crux) The lagging comparator's sensitivity, swept two ways:
           `output_knee` (the cohesion -> output transfer's knee, which sets the slope
           0.65/knee of output on cohesion) and `output_sigma` (the per-program reporting
           noise).  Both move the comparator's SIGNAL-TO-NOISE, which is the quantity
           that decides whether a 3.5-sigma rule can fire, and both are reported against
           a MEASURED gain axis rather than against the parameter: the sustained
           sigma-drop actually available in the lagging series.
  driver   (finding 3) `formation_driver = "degree"` -- the alternative the code retains
           "for the sensitivity check" and which was never run.  It wires the leading
           observable to mean repeat-tie degree instead of to mean component size, i.e.
           to a fast local quantity instead of to the slow order parameter.
  partial  (finding 4) `partial_obs_sigma` -- the assumed amplitude of the excess
           estimation error caused by seeing only part of a team.  The coverage floor is
           a direct function of this constant, so c*(partial_obs_sigma) is reported.
  size     (finding 5) N in {200, 800} with `programs_per_cycle` scaled to hold the
           mean-field control parameter K fixed (the Gate 1 recipe) AND `comp_ref`
           scaled with N, since a fixed comp_ref makes the leading observable's gain
           N-dependent by construction.

Every run is scored TWICE, under two pre-declared rules on the same series:
  * configs/detection_rule_recal.json      -- formation class only (the frozen rule)
  * configs/detection_rule_bothclass.json  -- formation AND escalation latency, threshold
                                              re-picked on the SAME 600-run null
                                              (calibrate_class.py), so the two rules are
                                              compared at a matched false-alarm target.
That is finding 11 answered by running the rule the manuscript describes rather than by
editing the sentence.

Run:  python run_remediation.py [--workers N] [--reps 30] [--gain-reps 40] [--quick]
Writes: aggregates/remediation_runs.json
"""

from __future__ import annotations

import argparse
import itertools
import math

import numpy as np

from indicators import fill_gaps
from calibrate_recal import max_sustained_min
from model import Config, Simulation
from scoring import score_run, bootstrap_mean
from runner import pmap, seed_for, write_agg, load_rule

C_GRID = (0.25, 0.333, 0.5, 0.75, 1.0)
SCENARIOS = ("degrade", "control")
ARM = "mnar"                 # the primary consent arm
COLOUR = "white"             # the headline noise colour
SEVERITY = "mid"             # the headline severity
RULES = ("detection_rule_recal.json", "detection_rule_bothclass.json")
BASE_VARIANT = "knee0.30_sig0.100"

# A(t) is a linear ramp from a_high at degrade_start to a_low at n_cycles, so the cycle
# at which the latent parameter passes a given A is analytic.  Used only to place the
# diagnostic windows in which cohesion and the two observables are measured.
def cycle_at_a(cfg: Config, a: float) -> int:
    span = cfg.n_cycles - cfg.degrade_start
    return int(round(cfg.degrade_start
                     + span * (cfg.a_high - a) / (cfg.a_high - cfg.a_low)))


# ---------------------------------------------------------------- the variant table
def variants() -> list[dict]:
    out: list[dict] = []
    for knee in (0.135, 0.15, 0.20, 0.30, 0.45, 0.60):
        out.append({"variant": f"knee{knee:.2f}_sig0.100", "leg": "gain",
                    "overrides": {"output_knee": knee}})
    for sig in (0.05, 0.025, 0.0125):
        out.append({"variant": f"knee0.30_sig{sig:.4f}", "leg": "gain",
                    "overrides": {"output_sigma": sig}})
    out.append({"variant": "driver_degree", "leg": "driver",
                "overrides": {"formation_driver": "degree"}})
    for pos in (0.55, 2.2):
        out.append({"variant": f"partial_obs_sigma{pos:.2f}", "leg": "partial",
                    "overrides": {"partial_obs_sigma": pos}})
    # N sweep: programs_per_cycle scaled to hold K = 2*W*R/N fixed (Gate 1's SIZES
    # table), and comp_ref scaled with N so the observable's gain is not N-dependent
    # by construction.
    for n, p in ((200, 6), (800, 24)):
        out.append({"variant": f"N{n}", "leg": "size",
                    "overrides": {"n_agents": n, "programs_per_cycle": p,
                                  "comp_ref": 100.0 * n / 400.0}})
    return out


class _CohSim(Simulation):
    """Simulation with a read-only cohesion tap.  `_cohesion` is a pure function of the
    tie graph and consumes no RNG, so logging its value cannot perturb a run: this class
    reproduces `run_config` bit-for-bit and the verifier checks that against the
    committed E-A rows."""

    def __init__(self, cfg: Config):
        super().__init__(cfg)
        self._coh: list[tuple[int, float]] = []
        self._t = -10 ** 6

    def step(self, t: int) -> None:            # noqa: D102
        self._t = t
        super().step(t)

    def _cohesion(self, team):                 # noqa: D102
        v = super()._cohesion(team)
        if self._t >= 0:
            self._coh.append((self._t, v))
        return v


def _window_mean(x: np.ndarray, lo: int, hi: int) -> float:
    seg = np.asarray(x, float)[lo:hi]
    seg = seg[np.isfinite(seg)]
    return float(seg.mean()) if seg.size else float("nan")


def _lag_diagnostics(series: dict[str, np.ndarray], rule: dict) -> dict:
    """The comparator's MEASURED gain: how many baseline sigmas of drop the lagging
    series actually offers, which is the quantity a k_sd rule is compared against.
    Computed with the frozen lagging rule's own baseline and gap-filling, so it is the
    same number the alarm loop sees."""
    G = rule["lagging"]
    base_end = int(G["baseline_end"])
    y = fill_gaps(np.asarray(series["output_close"], float), int(G["max_gap"]))
    base = y[:base_end]
    base = base[~np.isnan(base)]
    mu = float(base.mean())
    sd = float(base.std(ddof=1))
    drop = (mu - y) / (sd if sd > 1e-12 else 1e-12)
    tail = drop[base_end:]
    tail = tail[np.isfinite(tail)]
    return {
        "lag_baseline_mean": mu, "lag_baseline_sd": sd,
        "lag_sigma_drop_max": float(tail.max()) if tail.size else float("nan"),
        "lag_sigma_drop_sustained": max_sustained_min(drop, int(G["sustain"]), base_end),
        "lag_k_sd_required": float(G["k_sd"]),
    }


def _job(job):
    variant, overrides, c, scen, rep = job
    cfg = Config(consent_c=c, consent_arm=ARM, eps_colour=COLOUR,
                 eps_severity=SEVERITY, scenario=scen,
                 seed=seed_for(f"ea|{c}|{ARM}", rep), **overrides)
    sim = _CohSim(cfg)
    out = sim.run()
    s = out.series()
    row = {"variant": variant, "c": c, "scen": scen, "rep": rep, "seed": cfg.seed,
           "n_cycles": cfg.n_cycles,
           "achieved_edge_frac": float(s["achieved_edge_frac"][-1]),
           "obs_cycles_with_data": int(np.sum(~np.isnan(s["obs_formation"]))),
           "obs_esc_cycles_with_data": int(np.sum(~np.isnan(s["obs_escalation"])))}
    for rf in RULES:
        rule = load_rule(rf)
        sc = score_run(s, rule, cfg.n_cycles)
        pre = "" if rf == RULES[0] else "bc_"
        for k in ("t_lead", "t_lag", "dt_censored", "dt_paired", "lead_fired",
                  "lag_fired", "lead_max_score"):
            row[pre + k] = sc[k]
    row.update(_lag_diagnostics(s, load_rule(RULES[0])))

    # ---- amplitude diagnostics: the gain ratio review finding 2 is about ----
    t48 = cfg.degrade_start                      # A = a_high, still on the upper branch
    t40 = cycle_at_a(cfg, 0.40)                  # A = 0.40, still pre-fold
    tfold = cycle_at_a(cfg, 0.37)                # A = A_c (empirical fold, Gate 1)
    coh = np.asarray([v for (t, v) in sim._coh], float)
    coh_t = np.asarray([t for (t, v) in sim._coh], int)
    knee = cfg.output_knee

    def coh_win(lo, hi):
        m = (coh_t >= lo) & (coh_t < hi)
        return float(coh[m].mean()) if m.any() else float("nan")

    row.update({
        "fold_cycle": tfold,
        "coh_mean_a048": coh_win(t48 - 20, t48),
        "coh_mean_a040": coh_win(t40 - 10, t40 + 10),
        "coh_mean_prefold": coh_win(0, tfold),
        "coh_max": float(coh.max()) if coh.size else float("nan"),
        "coh_frac_at_or_above_knee": float(np.mean(coh >= knee)) if coh.size else
        float("nan"),
        "output_knee": knee, "output_sigma": cfg.output_sigma,
        "formation_driver": cfg.formation_driver,
        "partial_obs_sigma": cfg.partial_obs_sigma,
        "n_agents": cfg.n_agents, "comp_ref": cfg.comp_ref,
        "out_mean_a048": _window_mean(s["output_close"], t48 - 20, t48),
        "out_mean_a040": _window_mean(s["output_close"], t40 - 10, t40 + 10),
        "out_mean_last30": _window_mean(s["output_close"], cfg.n_cycles - 30,
                                        cfg.n_cycles),
        "form_true_mean_a048": _window_mean(s["true_formation"], t48 - 20, t48),
        "form_true_mean_a040": _window_mean(s["true_formation"], t40 - 10, t40 + 10),
        "form_true_mean_last30": _window_mean(s["true_formation"], cfg.n_cycles - 30,
                                             cfg.n_cycles),
        "form_obs_mean_a048": _window_mean(s["obs_formation"], t48 - 20, t48),
        "form_obs_mean_last30": _window_mean(s["obs_formation"], cfg.n_cycles - 30,
                                             cfg.n_cycles),
    })
    return row


# ------------------------------------------------------------------- reduction
def _iso(cs: list[float], ys: list[float]):
    for i in range(len(cs) - 1):
        if (ys[i] <= 0.0) != (ys[i + 1] <= 0.0):
            return cs[i] + (cs[i + 1] - cs[i]) * (0.0 - ys[i]) / (ys[i + 1] - ys[i])
    return None


def _rel(a: float, b: float) -> float:
    """Relative swing between two window means, as a fraction of the earlier one."""
    if not (np.isfinite(a) and np.isfinite(b)) or abs(a) < 1e-12:
        return float("nan")
    return (a - b) / a


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--gain-reps", type=int, default=40)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default="remediation_runs.json")
    args = ap.parse_args()

    vs = variants()
    if args.quick:
        vs = [v for v in vs if v["variant"] in (BASE_VARIANT, "driver_degree")]
    reps_of = {v["variant"]: (args.gain_reps if v["leg"] == "gain" else args.reps)
               for v in vs}
    if args.quick:
        reps_of = {k: 4 for k in reps_of}
    cs = C_GRID[:2] if args.quick else C_GRID

    jobs = [(v["variant"], v["overrides"], c, scen, r)
            for v in vs
            for c, scen in itertools.product(cs, SCENARIOS)
            for r in range(reps_of[v["variant"]])]
    rows = pmap(_job, jobs, args.workers, "remediation sweep")

    by_variant = {v["variant"]: v for v in vs}
    cells, gain_table, iso_table = [], [], {}
    for v in vs:
        name = v["variant"]
        vrows = [r for r in rows if r["variant"] == name]
        per_c = []
        for c in cs:
            deg = [r for r in vrows if r["c"] == c and r["scen"] == "degrade"]
            ctl = [r for r in vrows if r["c"] == c and r["scen"] == "control"]
            if not deg or not ctl:
                continue
            dt = bootstrap_mean([r["dt_censored"] for r in deg])
            dtp = bootstrap_mean([r["dt_paired"] for r in deg])
            bdt = bootstrap_mean([r["bc_dt_censored"] for r in deg])
            cell = {
                "variant": name, "leg": v["leg"], "c": c, "n_reps": len(deg),
                "dt_censored_mean": dt["mean"], "dt_censored_se": dt["se"],
                "dt_censored_lo95": dt["lo95"], "dt_censored_hi95": dt["hi95"],
                "dt_paired_mean": dtp["mean"], "dt_paired_se": dtp["se"],
                "dt_paired_n": dtp["n"],
                "lead_detect_rate": float(np.mean([r["lead_fired"] for r in deg])),
                "lag_detect_rate": float(np.mean([r["lag_fired"] for r in deg])),
                "fpr_leading": float(np.mean([r["lead_fired"] for r in ctl])),
                "fpr_lagging": float(np.mean([r["lag_fired"] for r in ctl])),
                "bc_dt_censored_mean": bdt["mean"], "bc_dt_censored_se": bdt["se"],
                "bc_lead_detect_rate": float(np.mean([r["bc_lead_fired"]
                                                      for r in deg])),
                "bc_fpr_leading": float(np.mean([r["bc_lead_fired"] for r in ctl])),
                "achieved_edge_frac_mean": float(np.mean(
                    [r["achieved_edge_frac"] for r in deg])),
                "lag_sigma_drop_sustained_mean": float(np.mean(
                    [r["lag_sigma_drop_sustained"] for r in deg])),
                "lag_sigma_drop_max_mean": float(np.mean(
                    [r["lag_sigma_drop_max"] for r in deg])),
                "lead_max_score_mean": float(np.mean(
                    [r["lead_max_score"] for r in deg])),
                "coh_mean_a048": float(np.mean([r["coh_mean_a048"] for r in deg])),
                "coh_mean_a040": float(np.mean([r["coh_mean_a040"] for r in deg])),
                "coh_frac_at_or_above_knee": float(np.mean(
                    [r["coh_frac_at_or_above_knee"] for r in deg])),
                "out_rel_swing_a048_to_a040": _rel(
                    float(np.mean([r["out_mean_a048"] for r in deg])),
                    float(np.mean([r["out_mean_a040"] for r in deg]))),
                "out_rel_swing_a048_to_end": _rel(
                    float(np.mean([r["out_mean_a048"] for r in deg])),
                    float(np.mean([r["out_mean_last30"] for r in deg]))),
                "form_true_rel_swing_a048_to_end": _rel(
                    float(np.mean([r["form_true_mean_a048"] for r in deg])),
                    float(np.mean([r["form_true_mean_last30"] for r in deg]))),
            }
            cells.append(cell)
            per_c.append(cell)
        per_c.sort(key=lambda x: x["c"])
        if len(per_c) >= 2:
            entry = {"variant": name, "leg": v["leg"],
                     "c": [x["c"] for x in per_c],
                     "dt": [x["dt_censored_mean"] for x in per_c],
                     "dt_se": [x["dt_censored_se"] for x in per_c],
                     "n_reps": [x["n_reps"] for x in per_c]}
            entry["c_at_dt_zero"] = _iso(entry["c"], entry["dt"])
            entry["c_at_dt_zero_band"] = [
                _iso(entry["c"], [x["dt_censored_mean"] + s * x["dt_censored_se"]
                                  for x in per_c]) for s in (-1.0, 1.0)]
            iso_table[name] = entry

        # ---- the measured comparator-gain axis, pooled over c ---------------------
        # The lagging metric is computed from program-close output, which no consent
        # dial touches, so its gain is the same at every c up to seed noise: pooling
        # over c is the tightest estimate of it and is reported as such.
        deg_all = [r for r in vrows if r["scen"] == "degrade"]
        top = [r for r in deg_all if r["c"] == max(cs)]
        gain_table.append({
            "variant": name, "leg": v["leg"],
            "overrides": v["overrides"],
            "n_degrade_runs": len(deg_all),
            "comparator_gain_sigma_drop": float(np.mean(
                [r["lag_sigma_drop_sustained"] for r in deg_all])),
            "comparator_gain_sigma_drop_se": float(
                np.std([r["lag_sigma_drop_sustained"] for r in deg_all], ddof=1)
                / math.sqrt(len(deg_all))),
            "comparator_gain_sigma_drop_max": float(np.mean(
                [r["lag_sigma_drop_max"] for r in deg_all])),
            "k_sd_required": float(load_rule(RULES[0])["lagging"]["k_sd"]),
            "comparator_can_fire": bool(np.mean(
                [r["lag_sigma_drop_sustained"] for r in deg_all])
                >= float(load_rule(RULES[0])["lagging"]["k_sd"])),
            "lag_detect_rate_pooled": float(np.mean([r["lag_fired"]
                                                     for r in deg_all])),
            "leading_gain_max_z_at_full_coverage": (
                float(np.mean([r["lead_max_score"] for r in top])) if top
                else float("nan")),
            "sensitivity_ratio_lead_over_lag_at_full_coverage": (
                float(np.mean([r["lead_max_score"] for r in top])
                      / np.mean([r["lag_sigma_drop_sustained"] for r in top]))
                if top else float("nan")),
            "dt_at_full_coverage": next(
                (x["dt_censored_mean"] for x in cells
                 if x["variant"] == name and x["c"] == max(cs)), None),
            "c_at_dt_zero": iso_table.get(name, {}).get("c_at_dt_zero"),
            "out_rel_swing_a048_to_a040": next(
                (x["out_rel_swing_a048_to_a040"] for x in cells
                 if x["variant"] == name and x["c"] == max(cs)), None),
            "coh_frac_at_or_above_knee": next(
                (x["coh_frac_at_or_above_knee"] for x in cells
                 if x["variant"] == name and x["c"] == max(cs)), None),
        })

    write_agg(args.out, {
        "_experiment": ("remediation campaign 2026-08-05 -- comparator-gain, "
                        "formation-driver, partial-observation-amplitude and "
                        "system-size sensitivity sweeps around the frozen headline "
                        "cell (review findings 2, 3, 4, 5, 11)"),
        "design": {
            "c_grid": list(cs), "arm": ARM, "colour": COLOUR, "severity": SEVERITY,
            "scenarios": list(SCENARIOS),
            "reps_default": args.reps, "reps_gain_leg": args.gain_reps,
            "base_variant": BASE_VARIANT,
            "rule_file": RULES[0], "rule_file_bothclass": RULES[1],
            "seed_construction": ("seed_for('ea|{c}|mnar', rep) -- the E-A seeds, so "
                                  "every variant is matched seed-for-seed to the "
                                  "headline run and the base variant reproduces "
                                  "ea_headline_recal.json where they overlap"),
            "n_runs": len(rows),
            "variants": [{"variant": v["variant"], "leg": v["leg"],
                          "overrides": v["overrides"]} for v in vs],
        },
        "detection_rule": load_rule(RULES[0]),
        "detection_rule_bothclass": load_rule(RULES[1]),
        "cells": cells,
        "comparator_gain_table": sorted(
            gain_table, key=lambda x: -x["comparator_gain_sigma_drop"]),
        "iso_line_dt_zero": iso_table,
        "per_run_rows": rows,
    })

    print("\n  variant                 gain(sigma)  lagDR   dt@c=1   c*(dt=0)")
    for x in sorted(gain_table, key=lambda y: -y["comparator_gain_sigma_drop"]):
        cstar = x["c_at_dt_zero"]
        print(f"  {x['variant']:22s} {x['comparator_gain_sigma_drop']:8.2f} "
              f"{x['lag_detect_rate_pooled']:7.2f} "
              f"{(x['dt_at_full_coverage'] or float('nan')):+8.1f}   "
              f"{'none' if cstar is None else f'{cstar:.3f}'}")


if __name__ == "__main__":
    main()
