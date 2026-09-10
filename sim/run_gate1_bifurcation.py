"""
GATE 1 -- the bifurcation obligation.

No early-warning signal may be reported before this gate passes.  Two obligations,
separately cited in the spec:
  (a) show the transition is bifurcation-like   (Boettiger & Hastings 2013: no EWS
      exist for noise-induced transitions);
  (b) show the indicator does not fire in the non-transitioning control arm
      (Looker, Rock & Dyson 2026: rising variance and return time can arise from
      non-normal geometry with no bifurcation) -- (b) is discharged by the paired
      control arm of E-A, not here.

Four independent lines of evidence for (a), all on the CUMULATIVE REPEAT-TIE GRAPH:

  E1  Equilibrium branch diagram S(A) at four system sizes, with the number of
      programs per cycle scaled with N so the mean-field control parameter
      K = 2 * W * P * <pairs> / N is held FIXED.  A fold shows a DISCONTINUOUS jump
      whose size does not vanish with N and a transition that SHARPENS with N; a
      continuous (percolation-only) transition shows S going to zero smoothly.

  E2  Hysteresis: the same A grid entered from the upper branch (burn-in at high A)
      and from the lower branch (burn-in at low A).  A non-zero-width bistable
      region is a saddle-node signature and cannot be produced by a continuous
      transition or by noise alone.

  E3  Perturb-and-relax: at equilibrium, delete a fraction of repeat ties and fit
      the exponential recovery rate of the mean degree.  The recovery rate must go
      to ZERO as A approaches the fold -- this is critical slowing down measured
      directly on the latent structure, independent of any observed indicator.

  E4  Equilibrium fluctuation and lag-1 autocorrelation of the order parameter
      peaking at the transition.

Run:  python run_gate1_bifurcation.py [--workers N] [--quick]
Writes: aggregates/gate1_bifurcation.json
"""

from __future__ import annotations

import argparse
import json
import os
import numpy as np

from model import Config, run_config
from runner import AGG, pmap, seed_for, write_agg

# N is swept with P scaled so K stays fixed at the N=400, P=12 value.
SIZES = ((200, 6), (400, 12), (800, 24), (1600, 48))
A_GRID = (0.48, 0.44, 0.42, 0.40, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.31, 0.28)
BRANCHES = {"upper": 0.90, "lower": 0.18}
N_SEEDS = 8
TAIL = 60            # cycles of the tail averaged as the equilibrium estimate

SHOCK_A = (0.48, 0.44, 0.42, 0.40, 0.39, 0.38, 0.37)
SHOCK_FRAC = 0.35
SHOCK_AT = 40
SHOCK_SEEDS = 12

# E5: long runs near the fold.  The relaxation time within the bistable window
# exceeds 300 cycles (see E3), so the branch diagram of E1 measures a slow transient
# rather than the stationary state.  These runs are long enough for each seed to have
# settled onto ONE branch, which is what makes the end-state distribution bimodal.
LONG_A = (0.40, 0.38, 0.37, 0.36, 0.35, 0.34)
LONG_SIZES = ((400, 12), (800, 24))
LONG_CYCLES = 700
LONG_SEEDS = 12


def mean_field_fold(cfg: Config) -> dict:
    """Fixed points of the substrate's own update rule, z = K*A*g(z), with
    K = 2*W*R/N and R the expected pair-events per cycle.  A saddle-node (fold) is a
    tangency: F(z) = K*A*g(z) - z has a double root.  This is an analytic property of
    the model as written, independent of any simulation."""
    m = np.asarray(cfg.team_size_choices, float)
    R = cfg.programs_per_cycle * float(np.mean(m * (m - 1) / 2.0))
    K = 2.0 * cfg.tie_window * R / cfg.n_agents
    z = np.linspace(1e-4, 8.0, 400_001)
    h2 = cfg.reengage_hill ** 2
    g = cfg.reengage_theta + (1.0 - cfg.reengage_theta) * z * z / (h2 + z * z)
    out = {"K_mean_field": K, "pair_events_per_cycle": R, "a_scan": {}}
    fold_a, fold_z = None, None
    for a in np.arange(0.20, 0.60, 0.0005):
        F = K * a * g - z
        sign = np.signbit(F)
        roots = z[:-1][sign[:-1] != sign[1:]]
        n_fp = int(roots.size)
        if n_fp >= 3 and fold_a is None:
            pass
        out["a_scan"][round(float(a), 4)] = n_fp
    # the fold is the largest A at which the upper pair of fixed points has merged,
    # i.e. the boundary between 1 fixed point (low only) and 3 fixed points
    keys = sorted(out["a_scan"])
    for a1, a2 in zip(keys[:-1], keys[1:]):
        if out["a_scan"][a1] == 1 and out["a_scan"][a2] >= 3:
            fold_a = 0.5 * (a1 + a2)
            F = K * fold_a * g - z
            fold_z = float(z[np.argmin(np.abs(F))])
            break
    n_bistable = sum(1 for a in keys if out["a_scan"][a] >= 3)
    out["a_scan"] = {k: v for k, v in list(out["a_scan"].items())[::40]}  # thin
    out["fold_a"] = fold_a
    out["fold_z"] = fold_z
    out["bistable_a_count"] = n_bistable
    out["bistable_a_width"] = n_bistable * 0.0005
    out["has_saddle_node"] = bool(fold_a is not None and n_bistable > 0)
    out["note"] = ("duplicate ties (a pair re-tied while already in the window) are "
                   "ignored in this estimate, so the empirical fold sits ABOVE this A")
    return out


def _long_job(job):
    n, p, a, rep = job
    cfg = Config(n_agents=n, programs_per_cycle=p, scenario="control", a_control=a,
                 burn_in_a=0.90, n_cycles=LONG_CYCLES,
                 seed=seed_for(f"gate1long|{n}|{a}", rep))
    s = run_config(cfg).series()
    return {"n": n, "a": a, "rep": rep,
            "z_end": float(s["mean_degree"][-60:].mean()),
            "S_end": float(s["giant_frac"][-60:].mean())}


def _equil_job(job):
    n, p, a, branch, rep, n_cycles = job
    cfg = Config(n_agents=n, programs_per_cycle=p, scenario="control", a_control=a,
                 burn_in_a=BRANCHES[branch], n_cycles=n_cycles,
                 seed=seed_for(f"gate1|{n}|{a}|{branch}", rep))
    s = run_config(cfg).series()
    S = s["giant_frac"][-TAIL:]
    z = s["mean_degree"][-TAIL:]
    dz = z - z.mean()
    ac = (float(np.corrcoef(dz[:-1], dz[1:])[0, 1])
          if dz.size > 3 and dz.std() > 1e-12 else float("nan"))
    return {"n": n, "a": a, "branch": branch, "rep": rep,
            "S": float(S.mean()), "S_sd": float(S.std()),
            "z": float(z.mean()), "z_sd": float(z.std()), "z_ac1": ac}


def _shock_job(job):
    a, rep, n_cycles = job
    cfg = Config(scenario="control", a_control=a, burn_in_a=0.90,
                 n_cycles=n_cycles, shock_cycle=SHOCK_AT, shock_frac=SHOCK_FRAC,
                 seed=seed_for(f"gate1shock|{a}", rep))
    s = run_config(cfg).series()
    z = s["mean_degree"]
    z_star = float(z[max(0, SHOCK_AT - 25):SHOCK_AT].mean())
    tail = z[-30:]
    z_end = float(tail.mean())
    # fit log-deviation decay over the recovery leg
    seg = z[SHOCK_AT:SHOCK_AT + 60]
    dev = np.abs(z_star - seg)
    ok = dev > 1e-3
    lam = float("nan")
    if ok.sum() >= 8:
        x = np.arange(seg.size, dtype=float)[ok]
        y = np.log(dev[ok])
        xm, ym = x.mean(), y.mean()
        sxx = float((x - xm) @ (x - xm))
        if sxx > 1e-9:
            lam = -float((x - xm) @ (y - ym)) / sxx
    return {"a": a, "rep": rep, "z_star": z_star, "z_min": float(seg.min()),
            "z_end": z_end, "recovery_rate": lam,
            "recovered": bool(z_end > 0.6 * z_star)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse the equilibrium/relaxation rows already committed")
    args = ap.parse_args()
    n_cycles = 90 if args.quick else 140
    sizes = SIZES[:2] if args.quick else SIZES
    seeds = 3 if args.quick else N_SEEDS

    mf = mean_field_fold(Config())
    print(f"  E0 mean-field: saddle-node={mf['has_saddle_node']} "
          f"fold_a={mf['fold_a']} fold_z={mf['fold_z']} "
          f"bistable width in A={mf['bistable_a_width']:.4f}", flush=True)

    prev = None
    if args.reuse and os.path.exists(os.path.join(AGG, "gate1_bifurcation.json")):
        with open(os.path.join(AGG, "gate1_bifurcation.json")) as f:
            prev = json.load(f)
    if prev and prev.get("per_seed_rows"):
        rows = prev["per_seed_rows"]
        shock_prev = prev.get("relaxation_table")
        print(f"  reusing {len(rows)} equilibrium rows from the previous aggregate",
              flush=True)
    else:
        shock_prev = None
        jobs = [(n, p, a, b, r, n_cycles)
                for (n, p) in sizes for a in A_GRID for b in BRANCHES
                for r in range(seeds)]
        rows = pmap(_equil_job, jobs, args.workers, "E1/E2/E4 equilibrium sweep")

    # aggregate per (n, a, branch)
    cells: dict[tuple, list] = {}
    for r in rows:
        cells.setdefault((r["n"], r["a"], r["branch"]), []).append(r)
    branch_tab = []
    for (n, a, b), rs in sorted(cells.items()):
        branch_tab.append({
            "n": n, "a": a, "branch": b, "n_seeds": len(rs),
            "S_mean": float(np.mean([x["S"] for x in rs])),
            "S_between_seed_sd": float(np.std([x["S"] for x in rs])),
            "z_mean": float(np.mean([x["z"] for x in rs])),
            "z_temporal_sd": float(np.mean([x["z_sd"] for x in rs])),
            "z_ac1": float(np.nanmean([x["z_ac1"] for x in rs])),
        })

    # ---- fold location, discontinuity, hysteresis width, per N ---------------
    # The order parameter of the FOLD is the mean degree z of the cumulative
    # repeat-tie graph (the giant-component fraction S is a percolation function of
    # z and is already small at the fold, so the jump must be read on z).
    finite = {}
    for n, _p in sizes:
        up = {r["a"]: r for r in branch_tab if r["n"] == n and r["branch"] == "upper"}
        lo = {r["a"]: r for r in branch_tab if r["n"] == n and r["branch"] == "lower"}
        asc = sorted(up, reverse=True)
        z_hi = up[asc[0]]["z_mean"]
        z_lo = up[asc[-1]]["z_mean"]
        rng_z = max(1e-9, z_hi - z_lo)
        drops = [(up[a1]["z_mean"] - up[a2]["z_mean"], a1, a2)
                 for a1, a2 in zip(asc[:-1], asc[1:])]
        jump_z, a_before, a_after = max(drops)
        # finite-size signature of a DISCONTINUOUS transition: at the A where the
        # ensemble mean sits mid-way between branches, individual seeds do not --
        # they are on one branch or the other, so the middle of the range is depleted.
        frac = {a: (up[a]["z_mean"] - z_lo) / rng_z for a in asc}
        a_mid = min(asc, key=lambda a: abs(frac[a] - 0.5))
        seed_z = np.asarray([r["z"] for r in rows
                             if r["n"] == n and r["a"] == a_mid
                             and r["branch"] == "upper"], float)
        mid_z = 0.5 * (z_hi + z_lo)
        band = 0.20 * rng_z
        occupancy = (float(np.mean(np.abs(seed_z - mid_z) < band))
                     if seed_z.size else float("nan"))
        # transition width: A window over which the ensemble mean crosses 0.15->0.85
        inner = [a for a in asc if 0.15 < frac[a] < 0.85]
        width = (max(inner) - min(inner)) if len(inner) > 1 else 0.0
        bist = [a for a in asc if (up[a]["z_mean"] - lo[a]["z_mean"]) > 0.30]
        finite[str(n)] = {
            "fold_between_a": [a_before, a_after],
            "fold_a_midpoint": 0.5 * (a_before + a_after),
            "jump_in_z": jump_z,
            "jump_in_z_relative": jump_z / rng_z,
            "z_above_fold": up[a_before]["z_mean"],
            "z_below_fold": up[a_after]["z_mean"],
            "jump_in_S": up[a_before]["S_mean"] - up[a_after]["S_mean"],
            "a_midpoint_of_crossing": a_mid,
            "mid_band_seed_occupancy": occupancy,
            "transition_width_a": width,
            "bistable_a_values": bist,
            "bistable_width_a": (max(bist) - min(bist)) if len(bist) > 1 else 0.0,
            "max_upper_lower_gap_z": max(up[a]["z_mean"] - lo[a]["z_mean"]
                                         for a in asc),
            "peak_z_temporal_sd_at_a": max(
                (r["z_temporal_sd"], r["a"]) for r in branch_tab
                if r["n"] == n and r["branch"] == "upper")[1],
            "peak_z_ac1_at_a": max(
                (r["z_ac1"], r["a"]) for r in branch_tab
                if r["n"] == n and r["branch"] == "upper")[1],
        }

    # ---- E3 perturb-and-relax ----------------------------------------------
    if shock_prev:
        shock_tab = shock_prev
        print(f"  reusing {len(shock_tab)} relaxation rows", flush=True)
    else:
        sjobs = [(a, r, n_cycles) for a in SHOCK_A
                 for r in range(SHOCK_SEEDS if not args.quick else 4)]
        srows = pmap(_shock_job, sjobs, args.workers, "E3 perturb-and-relax")
        shock_tab = []
        for a in SHOCK_A:
            rs = [r for r in srows if r["a"] == a]
            lam = [r["recovery_rate"] for r in rs if np.isfinite(r["recovery_rate"])]
            shock_tab.append({
                "a": a, "n_seeds": len(rs),
                "recovery_rate_mean": float(np.mean(lam)) if lam else None,
                "recovery_rate_sd": float(np.std(lam)) if lam else None,
                "recovery_time_cycles": (float(1.0 / np.mean(lam))
                                         if lam and np.mean(lam) > 1e-6 else None),
                "recovered_frac": float(np.mean([r["recovered"] for r in rs])),
                "z_star_mean": float(np.mean([r["z_star"] for r in rs])),
            })

    # ---- E5 long-run end-state distribution near the fold -------------------
    lsizes = LONG_SIZES[:1] if args.quick else LONG_SIZES
    lseeds = 4 if args.quick else LONG_SEEDS
    ljobs = [(n, p, a, r) for (n, p) in lsizes for a in LONG_A
             for r in range(lseeds)]
    lrows = pmap(_long_job, ljobs, args.workers, "E5 long runs near the fold")
    long_tab, bimodal = [], {}
    for n, _p in lsizes:
        zs_by_a = {a: np.asarray([r["z_end"] for r in lrows
                                  if r["n"] == n and r["a"] == a], float)
                   for a in LONG_A}
        z_hi = float(zs_by_a[LONG_A[0]].mean())
        z_lo = float(zs_by_a[LONG_A[-1]].mean())
        rngz = max(1e-9, z_hi - z_lo)
        for a in LONG_A:
            zs = zs_by_a[a]
            long_tab.append({
                "n": n, "a": a, "n_seeds": int(zs.size),
                "z_end_mean": float(zs.mean()), "z_end_sd": float(zs.std()),
                "frac_on_upper_branch": float(np.mean(zs > z_lo + 0.6 * rngz)),
                "frac_on_lower_branch": float(np.mean(zs < z_lo + 0.2 * rngz)),
            })
        # crossing point: A whose mean sits nearest mid-range; middle-band occupancy
        fr = {a: (float(zs_by_a[a].mean()) - z_lo) / rngz for a in LONG_A}
        a_mid = min(LONG_A, key=lambda a: abs(fr[a] - 0.5))
        zs = zs_by_a[a_mid]
        mid = 0.5 * (z_hi + z_lo)
        occ = float(np.mean(np.abs(zs - mid) < 0.20 * rngz))
        # largest gap in the pooled end-state distribution over the fold window
        pooled = np.sort(np.concatenate([zs_by_a[a] for a in LONG_A]))
        gaps = np.diff(pooled)
        bimodal[str(n)] = {
            "a_crossing": a_mid, "mid_band_occupancy": occ,
            "z_upper_ref": z_hi, "z_lower_ref": z_lo,
            "largest_gap_in_pooled_endstates": float(gaps.max()) if gaps.size else 0.0,
            "largest_gap_relative": (float(gaps.max()) / rngz) if gaps.size else 0.0,
        }

    # ---- verdict ------------------------------------------------------------
    n_keys = sorted(finite, key=lambda s: int(s))
    rel_jumps = [finite[k]["jump_in_z_relative"] for k in n_keys]
    occ = [finite[k]["mid_band_seed_occupancy"] for k in n_keys]
    tw = [finite[k]["transition_width_a"] for k in n_keys]
    widths = [finite[k]["bistable_width_a"] for k in n_keys]
    rates = [r["recovery_rate_mean"] for r in shock_tab
             if r["recovery_rate_mean"] is not None]
    slowing = (len(rates) >= 3 and rates[0] > 0 and rates[-1] < 0.5 * rates[0])
    tests = {
        "T1_mean_field_saddle_node": {
            "pass": bool(mf["has_saddle_node"]),
            "fold_a_mean_field": mf["fold_a"],
            "fold_z_mean_field": mf["fold_z"],
            "bistable_a_width_mean_field": mf["bistable_a_width"],
            "criterion": "the substrate's own update rule z = K*A*g(z) has a "
                         "tangency (1 fixed point below the fold, 3 above) -- an "
                         "analytic saddle-node, not a noise-induced transition",
        },
        "T2_endstate_bimodality": {
            "pass": bool(all(v["mid_band_occupancy"] < 0.5
                             and v["largest_gap_relative"] > 0.25
                             for v in bimodal.values())),
            "by_N": bimodal,
            "criterion": "in runs long enough to outlast the divergent relaxation "
                         "time, end states are bimodal: fewer than half the seeds in "
                         "the middle 40% of the branch range at the crossing A, and a "
                         "gap of >25% of the range in the pooled end-state "
                         "distribution",
        },
        "T2b_transient_branch_diagram": {
            "pass": None,
            "note": ("DESCRIPTIVE ONLY, not a gate test. Within the bistable window "
                     "the relaxation time (E3) exceeds the run length, so the E1 "
                     "branch diagram measures a slow transient and its ensemble mean "
                     "interpolates the two branches. This is why the discontinuity "
                     "test moved to E5 long runs (T2)."),
            "relative_jump_in_z_by_N": dict(zip(n_keys, rel_jumps)),
            "mid_band_seed_occupancy_by_N": dict(zip(n_keys, occ)),
            "transition_width_a_by_N": dict(zip(n_keys, tw)),
        },
        "T3_hysteresis": {
            "pass": bool(min(widths) > 0.0),
            "bistable_width_by_N": dict(zip(n_keys, widths)),
            "max_upper_lower_gap_z_by_N": {k: finite[k]["max_upper_lower_gap_z"]
                                           for k in n_keys},
            "criterion": "non-zero bistable A-window (upper and lower branch z "
                         "differ by >0.30) at every N",
        },
        "T4_critical_slowing_down": {
            "pass": bool(slowing),
            "recovery_rate_by_a": {r["a"]: r["recovery_rate_mean"]
                                   for r in shock_tab},
            "recovery_time_cycles_by_a": {r["a"]: r["recovery_time_cycles"]
                                          for r in shock_tab},
            "criterion": "post-shock recovery rate of z falls by more than 2x as A "
                         "approaches the fold",
        },
    }
    graded = {k: v for k, v in tests.items() if v.get("pass") is not None}
    verdict = "PASS" if all(v["pass"] for v in graded.values()) else "FAIL"

    write_agg("gate1_bifurcation.json", {
        "_gate": "GATE 1 -- bifurcation obligation",
        "verdict": verdict,
        "tests": tests,
        "design": {"sizes_N_and_P": [list(x) for x in sizes], "a_grid": list(A_GRID),
                   "branches_burn_in_a": BRANCHES, "n_seeds": seeds,
                   "equilibrium_tail_cycles": TAIL, "n_cycles": n_cycles,
                   "shock_frac": SHOCK_FRAC, "shock_cycle": SHOCK_AT,
                   "long_run_a": list(LONG_A), "long_run_cycles": LONG_CYCLES,
                   "long_run_sizes": [list(x) for x in LONG_SIZES],
                   "long_run_seeds": LONG_SEEDS,
                   "K_held_fixed": "P scaled proportionally to N"},
        "finite_size": finite,
        "branch_table": branch_tab,
        "relaxation_table": shock_tab,
        "mean_field": mf,
        "long_run_table": long_tab,
        "per_seed_rows": rows,
        "long_run_rows": lrows,
    })
    print(f"GATE 1 verdict: {verdict}")
    for k, v in graded.items():
        print(f"  {k}: {'PASS' if v['pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
