"""Recompute the false-alarm MATCHING tests with the clustering handled.

The heterogeneity recheck (cluster_recheck.py) showed the three observed-subgraph detectors
give one verdict per (coverage, arm, replicate) group, replicated nine times by the noise
grid, so every test built on runs has a design effect of nine.  That applies to the MATCHING
tests as well, and those carry more weight in the campaign's conclusion: "all 3 of 3
observables are unmatched and no lead comparison from this campaign is reportable" rests on
two-proportion z-tests against the scalar comparator, pooled and per coverage cell, computed
over runs.

A two-proportion z-test on nine copies of each observation understates both standard errors by
a factor of three, so it rejects far too readily.  This script recomputes:

  1. the pooled false-alarm rate per detector, with a cluster-corrected Wilson interval,
  2. the pooled two-proportion test against the comparator at the group level,
  3. the per-cell tests at the group level, and how many cells fail.

The point of the exercise is that "unmatched" is the finding that suppressed every lead-time
result in the campaign.  If the detectors are in fact matched, those results come back.

Run:  python cluster_matching.py
"""

from __future__ import annotations

import collections
import json
import math
import os

from runner import HERE

FINE = "coverage_finegrid_block_b_20260904.json"
COARSE = "three_observable_block_b_results_20260903.json"
OBSERVABLES = ("obs_giant_frac", "obs_mean_degree", "obs_susceptibility")
Z95 = 1.959963985
ALPHA = 0.05


def phi(z):
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def ztest(k1, n1, k2, n2):
    """Two-proportion z-test.  k may be fractional after collapsing to groups."""
    if n1 <= 0 or n2 <= 0:
        return {"z": float("nan"), "p": float("nan")}
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2)) if 0 < p < 1 else 0.0
    if se == 0.0:
        return {"z": 0.0, "p": 1.0}
    z = (k1 / n1 - k2 / n2) / se
    return {"z": z, "p": 2.0 * (1.0 - phi(abs(z)))}


def wilson(k, n):
    if n <= 0:
        return [float("nan"), float("nan")]
    p = k / n
    c = (p + Z95 * Z95 / (2 * n)) / (1 + Z95 * Z95 / n)
    h = Z95 * math.sqrt(max(p * (1 - p), 0.0) / n + Z95 * Z95 / (4 * n * n)) / (1 + Z95 * Z95 / n)
    return [c - h, c + h]


def load(name):
    with open(os.path.join(HERE, "aggregates", name), encoding="utf-8") as f:
        return json.load(f)


def analyse(label, agg_name):
    d = load(agg_name)
    rows = [r for r in d["per_run_rows"] if r["scen"] == "control"]
    groups = collections.defaultdict(list)
    for r in rows:
        groups[(r["c"], r["arm"], r["rep"])].append(r)
    coverages = sorted({r["c"] for r in rows})

    def grp(det, c=None):
        """(sum of group fired-fractions, number of groups), optionally within one cell."""
        gs = [v for k, v in groups.items() if c is None or k[0] == c]
        return (sum(sum(x[det]["lead_fired"] for x in v) / len(v) for v in gs), len(gs))

    def run(det, c=None):
        rr = [r for r in rows if c is None or r["c"] == c]
        return (sum(r[det]["lead_fired"] for r in rr), len(rr))

    print(f"\n=== {label} ===")
    print(f"{len(rows)} control runs -> {len(groups)} independent groups, "
          f"{len(coverages)} coverage cells\n")

    lk_r, ln_r = run("latency")
    lk_g, ln_g = grp("latency")
    print(f"comparator pooled false-alarm rate: {100 * lk_r / ln_r:.2f}% "
          f"({lk_r:.0f}/{ln_r} runs) = {100 * lk_g / ln_g:.2f}% "
          f"({lk_g:.2f}/{ln_g} groups); cluster-corrected Wilson 95% "
          f"[{100 * wilson(lk_g, ln_g)[0]:.2f}%, {100 * wilson(lk_g, ln_g)[1]:.2f}%]")
    print()
    header = (f"{'detector':22s} {'rate':>7s}  {'RUN-level p':>12s} {'verdict':>10s}   "
              f"{'GROUP-level p':>14s} {'verdict':>10s}   {'cells failing':>14s}")
    print(header)
    out = {}
    for det in OBSERVABLES:
        kr, nr = run(det)
        kg, ng = grp(det)
        a = ztest(kr, nr, lk_r, ln_r)
        b = ztest(kg, ng, lk_g, ln_g)
        fail_run = fail_grp = 0
        cellrows = []
        for c in coverages:
            k1, n1 = run(det, c)
            k2, n2 = run("latency", c)
            g1, m1 = grp(det, c)
            g2, m2 = grp("latency", c)
            ta, tb = ztest(k1, n1, k2, n2), ztest(g1, m1, g2, m2)
            fail_run += ta["p"] < ALPHA
            fail_grp += tb["p"] < ALPHA
            cellrows.append({"c": c, "run_p": ta["p"], "group_p": tb["p"],
                             "detector_groups": [g1, m1], "latency_groups": [g2, m2]})
        va = "unmatched" if a["p"] < ALPHA else "MATCHED"
        vb = "unmatched" if b["p"] < ALPHA else "MATCHED"
        print(f"{det:22s} {100 * kg / ng:6.2f}%  {a['p']:12.2e} {va:>10s}   "
              f"{b['p']:14.4f} {vb:>10s}   {fail_grp:2d} of {len(coverages)} "
              f"(was {fail_run})")
        out[det] = {"pooled_rate_groups": kg / ng,
                    "wilson95_groups": wilson(kg, ng),
                    "pooled_run_level": a, "pooled_group_level": b,
                    "matched_at_group_level": b["p"] >= ALPHA,
                    "cells_failing_group_level": fail_grp,
                    "cells_failing_run_level": fail_run,
                    "n_cells": len(coverages), "per_cell": cellrows}
    return out


if __name__ == "__main__":
    result = {"_analysis": ("false-alarm matching tests recomputed at the (coverage, arm, "
                            "replicate) group level; run-level results shown alongside"),
              "coarse_4_cell_20260903": analyse("4-cell grid, 2026-09-03 campaign", COARSE),
              "fine_16_cell_20260904": analyse("16-cell grid, 2026-09-04 rerun", FINE)}
    path = os.path.join(HERE, "aggregates", "cluster_matching_20260904.json")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=1)
        f.write("\n")
    print(f"\nwritten: aggregates/cluster_matching_20260904.json")
