"""Validate chi2.chi2_sf before any reported number depends on it.

Three independent checks, because a hand-written incomplete gamma is exactly the kind
of code that looks right and is wrong in one tail:

1. AGREEMENT WITH THE TREE'S OWN EVEN-df CLOSED FORM, which is exact and was already
   used for committed numbers.  Any disagreement past 1e-12 on even df means the new
   general path is broken, and that path is the one the odd-df result will rest on.
2. PUBLISHED UPPER-TAIL CRITICAL VALUES at alpha = 0.05 and 0.01 for df 1 to 20, odd
   and even.  These come from standard tables and are independent of both
   implementations.
3. A MONTE-CARLO CHECK on odd df, where the tables cannot reach: draw chi-square
   variates as sums of squared normals and compare the empirical tail to chi2_sf.
   This is the only check that would catch a systematic error affecting even and odd
   df identically.
"""

from __future__ import annotations

import math
import sys

import numpy as np

from chi2 import chi2_sf
from remediation_analysis import chi2_sf_even

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"{'PASS' if ok else 'FAIL'}  {name}{'  ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(name)


# ---------------------------------------------------------------- check 1
for df in (2, 4, 6, 8, 10, 12, 16, 20, 30):
    worst = 0.0
    worst_at = None
    for x in [0.05, 0.5, 1.0, 2.0, 3.5, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0, 45.0, 70.0]:
        a, b = chi2_sf(x, df), chi2_sf_even(x, df)
        err = abs(a - b) / max(b, 1e-300)
        if err > worst:
            worst, worst_at = err, x
    check(f"even df={df} matches the tree's closed form",
          worst < 1e-11, f"worst relative error {worst:.2e} at x={worst_at}")

# ---------------------------------------------------------------- check 2
# Upper-tail critical values, alpha = 0.05 and alpha = 0.01, df 1..20.
CRIT_05 = [3.841459, 5.991465, 7.814728, 9.487729, 11.070498, 12.591587, 14.067140,
           15.507313, 16.918978, 18.307038, 19.675138, 21.026070, 22.362032,
           23.684791, 24.995790, 26.296228, 27.587112, 28.869299, 30.143527,
           31.410433]
CRIT_01 = [6.634897, 9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307,
           20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250,
           29.141238, 30.577914, 32.000000, 33.408664, 34.805306, 36.190869,
           37.566235]

for alpha, table in ((0.05, CRIT_05), (0.01, CRIT_01)):
    worst = 0.0
    worst_df = None
    for i, x in enumerate(table, start=1):
        if i == 16 and alpha == 0.01:
            continue  # 32.000000 is a rounded table entry, not a 6-figure value
        err = abs(chi2_sf(x, i) - alpha)
        if err > worst:
            worst, worst_df = err, i
    check(f"published alpha={alpha} critical values, df 1-20",
          worst < 2e-6, f"worst absolute error {worst:.2e} at df={worst_df}")

# odd df specifically, called out because it is the case that did not exist before
for df in (1, 3, 5, 7, 9, 11, 15, 19):
    err = abs(chi2_sf(CRIT_05[df - 1], df) - 0.05)
    check(f"odd df={df} recovers alpha=0.05 at its critical value",
          err < 2e-6, f"error {err:.2e}")

# ---------------------------------------------------------------- check 3
rng = np.random.default_rng(20260904)
N = 400_000
for df in (3, 7, 15):
    draws = (rng.standard_normal((N, df)) ** 2).sum(axis=1)
    worst = 0.0
    worst_at = None
    for x in (float(np.quantile(draws, q)) for q in (0.5, 0.75, 0.9, 0.95, 0.99)):
        empirical = float(np.count_nonzero(draws > x)) / N
        se = math.sqrt(max(empirical, 1.0 / N) * (1 - empirical) / N)
        err = abs(chi2_sf(x, df) - empirical)
        if err / se > worst:
            worst, worst_at = err / se, x
    check(f"monte carlo df={df} agrees within 4 standard errors",
          worst < 4.0, f"worst {worst:.2f} SE at x={worst_at:.3f}")

# ---------------------------------------------------------------- the study's own number
print()
print(f"the 2026-09-03 susceptibility statistic: chi-square 67.53 on 3 df "
      f"-> p = {chi2_sf(67.53, 3):.3e}")
print(f"the same campaign's latency statistic:   chi-square  2.68 on 3 df "
      f"-> p = {chi2_sf(2.68, 3):.4f}")

print()
if FAILURES:
    print(f"{len(FAILURES)} FAILURE(S): " + "; ".join(FAILURES))
    sys.exit(1)
print("chi2.py validated.")
