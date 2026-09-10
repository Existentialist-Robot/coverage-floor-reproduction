"""
Early-warning indicators and the pre-declared detection rule.

Frozen-spec constraints honoured here:
  * VARIANCE FAMILY ONLY (spec section 2 / Yu, MacLaren & Masuda 2026):
    rolling variance and rolling lag-1 autocorrelation.  There is no covariance
    term and no pair-combination rule anywhere in this file.
  * Node/class mixing is permitted for variance-family signals (Masuda 2024): we
    compute each statistic SEPARATELY per signal class
    (formation latency, escalation latency) and mix at the level of the trend
    statistic, never by combining a pair of nodes into one statistic.
  * The node set is never optimised -- consent decides it.  Nothing in this file
    selects which agents to observe.
  * The detection rule is PRE-DECLARED in the JSON configs (see
    configs/detection_rule.json); this module only executes a rule handed to it.

Everything is stdlib + numpy (no scipy): Kendall tau, the Mann-Whitney/AUC
discrimination statistic and the logistic propensity fit are implemented here.
"""

from __future__ import annotations

import math
import warnings
from typing import Any, Sequence

import numpy as np


# --------------------------------------------------------------------------
# missing-data handling (pre-declared)
# --------------------------------------------------------------------------
def fill_gaps(x: np.ndarray, max_gap: int) -> np.ndarray:
    """Linear interpolation across gaps of at most `max_gap` cycles; longer gaps
    stay nan.  Leading/trailing nans stay nan.  Declared rule, not a tuning knob."""
    y = np.asarray(x, float).copy()
    n = y.size
    good = np.flatnonzero(~np.isnan(y))
    if good.size == 0:
        return y
    for a, b in zip(good[:-1], good[1:]):
        gap = b - a - 1
        if 0 < gap <= max_gap:
            y[a + 1:b] = np.linspace(y[a], y[b], gap + 2)[1:-1]
    return y


# --------------------------------------------------------------------------
# rolling variance-family statistics with in-window linear detrending
# --------------------------------------------------------------------------
def _detrend(seg: np.ndarray) -> np.ndarray:
    """In-window linear detrending, closed form (no polyfit -- this is the hot loop)."""
    m = ~np.isnan(seg)
    k = int(m.sum())
    if k < 3:
        return np.full_like(seg, np.nan)
    idx = np.arange(seg.size, dtype=float)
    x = idx[m]
    y = seg[m]
    xm = x.mean()
    ym = y.mean()
    dx = x - xm
    sxx = float(dx @ dx)
    slope = float(dx @ (y - ym)) / sxx if sxx > 1e-12 else 0.0
    return seg - (slope * (idx - xm) + ym)


def rolling_var_ac(x: np.ndarray, window: int, min_frac: float = 0.5,
                   max_gap: int = 3) -> tuple[np.ndarray, np.ndarray]:
    """Rolling variance and rolling lag-1 autocorrelation of the in-window
    linearly-detrended series.  Returns arrays aligned to the window END index;
    entries before `window-1` are nan."""
    y = fill_gaps(x, max_gap)
    n = y.size
    var = np.full(n, np.nan)
    ac = np.full(n, np.nan)
    need = int(math.ceil(min_frac * window))
    for end in range(window - 1, n):
        seg = y[end - window + 1:end + 1]
        if np.count_nonzero(~np.isnan(seg)) < need:
            continue
        d = _detrend(seg)
        m = ~np.isnan(d)
        dm = d[m]
        if dm.size < 4:
            continue
        var[end] = float(np.var(dm, ddof=1))
        # lag-1 AC on the detrended segment, using consecutive valid pairs only
        pair = m[:-1] & m[1:]
        if np.count_nonzero(pair) < 4:
            continue
        a = d[:-1][pair]
        b = d[1:][pair]
        sa, sb = a.std(), b.std()
        if sa < 1e-12 or sb < 1e-12:
            continue
        ac[end] = float(np.mean((a - a.mean()) * (b - b.mean())) / (sa * sb))
    return var, ac


# --------------------------------------------------------------------------
# Kendall tau (tau-b), no scipy
# --------------------------------------------------------------------------
def kendall_tau(x: Sequence[float], y: Sequence[float]) -> float:
    xa = np.asarray(x, float)
    ya = np.asarray(y, float)
    m = ~(np.isnan(xa) | np.isnan(ya))
    xa, ya = xa[m], ya[m]
    n = xa.size
    if n < 4:
        return float("nan")
    conc = disc = tx = ty = 0
    for i in range(n - 1):
        dx = xa[i + 1:] - xa[i]
        dy = ya[i + 1:] - ya[i]
        s = np.sign(dx) * np.sign(dy)
        conc += int(np.count_nonzero(s > 0))
        disc += int(np.count_nonzero(s < 0))
        tx += int(np.count_nonzero((dx == 0) & (dy != 0)))
        ty += int(np.count_nonzero((dy == 0) & (dx != 0)))
    d0 = math.sqrt((conc + disc + tx) * (conc + disc + ty))
    if d0 <= 0:
        return float("nan")
    return (conc - disc) / d0


def rolling_tau(stat: np.ndarray, span: int) -> np.ndarray:
    """Kendall tau of the trailing `span` values of a statistic series, aligned
    to the trailing index."""
    n = stat.size
    out = np.full(n, np.nan)
    t = np.arange(span, dtype=float)
    for end in range(span - 1, n):
        seg = stat[end - span + 1:end + 1]
        if np.count_nonzero(~np.isnan(seg)) < max(4, int(0.6 * span)):
            continue
        out[end] = kendall_tau(t, seg)
    return out


# --------------------------------------------------------------------------
# the pre-declared detection rule
# --------------------------------------------------------------------------
def class_trends(series: dict[str, np.ndarray], window: int, span: int,
                 classes: Sequence[str], min_frac: float = 0.5,
                 max_gap: int = 3) -> dict[str, Any]:
    """Per-class rolling variance / lag-1 AC and their Kendall-tau trend series,
    then the class mix (mean tau across classes -- variance-family mixing at the
    level of the trend statistic, never a covariance across a node pair)."""
    tv_list, ta_list, per_class = [], [], {}
    for cls in classes:
        var, ac = rolling_var_ac(series[cls], window, min_frac, max_gap)
        tv = rolling_tau(var, span)
        ta = rolling_tau(ac, span)
        tv_list.append(tv)
        ta_list.append(ta)
        per_class[cls] = {"var": var, "ac": ac, "tau_var": tv, "tau_ac": ta}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        tau_var = np.nanmean(np.vstack(tv_list), axis=0)
        tau_ac = np.nanmean(np.vstack(ta_list), axis=0)
    return {"tau_var": tau_var, "tau_ac": tau_ac, "per_class": per_class}


def baseline_z(series: dict[str, np.ndarray], window: int, classes: Sequence[str],
               baseline_end: int, min_frac: float = 0.5,
               max_gap: int = 3) -> dict[str, Any]:
    """Rolling variance and lag-1 AC standardised against their own BASELINE-PERIOD
    distribution, then mixed across signal classes.

    The baseline period is the pre-degradation stretch of the same run, which the
    instrument observes in-run -- so this is a rule an operator could actually
    execute, and it is the same form as the lagging comparator (a baseline-sigma
    rule), which is what makes the lead-time comparison fair.

    Mixing is across signal CLASSES of a variance-family statistic (Masuda 2024).
    No covariance across node pairs appears anywhere (Yu et al. 2026)."""
    zv_list, za_list, per_class = [], [], {}
    for cls in classes:
        var, ac = rolling_var_ac(series[cls], window, min_frac, max_gap)
        zs = {}
        for name, x in (("var", var), ("ac", ac)):
            base = x[window - 1:baseline_end]
            base = base[~np.isnan(base)]
            if base.size < 8:
                zs[name] = np.full_like(x, np.nan)
                continue
            mu = float(base.mean())
            sd = float(base.std(ddof=1))
            zs[name] = (x - mu) / (sd if sd > 1e-9 else 1e-9)
        zv_list.append(zs["var"])
        za_list.append(zs["ac"])
        per_class[cls] = {"var": var, "ac": ac, "z_var": zs["var"], "z_ac": zs["ac"]}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        z_var = np.nanmean(np.vstack(zv_list), axis=0)
        z_ac = np.nanmean(np.vstack(za_list), axis=0)
    return {"z_var": z_var, "z_ac": z_ac, "per_class": per_class}


def alarm_from_z(z_var: np.ndarray, z_ac: np.ndarray, k: float, sustain: int,
                 combine: str = "and", start: int = 0) -> dict[str, Any]:
    """Execute the level rule.  combine: "and" (both legs), "mean" (composite
    index (z_var+z_ac)/2), "var_only" (variance leg alone)."""
    with np.errstate(invalid="ignore"):
        if combine == "var_only":
            score = z_var
        elif combine == "mean":
            score = 0.5 * (z_var + z_ac)
        else:
            score = np.minimum(z_var, z_ac)
    cond = np.asarray(score >= k, bool).copy()
    cond[np.isnan(score)] = False
    if start:
        cond[:start] = False
    alarm = None
    run = 0
    for i, ok in enumerate(cond):
        run = run + 1 if ok else 0
        if run >= sustain:
            alarm = i
            break
    valid = score[start:]
    valid = valid[np.isfinite(valid)]
    return {"alarm_cycle": alarm, "score": score,
            "max_score": float(valid.max()) if valid.size else float("nan")}


def alarm_from_trends(tau_var: np.ndarray, tau_ac: np.ndarray, thr: float,
                      sustain: int, combine: str = "and",
                      burn_in: int = 0) -> dict[str, Any]:
    with np.errstate(invalid="ignore"):
        if combine == "var_only":
            cond = tau_var >= thr
            score = np.where(np.isnan(tau_var), -np.inf, tau_var)
        else:
            cond = (tau_var >= thr) & (tau_ac >= thr)
            score = np.where(np.isnan(tau_var) | np.isnan(tau_ac), -np.inf,
                             np.minimum(tau_var, tau_ac))
    cond = np.asarray(cond, bool).copy()
    if burn_in:
        cond[:burn_in] = False
    alarm = None
    run = 0
    for i, ok in enumerate(cond):
        run = run + 1 if ok else 0
        if run >= sustain:
            alarm = i
            break
    valid = score[burn_in:]
    valid = valid[np.isfinite(valid)]
    return {"alarm_cycle": alarm, "score": score,
            "max_score": float(valid.max()) if valid.size else float("nan")}


def leading_alarm(series: dict[str, np.ndarray], rule: dict[str, Any]) -> dict[str, Any]:
    """Execute the pre-declared leading-indicator rule on the OBSERVED series.

    rule keys (all pre-declared in configs/detection_rule.json):
      window        rolling-window length for var / AR(1)
      tau_span      number of trailing window-statistics entered into Kendall tau
      tau_threshold tau* the mixed trend statistic must reach
      sustain       consecutive cycles the condition must hold
      combine       "and" (both variance and AR(1) trends) | "var_only"
      classes       signal classes mixed at the trend-statistic level
      min_frac      minimum fraction of non-missing cycles in a window
      max_gap       maximum interpolated gap
      burn_in       no alarm may be raised before this cycle
    """
    tr = class_trends(series, int(rule["window"]), int(rule["tau_span"]),
                      list(rule["classes"]), float(rule.get("min_frac", 0.5)),
                      int(rule.get("max_gap", 3)))
    res = alarm_from_trends(tr["tau_var"], tr["tau_ac"],
                            float(rule["tau_threshold"]), int(rule["sustain"]),
                            rule.get("combine", "and"),
                            int(rule.get("burn_in", 0)))
    res["tau_var"] = tr["tau_var"]
    res["tau_ac"] = tr["tau_ac"]
    res["per_class"] = tr["per_class"]
    return res


def lagging_alarm(output: np.ndarray, rule: dict[str, Any]) -> dict[str, Any]:
    """Pre-declared lagging rule on the per-program output metric recorded at
    program close: alarm when the cycle-mean output falls more than `k_sd`
    baseline standard deviations below the baseline mean, sustained."""
    base_end = int(rule["baseline_end"])
    k_sd = float(rule["k_sd"])
    sustain = int(rule["sustain"])
    y = fill_gaps(np.asarray(output, float), int(rule.get("max_gap", 3)))
    base = y[:base_end]
    base = base[~np.isnan(base)]
    if base.size < 10:
        return {"alarm_cycle": None, "threshold": float("nan")}
    mu, sd = float(base.mean()), float(base.std(ddof=1))
    thr = mu - k_sd * sd
    alarm = None
    run = 0
    for i in range(base_end, y.size):
        v = y[i]
        run = run + 1 if (not np.isnan(v) and v < thr) else 0
        if run >= sustain:
            alarm = i
            break
    return {"alarm_cycle": alarm, "threshold": thr, "baseline_mean": mu,
            "baseline_sd": sd}


# --------------------------------------------------------------------------
# discrimination statistic (AUC via Mann-Whitney U), no scipy
# --------------------------------------------------------------------------
def auc(pos: Sequence[float], neg: Sequence[float]) -> float:
    p = np.asarray([v for v in pos if np.isfinite(v)], float)
    q = np.asarray([v for v in neg if np.isfinite(v)], float)
    if p.size == 0 or q.size == 0:
        return float("nan")
    allv = np.concatenate([p, q])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(allv.size, float)
    sorted_v = allv[order]
    i = 0
    while i < sorted_v.size:
        j = i
        while j + 1 < sorted_v.size and sorted_v[j + 1] == sorted_v[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    r_pos = ranks[:p.size].sum()
    u = r_pos - p.size * (p.size + 1) / 2.0
    return float(u / (p.size * q.size))


# --------------------------------------------------------------------------
# logistic propensity fit (Ensign-fix arm), no scipy
# --------------------------------------------------------------------------
def logistic_fit(X: np.ndarray, y: np.ndarray, iters: int = 200,
                 lam: float = 1e-3) -> np.ndarray:
    """Newton-IRLS logistic regression with a small ridge penalty.
    X already includes an intercept column."""
    n, d = X.shape
    w = np.zeros(d)
    for _ in range(iters):
        eta = X @ w
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        g = X.T @ (y - p) - lam * w
        W = p * (1.0 - p) + 1e-9
        H = (X * W[:, None]).T @ X + lam * np.eye(d)
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            break
        w = w + step
        if np.max(np.abs(step)) < 1e-8:
            break
    return w


def logistic_predict(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(X @ w, -30, 30)))
