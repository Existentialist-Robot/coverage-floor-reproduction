"""
Scoring: apply the frozen detection rule to a run and reduce it to scalars.

Shared by the coverage-and-noise campaign headline sweep, the auxiliary arms and coverage-and-coupling campaign, so all three score
runs identically and every reported number comes from one code path.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from indicators import baseline_z, alarm_from_z, lagging_alarm

CENSOR = "n_cycles"      # a detector that never fires is censored at the run length


def score_run(series: dict[str, np.ndarray], rule: dict[str, Any],
              n_cycles: int) -> dict[str, Any]:
    """Leading alarm from the OBSERVED subgraph only; lagging alarm from the
    program-close output metric.  Both rules are the frozen ones -- nothing here
    is fitted per run."""
    L = rule["leading"]
    z = baseline_z(series, int(L["window"]), tuple(L["classes"]),
                   int(L["start"]), float(L["min_frac"]), int(L["max_gap"]))
    lead = alarm_from_z(z["z_var"], z["z_ac"], float(L["k"]), int(L["sustain"]),
                        L["combine"], int(L["start"]))
    lag = lagging_alarm(series["output_close"], rule["lagging"])
    t_lead = lead["alarm_cycle"]
    t_lag = lag["alarm_cycle"]
    eff_lead = n_cycles if t_lead is None else t_lead
    eff_lag = n_cycles if t_lag is None else t_lag
    return {
        "t_lead": t_lead, "t_lag": t_lag,
        "dt_censored": eff_lag - eff_lead,
        "dt_paired": (t_lag - t_lead) if (t_lead is not None
                                         and t_lag is not None) else None,
        "lead_max_score": lead["max_score"],
        "lead_fired": t_lead is not None,
        "lag_fired": t_lag is not None,
    }


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    def rank(v):
        order = v.argsort(kind="mergesort")
        r = np.empty(v.size, float)
        sv = v[order]
        i = 0
        while i < sv.size:
            j = i
            while j + 1 < sv.size and sv[j + 1] == sv[i]:
                j += 1
            r[order[i:j + 1]] = 0.5 * (i + j) + 1.0
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    if rx.std() < 1e-12 or ry.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def _smooth(x: np.ndarray, w: int) -> np.ndarray:
    """Trailing mean over w cycles -- the form in which an indicator is published."""
    out = np.full(x.size, np.nan)
    for i in range(x.size):
        seg = x[max(0, i - w + 1):i + 1]
        seg = seg[np.isfinite(seg)]
        if seg.size >= max(3, w // 3):
            out[i] = seg.mean()
    return out


EB_SMOOTH = 20      # cycles; declared with the coverage-and-coupling campaign statistic, not tuned per cell
# The correlation is evaluated over the PRE-FOLD stretch, from the end of the reference
# period to the cycle at which A(t) crosses the empirically measured fold (A_c = 0.37,
# Bifurcation validation).  Past the fold the indicator falls again -- the system has already
# transitioned -- so a correlation taken across the transition is non-monotone and
# measures nothing.  A_c is experimenter-side knowledge of the scenario design, not
# something the in-model instrument uses; the window is identical in every cell.
EB_FOLD_A = 0.37


def indicator_vs_latent(series: dict[str, np.ndarray], rule: dict[str, Any]
                        ) -> dict[str, float]:
    """coverage-and-coupling campaign statistic: correlation between the OBSERVED leading indicator and the
    latent adaptive capacity A(t), over the scored stretch of the run.

    Sign convention: the indicator RISES as A(t) FALLS, so a well-tracking
    instrument gives a NEGATIVE correlation; corr -> 0 means the indicator has
    stopped tracking the latent state."""
    L = rule["leading"]
    z = baseline_z(series, int(L["window"]), tuple(L["classes"]),
                   int(L["start"]), float(L["min_frac"]), int(L["max_gap"]))
    if L["combine"] == "var_only":
        ind = z["z_var"]
    elif L["combine"] == "mean":
        ind = 0.5 * (z["z_var"] + z["z_ac"])
    else:
        ind = np.minimum(z["z_var"], z["z_ac"])
    a = series["a_t"]
    start = int(L["start"])
    below = np.flatnonzero(a <= EB_FOLD_A)
    stop = int(below[0]) if below.size else a.size
    stop = max(start + 40, stop)
    x, y = ind[start:stop], a[start:stop]
    m = np.isfinite(x) & np.isfinite(y)
    out = {"n_points": int(m.sum())}
    if m.sum() < 20 or np.std(x[m]) < 1e-12 or np.std(y[m]) < 1e-12:
        out["corr_indicator_latent"] = float("nan")
        out["corr_var_leg"] = float("nan")
        return out
    out["corr_indicator_latent_pearson_raw"] = float(np.corrcoef(x[m], y[m])[0, 1])
    # Reported statistic: Spearman rank correlation of the PUBLISHED (trailing-mean)
    # indicator against A(t).  Rank correlation because the indicator-vs-latent
    # relationship is monotone but strongly non-linear (fluctuations diverge only near
    # the fold), and smoothed because a published indicator is an aggregate.  Both
    # choices are declared here and applied identically in every cell.
    xs = _smooth(x, EB_SMOOTH)
    ms = np.isfinite(xs) & np.isfinite(y)
    out["corr_indicator_latent"] = (_spearman(xs[ms], y[ms]) if ms.sum() >= 20
                                    else float("nan"))
    xv = _smooth(z["z_var"][start:stop], EB_SMOOTH)
    mv = np.isfinite(xv) & np.isfinite(y)
    out["corr_var_leg"] = (_spearman(xv[mv], y[mv]) if mv.sum() >= 20
                           else float("nan"))
    return out


def bootstrap_mean(x, n_boot: int = 2000, seed: int = 12345) -> dict[str, float]:
    a = np.asarray([v for v in x if v is not None and np.isfinite(v)], float)
    if a.size == 0:
        return {"mean": float("nan"), "se": float("nan"),
                "lo95": float("nan"), "hi95": float("nan"), "n": 0}
    rng = np.random.Generator(np.random.PCG64(seed))
    idx = rng.integers(0, a.size, size=(n_boot, a.size))
    bs = a[idx].mean(axis=1)
    return {"mean": float(a.mean()), "se": float(bs.std(ddof=1)),
            "lo95": float(np.quantile(bs, 0.025)),
            "hi95": float(np.quantile(bs, 0.975)), "n": int(a.size)}
