"""Chi-square upper tail for ANY degrees of freedom, stdlib only.

`remediation_analysis.chi2_sf_even` is a closed form that is correct and exact for
even df and returns NaN for odd df.  Every coverage-heterogeneity test in this study
runs on len(C_GRID) - 1 degrees of freedom, so a four-cell grid gives df = 3 and the
p-value came back NaN: the 2026-09-03 campaign reported a chi-square statistic with no
tail probability attached to it.  This module supplies the general case so a
heterogeneity test on any grid size reports a p-value.

Q(a, x) is the regularised upper incomplete gamma function.  Series below x = a + 1,
Lentz continued fraction above it -- the standard split, because each form converges
quickly only on its own side.  chi2_sf(x, df) = Q(df/2, x/2).
"""

from __future__ import annotations

import math

_EPS = 3.0e-16
_TINY = 1.0e-300
_MAXIT = 1000


def _gamma_p_series(a: float, x: float) -> float:
    """Lower regularised P(a, x) by its power series.  Converges fast for x < a + 1."""
    term = 1.0 / a
    total = term
    n = a
    for _ in range(_MAXIT):
        n += 1.0
        term *= x / n
        total += term
        if abs(term) < abs(total) * _EPS:
            break
    else:
        raise RuntimeError(f"gamma series did not converge at a={a}, x={x}")
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gamma_q_cf(a: float, x: float) -> float:
    """Upper regularised Q(a, x) by modified Lentz continued fraction, x >= a + 1."""
    b = x + 1.0 - a
    c = 1.0 / _TINY
    d = 1.0 / b
    h = d
    for i in range(1, _MAXIT + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < _TINY:
            d = _TINY
        c = b + an / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _EPS:
            break
    else:
        raise RuntimeError(f"gamma continued fraction did not converge at a={a}, x={x}")
    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


def gamma_q(a: float, x: float) -> float:
    """Regularised upper incomplete gamma Q(a, x) = 1 - P(a, x)."""
    if not (a > 0.0):
        raise ValueError(f"a must be positive, got {a}")
    if x < 0.0:
        raise ValueError(f"x must be non-negative, got {x}")
    if x == 0.0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - _gamma_p_series(a, x)
    return _gamma_q_cf(a, x)


def chi2_sf(x: float, df: int) -> float:
    """Upper tail P(X > x) for a chi-square variate on `df` degrees of freedom.

    Exact for even and odd df alike.  Returns NaN for a non-finite statistic rather
    than raising, because a heterogeneity test on a degenerate pooled rate legitimately
    produces one and the caller reports it as undefined.
    """
    if df < 1:
        raise ValueError(f"df must be at least 1, got {df}")
    if not math.isfinite(x):
        return float("nan")
    if x <= 0.0:
        return 1.0
    return gamma_q(df / 2.0, x / 2.0)


def het_chi2(counts: list[tuple[int, int]]) -> dict:
    """Test that one binomial rate generated every (k, n) group.

    Same statistic as `remediation_analysis.het_chi2`; the only difference is that the
    p-value is defined for odd df.  Kept here rather than patched there so no committed
    aggregate changes as a side effect of this module landing.
    """
    tk = sum(k for k, _ in counts)
    tn = sum(n for _, n in counts)
    p = tk / tn if tn else float("nan")
    df = len(counts) - 1
    if not (0.0 < p < 1.0):
        return {"chi2": float("nan"), "df": df, "p_value": float("nan"),
                "pooled_rate": p, "rates": [k / n if n else float("nan")
                                            for k, n in counts]}
    chi = sum((k - n * p) ** 2 / (n * p * (1 - p)) for k, n in counts)
    return {"chi2": chi, "df": df, "p_value": chi2_sf(chi, df),
            "pooled_rate": p, "rates": [k / n for k, n in counts]}
