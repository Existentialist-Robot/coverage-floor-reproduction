"""
Figure production for the three frozen figures (described below), plus a
Gate-1 appendix figure.

  fig1  E-A headline, Ferretti Fig. 3 form: the (c x eps) plane with Delta t
        heat-mapped, a SOLID Delta t = 0 iso-line (the coverage floor), a DASHED
        seed-variation band on that iso-line, and PANELS ALONG NOISE COLOUR
        (white / red) -- the panel allocation is fixed across the figure set.
  fig2  Confusion as a TIME SERIES with a cost axis (Cencetti Fig. 6a/6b form):
        false positives and false negatives against cycle, per coverage level.
  fig3  E-B: the (c x g) grid with corr(observed indicator, latent A(t))
        heat-mapped and the coupling-severed control as a reference contour.
  figA1 Gate 1: branch diagram with hysteresis, and the divergence of the
        post-shock recovery time.

Follow-up campaign 2026-08-05 adds three figures and re-emits fig1-fig3 from the
recalibrated aggregates (suffix `_recal` / `_n500`, the campaign-1 versions are kept):
  fig4  the Delta t = 0 iso-line firmed up: Delta t against c per noise colour with
        +-1 s.e. bands at ~100 seeds/cell on the iso-line neighbourhood.
  fig5  the calibration null: FPR against threshold on 600 fresh control runs, with
        the declared target, the external bar, and where campaign 1's threshold sat.
  fig6  E-B coupling-loss forest: the paired (severed - coupled) tracking difference
        with 95% CIs per (g, c) and pooled per g at >= 500 matched pairs.

Uses matplotlib when it is importable, and otherwise falls back to the
dependency-free SVG writer in svgwrite_min.py.  Both paths emit SVG; the matplotlib
path also emits PNG.

Run:  python figures.py
Writes: figures/*.png, figures/*.svg
"""

from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AGG = os.path.join(HERE, "aggregates")
FIG = os.path.join(HERE, "figures")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from matplotlib.colors import TwoSlopeNorm
    HAVE_MPL = True
except Exception:                                     # pragma: no cover
    HAVE_MPL = False
    import svgwrite_min as svgw

# --- print sizing (2026-08-22) ---------------------------------------------
# Every figure is now DRAWN at the width it occupies on the printed page
# (Springer svproc \textwidth = 122 mm = 4.803 in) instead of being drawn
# 1.5-2.3x oversize and shrunk by LaTeX -- the shrink was landing text authored
# at 7.4-10.5 pt at an effective 1.4-2 pt in the printed paper (fig4 at
# 0.72\linewidth was 18% of authored size).  Fonts are chosen so nothing lands
# below ~7 pt at final size.  This is typography and geometry only: the data,
# panel allocation and every plotted number are unchanged.
PRINT_W = 4.803                       # inches; = 122 mm, the svproc \textwidth
FS_TITLE, FS_LABEL, FS_TICK, FS_ANN = 8.0, 7.5, 7.0, 7.0
if HAVE_MPL:
    plt.rcParams.update({
        "font.size": 7.5,
        "axes.titlesize": FS_TITLE,
        "axes.labelsize": FS_LABEL,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_ANN,
        "figure.titlesize": FS_TITLE,
    })

SEVERITIES = ("low", "mid", "high")
PANEL_COLOURS = ("white", "red")


def load(name):
    with open(os.path.join(AGG, name)) as f:
        return json.load(f)


def _save(fig, stem):
    os.makedirs(FIG, exist_ok=True)
    # dpi 300 (was 200): the figures now print at 1:1, so the raster IS the
    # print resolution.  pad_inches 0.02 (was the 0.1 default): the pad is dead
    # width, and a figure authored at \textwidth that ships wider than
    # \textwidth gets scaled back down by LaTeX, shrinking its fonts below the
    # size they were set at.
    fig.savefig(os.path.join(FIG, stem + ".png"), dpi=300, bbox_inches="tight",
                pad_inches=0.02)
    fig.savefig(os.path.join(FIG, stem + ".svg"), bbox_inches="tight",
                pad_inches=0.02)
    plt.close(fig)
    print(f"  wrote figures/{stem}.png + .svg")


# ---------------------------------------------------------------- figure 1
def fig1_headline(arm: str = "mnar", src: str = "ea_headline.json",
                  suffix: str = "") -> None:
    d = load(src)
    cs = sorted({c["c"] for c in d["cells"]})
    grids, isos = {}, {}
    for col in PANEL_COLOURS:
        g = np.full((len(SEVERITIES), len(cs)), np.nan)
        se = np.full_like(g, np.nan)
        for x in d["cells"]:
            if x["arm"] != arm or x["colour"] != col:
                continue
            i = SEVERITIES.index(x["severity"])
            j = cs.index(x["c"])
            g[i, j] = x["dt_censored_mean"]
            se[i, j] = x["dt_censored_se"]
        grids[col] = g
        isos[col] = se
    vals = np.concatenate([grids[c].ravel() for c in PANEL_COLOURS])
    vmax = float(np.nanmax(np.abs(vals))) or 1.0

    if not HAVE_MPL:                                   # pragma: no cover
        panels = []
        for col in PANEL_COLOURS:
            panels.append({"title": f"{col} noise",
                           "grid": [[None if np.isnan(v) else float(v)
                                     for v in row] for row in grids[col]],
                           "isoline": [], "band": []})
        svgw.heatmap_panels(
            os.path.join(FIG, "fig1_ea_headline.svg"), panels,
            [f"c={c:g}" for c in cs], list(SEVERITIES),
            "consent coverage c", "proxy-error severity",
            "lead time dt (cycles)",
            f"E-A: lead-time advantage, {arm.upper()} consent", -vmax, vmax)
        print("  wrote figures/fig1_ea_headline.svg (SVG fallback)")
        return

    fig, axes = plt.subplots(1, len(PANEL_COLOURS), figsize=(PRINT_W, 2.4),
                             sharey=True)
    fig.subplots_adjust(top=0.70)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    for ax, col in zip(np.atleast_1d(axes), PANEL_COLOURS):
        g = grids[col]
        im = ax.pcolormesh(np.arange(len(cs) + 1), np.arange(len(SEVERITIES) + 1),
                           g, cmap="RdBu_r", norm=norm, edgecolors="w", lw=0.6)
        for i in range(g.shape[0]):
            for j in range(g.shape[1]):
                if np.isfinite(g[i, j]):
                    ax.text(j + 0.5, i + 0.5, f"{g[i, j]:+.0f}", ha="center",
                            va="center", fontsize=FS_ANN, color="black")
        # Delta t = 0 iso-line, plus dashed band at +-1 seed-variation se
        xc = np.arange(len(cs)) + 0.5
        yc = np.arange(len(SEVERITIES)) + 0.5
        X, Y = np.meshgrid(xc, yc)
        if np.isfinite(g).sum() >= 4 and np.nanmin(g) < 0 < np.nanmax(g):
            ax.contour(X, Y, g, levels=[0.0], colors="k", linewidths=1.5)
            for sgn, ls in ((-1, "--"), (+1, "--")):
                ax.contour(X, Y, g + sgn * isos[col], levels=[0.0], colors="k",
                           linewidths=0.8, linestyles=ls)
        ax.set_xticks(xc)
        ax.set_xticklabels([f"{c:g}" for c in cs])
        ax.set_yticks(yc)
        ax.set_yticklabels(list(SEVERITIES))
        ax.set_xlabel("consent coverage $c$")
        ax.set_title(f"{col} noise")
    np.atleast_1d(axes)[0].set_ylabel("proxy-error severity $\\epsilon$")
    cb = fig.colorbar(im, ax=np.atleast_1d(axes).tolist(), pad=0.02)
    cb.set_label("lead time $\\Delta t$ (cycles)")
    extra = (d.get("design") or {}).get("c_extra_points") or []
    note = ("\nblank column at $c$=" + ", ".join(f"{x:g}" for x in extra)
            + ": added iso-line point, headline severity only (amendment A-9)"
            if extra else "")
    fig.suptitle(f"E-A: lead-time advantage over the lagging output metric "
                 f"({arm.upper()} consent)\n"
                 "solid = $\\Delta t = 0$ iso-line (the coverage floor); "
                 "dashed = $\\pm1$ s.e. seed-variation band" + note,
                 fontsize=FS_TITLE)
    _save(fig, f"fig1_ea_headline_{arm}{suffix}")


# ---------------------------------------------------------------- figure 2
def fig2_confusion(arm: str = "mnar", colour: str = "red",
                   src: str = "ea_headline.json", suffix: str = "") -> None:
    d = load(src)
    n_cycles = d["design"]["n_cycles"]
    keys = [k for k in d["confusion_series"]
            if k.startswith(f"{arm}|{colour}|")]
    keys.sort(key=lambda k: float(k.split("|")[-1]))
    t = np.arange(n_cycles)
    series = []
    for k in keys:
        v = d["confusion_series"][k]
        c = float(k.split("|")[-1])
        tr = [x for x in v["alarm_cycles_transition"]]
        ct = [x for x in v["alarm_cycles_control"]]
        fn = np.array([np.mean([1.0 if (x is None or x > ti) else 0.0 for x in tr])
                       for ti in t])
        fp = np.array([np.mean([1.0 if (x is not None and x <= ti) else 0.0
                                for x in ct])
                       for ti in t])
        series.append((c, fp, fn))

    if not HAVE_MPL:                                   # pragma: no cover
        ss = []
        pal = ["#440154", "#31688e", "#35b779", "#fde725"]
        for i, (c, fp, fn) in enumerate(series):
            ss.append({"x": t.tolist(), "y": fn.tolist(),
                       "label": f"FN c={c:g}", "colour": pal[i % 4]})
            ss.append({"x": t.tolist(), "y": fp.tolist(),
                       "label": f"FP c={c:g}", "colour": pal[i % 4], "dash": "5,4"})
        svgw.lines(os.path.join(FIG, "fig2_confusion.svg"), ss, "cycle",
                   "rate", f"E-A confusion over time ({arm.upper()}, {colour} noise)",
                   ylim=(0, 1))
        print("  wrote figures/fig2_confusion.svg (SVG fallback)")
        return

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(PRINT_W, 3.7), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1]})
    cmap = plt.get_cmap("viridis")
    for i, (c, fp, fn) in enumerate(series):
        col = cmap(i / max(1, len(series) - 1))
        ax.plot(t, fn, color=col, lw=1.3, label=f"false negative, $c={c:g}$")
        ax.plot(t, fp, color=col, lw=0.9, ls="--",
                label=f"false positive, $c={c:g}$")
    ax.axhline(0.036, color="crimson", lw=0.9, ls=":")
    ax.text(4, 0.055, "external FPR bar 3.6% (Falmagne et al. 2026)",
            fontsize=FS_ANN, color="crimson")
    ax.set_ylabel("rate")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=FS_ANN, ncol=2, loc="center left")
    ax.set_title(f"E-A confusion as a time series ({arm.upper()} consent, "
                 f"{colour} noise, mid severity)", fontsize=FS_TITLE)
    # cost axis: unit cost per false alarm, unit cost per cycle of missed warning
    for i, (c, fp, fn) in enumerate(series):
        col = cmap(i / max(1, len(series) - 1))
        ax2.plot(t, np.cumsum(fn) / len(t) + fp, color=col, lw=1.2,
                 label=f"$c={c:g}$")
    ax2.set_xlabel("cycle")
    ax2.set_ylabel("cost\n(1/alarm + 1/cycle missed)", fontsize=FS_LABEL)
    ax2.legend(fontsize=FS_ANN, ncol=4)
    _save(fig, f"fig2_confusion_{arm}_{colour}{suffix}")


# ---------------------------------------------------------------- figure 3
def fig3_eb(src: str = "eb_joint_sweep.json", stem: str = "fig3_eb_joint") -> None:
    d = load(src)
    cs = sorted({x["c"] for x in d["cells"]})
    gs = ["0", "mid", "high"]
    cp = np.full((len(gs), len(cs)), np.nan)
    sv = np.full_like(cp, np.nan)
    for x in d["cells"]:
        i, j = gs.index(x["g_name"]), cs.index(x["c"])
        (cp if x["mode"] == "coupled" else sv)[i, j] = x["corr_mean"]

    if not HAVE_MPL:                                   # pragma: no cover
        svgw.heatmap_panels(
            os.path.join(FIG, "fig3_eb_joint.svg"),
            [{"title": "coupled", "grid": [[None if np.isnan(v) else float(v * 100)
                                            for v in r] for r in cp]},
             {"title": "coupling severed (control)",
              "grid": [[None if np.isnan(v) else float(v * 100) for v in r]
                       for r in sv]}],
            [f"c={c:g}" for c in cs], gs, "consent coverage c", "reward coupling g",
            "corr x100", "E-B: corr(observed indicator, latent A(t))", -100, 0)
        print("  wrote figures/fig3_eb_joint.svg (SVG fallback)")
        return

    fig, ax = plt.subplots(figsize=(PRINT_W, 3.1))
    im = ax.pcolormesh(np.arange(len(cs) + 1), np.arange(len(gs) + 1), cp,
                       cmap="magma", vmin=float(np.nanmin([cp.min(), sv.min()])),
                       vmax=0.0, edgecolors="w", lw=0.6)
    for i in range(cp.shape[0]):
        for j in range(cp.shape[1]):
            ax.text(j + 0.5, i + 0.62, f"{cp[i, j]:+.2f}", ha="center",
                    va="center", fontsize=FS_LABEL, color="white")
            ax.text(j + 0.5, i + 0.30, f"(severed {sv[i, j]:+.2f})", ha="center",
                    va="center", fontsize=FS_ANN, color="0.85")
    xc, yc = np.arange(len(cs)) + 0.5, np.arange(len(gs)) + 0.5
    X, Y = np.meshgrid(xc, yc)
    levels = np.round(np.linspace(np.nanmin(sv), np.nanmax(sv), 4), 2)
    cn = ax.contour(X, Y, sv, levels=sorted(set(levels)), colors="cyan",
                    linewidths=0.9, linestyles="--")
    ax.clabel(cn, fontsize=FS_ANN, fmt="%.2f")
    ax.set_xticks(xc)
    ax.set_xticklabels([f"{c:g}" for c in cs])
    ax.set_yticks(yc)
    ax.set_yticklabels(gs)
    ax.set_xlabel("consent coverage $c$")
    ax.set_ylabel("reward coupling $g$")
    ax.set_title("E-B: corr(observed indicator, latent $A(t)$)\n"
                 "cell = coupled; dashed cyan = coupling-severed control contour "
                 "(more negative = better tracking)", fontsize=FS_TITLE)
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("corr(indicator, $A(t)$)")
    _save(fig, stem)


# ------------------------------------------------- figure 4 (follow-up campaign)
def fig4_isoline(src: str = "ea_headline_recal.json",
                 stem: str = "fig4_isoline_recal") -> None:
    """Delta t against coverage per noise colour, with the +-1 s.e. band that the
    iso-line's uncertainty is read off.  One panel per consent arm."""
    d = load(src)
    if not HAVE_MPL:                                       # pragma: no cover
        ss = []
        pal = {"white": "#31688e", "red": "#b2182b", "systematic": "#35b779"}
        for key, e in d["iso_line_dt_zero"].items():
            if not key.startswith("mnar"):
                continue
            ss.append({"x": e["c"], "y": e["dt"], "label": key,
                       "colour": pal.get(key.split("|")[1], "#333")})
        svgw.lines(os.path.join(FIG, stem + ".svg"), ss, "consent coverage c",
                   "lead time dt (cycles)",
                   "E-A: lead time against coverage (Delta-RMTA = 0 crossing)")
        print(f"  wrote figures/{stem}.svg (SVG fallback)")
        return
    arms = ("mnar", "mcar")
    # Drawn AT print size (see the print-sizing block at the top of this file):
    # this figure is Fig. 2 of the paper, included at \linewidth = 122 mm, so
    # 4.803 x 1.76 in here is 1:1 with the page and every label is a true >= 7 pt.
    # The suptitle ("E-A: the Delta t = 0 coverage floor, firmed up" -- a name
    # the 2026-08-15 number audit records as DEPRECATED; the caption calls the
    # same line the Delta-RMTA = 0 crossing) and the subtitle line are gone from
    # the artwork: the LaTeX caption already carries all of it, and dropping
    # them buys back page height.
    fig, axes = plt.subplots(1, 2, figsize=(PRINT_W, 1.76), sharey=True)
    # top lowered 0.82 -> 0.76 to fix a legend/title overlap: the shared legend at bbox y=1.02
    # was crowding the panel titles. Taken out of the plot area rather than out of the page
    # budget -- the figure height is unchanged, so pagination cannot move.
    fig.subplots_adjust(left=0.09, right=0.99, top=0.76, bottom=0.22, wspace=0.08)
    pal = {"white": "#2166ac", "red": "#b2182b", "systematic": "#1b7837"}
    # Greyscale differentiation (print-polish-report.md item 11): the three
    # noise-colour series had identical `o-` markers and near-identical greyscale
    # luminance (0.35 / 0.28 / 0.33), so they collapsed to one grey in print.
    # Marker shape AND linestyle now carry the series identity; colour is
    # redundant reinforcement rather than the sole encoding.
    style = {"white": ("o", "-"), "red": ("^", "--"), "systematic": ("s", ":")}
    hatch = {"white": None, "red": "//", "systematic": ".."}  # sparse hatching
    # stays legible now that the figure prints at 1:1 rather than ~0.3x.
    for ax, arm in zip(axes, arms):
        for iw, col in enumerate(("white", "red", "systematic")):
            e = d["iso_line_dt_zero"].get(f"{arm}|{col}")
            if not e:
                continue
            c = np.array(e["c"], float)
            y = np.array(e["dt"], float)
            se = np.array(e["dt_se"], float)
            mk, ls = style[col]
            ax.plot(c, y, marker=mk, linestyle=ls, color=pal[col], lw=1.1, ms=2.8,
                    label=col)
            ax.fill_between(c, y - se, y + se, color=pal[col], alpha=0.16, lw=0,
                            hatch=hatch[col], edgecolor=pal[col])
            x0 = e.get("c_at_dt_zero")
            if x0 is not None:
                ax.axvline(x0, color=pal[col], lw=0.8, ls=":")
                ax.annotate(f"$c_0$={x0:.3f}", (x0, 0), textcoords="offset points",
                            xytext=(3, 10 + 12 * iw), fontsize=FS_ANN, color=pal[col],
                            bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.4))
        ax.axhline(0.0, color="k", lw=0.8)
        ax.set_xlabel("consent coverage $c$")
        ax.set_title(f"{arm.upper()} consent", fontsize=FS_TITLE)
    axes[0].set_ylabel("lead time $\\Delta t$ (cycles)")
    # Explicit tick count, added 2026-08-22: lowering the axes top to open the legend gap
    # shortened the plot area enough that the auto-locator dropped to two labelled ticks, on a
    # panel whose whole subject is where the curve crosses zero. Four keeps the scale readable
    # at the same height.
    for _ax in axes:
        _ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
    # One shared legend for both panels (the series are identical); the per-cell
    # replication ("30-100 runs per cell") is stated in the LaTeX caption, which
    # is why the old per-series "(n=30-100/cell)" suffix -- identical for all
    # six (arm, colour) series in the aggregate -- is not repeated here.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=FS_ANN,
               frameon=False, bbox_to_anchor=(0.5, 1.02), columnspacing=1.4,
               handlelength=1.8)
    _save(fig, stem)


# ------------------------------------------------- figure 5 (follow-up campaign)
def fig5_calibration(src: str = "calibration_null_large.json",
                     stem: str = "fig5_calibration_null") -> None:
    d = load(src)
    if not HAVE_MPL:                                       # pragma: no cover
        cur = d["leading_fpr_curve"]
        svgw.lines(os.path.join(FIG, stem + ".svg"),
                   [{"x": [x["k"] for x in cur], "y": [x["fpr"] for x in cur],
                     "label": "null FPR", "colour": "#2166ac"}],
                   "threshold k", "false-alarm rate",
                   "threshold recalibration on a large null")
        print(f"  wrote figures/{stem}.svg (SVG fallback)")
        return
    cur = d["leading_fpr_curve"]
    k = np.array([x["k"] for x in cur], float)
    f = np.array([x["fpr"] for x in cur], float)
    n = d["design"]["n_null_control_runs"]
    fig, ax = plt.subplots(figsize=(PRINT_W, 2.9))
    ax.plot(k, f, "-", color="#2166ac", lw=1.3,
            label=f"leading-rule FPR on {n} fresh null runs")
    ax.axhline(d["design"]["previous_nominal_fpr"], color="0.5", lw=0.8, ls="--")
    ax.axhline(d["design"]["target_fpr"], color="#1b7837", lw=0.9, ls="-.")
    ax.axhline(d["design"]["external_bar"], color="crimson", lw=0.9, ls=":")
    k_old = d["campaign1_threshold_on_this_null"]["k"]
    f_old = d["campaign1_threshold_on_this_null"]["fpr"]
    k_new = d["chosen"]["leading_k"]
    f_new = d["chosen"]["leading_calibration_fpr"]
    ax.plot([k_old], [f_old], "o", ms=5, color="crimson")
    ax.annotate(f"campaign 1: $k$={k_old:g}\ntrue null FPR {f_old:.1%}\n"
                f"(picked on 40 runs at a nominal 5%)",
                (k_old, f_old), textcoords="offset points", xytext=(8, 5),
                fontsize=FS_ANN, color="crimson")
    ax.plot([k_new], [f_new], "s", ms=5, color="#1b7837")
    ax.annotate(f"recalibrated: $k$={k_new:g}\nnull FPR {f_new:.1%}",
                (k_new, f_new), textcoords="offset points", xytext=(9, 13),
                fontsize=FS_ANN, color="#1b7837")
    kmax = float(k[f > 0].max()) + 0.4 if np.any(f > 0) else float(k[-1])
    ax.set_xlim(float(k[0]) - 0.2, kmax)
    ax.text(kmax, d["design"]["external_bar"], "external bar 3.6% ", va="bottom",
            ha="right", fontsize=FS_ANN, color="crimson")
    ax.text(kmax, d["design"]["previous_nominal_fpr"], "campaign-1 nominal 5% ",
            va="bottom", ha="right", fontsize=FS_ANN, color="0.4")
    ax.text(kmax, d["design"]["target_fpr"], "declared target 3.5% ", va="top",
            ha="right", fontsize=FS_ANN, color="#1b7837")
    ax.set_yscale("log")
    ax.set_xlabel("composite z-score threshold $k$ (sustained 5 cycles)")
    ax.set_ylabel("false-alarm rate on the calibration null")
    ax.set_title("Threshold recalibration: the 40-run threshold was overfitted\n"
                 "rule form unchanged; only $k$ re-picked, target 3.5%",
                 fontsize=FS_TITLE)
    ax.legend(fontsize=FS_ANN, loc="lower left")
    _save(fig, stem)


# ------------------------------------------------- figure 6 (follow-up campaign)
def fig6_eb_forest(src: str = "eb_joint_sweep_n500.json",
                   stem: str = "fig6_eb_coupling_forest") -> None:
    """The E-B verdict as a forest plot: paired (severed - coupled) tracking loss with
    95% CIs, per (g, c) and pooled per g.  Zero line = no coupling effect."""
    d = load(src)
    dec = d["decision"]
    labels, mean, lo, hi, kind = [], [], [], [], []
    for g in ("mid", "high"):
        v = dec["per_g"].get(g)
        if not v:
            continue
        for cs, pc in sorted(v["per_c_tracking_loss"].items(), key=lambda t: float(t[0])):
            labels.append(f"g={g}, c={float(cs):g}  (n={pc['n_pairs']})")
            mean.append(pc["tracking_loss_mean"])
            lo.append(pc["lo95"])
            hi.append(pc["hi95"])
            kind.append("cell")
        labels.append(f"g={g}, pooled over c  (n={v['n_pairs_pooled']})")
        mean.append(v["tracking_loss_mean"])
        lo.append(v["lo95"])
        hi.append(v["hi95"])
        kind.append("pooled")
    if not HAVE_MPL:                                       # pragma: no cover
        svgw.lines(os.path.join(FIG, stem + ".svg"),
                   [{"x": mean, "y": list(range(len(mean))), "label": "loss",
                     "colour": "#2166ac"}], "tracking loss", "cell",
                   "E-B coupling-loss CIs")
        print(f"  wrote figures/{stem}.svg (SVG fallback)")
        return
    y = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(PRINT_W, 0.26 * len(labels) + 1.2))
    for yi, m, a, b, kd in zip(y, mean, lo, hi, kind):
        col = "#111111" if kd == "pooled" else "#2166ac"
        ax.plot([a, b], [yi, yi], "-", color=col, lw=1.8 if kd == "pooled" else 1.1)
        ax.plot([m], [yi], "D" if kd == "pooled" else "o", color=col,
                ms=4.5 if kd == "pooled" else 3.2)
    ax.axvline(0.0, color="crimson", lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=FS_ANN)
    ax.set_xlabel("paired tracking loss from coupling  "
                  "$|corr|_{severed} - |corr|_{coupled}$   (>0 = coupling costs "
                  "tracking)")
    ax.set_title(f"E-B at {d['design']['reps']} matched pairs/cell: "
                 f"{dec['verdict']}\nbootstrap 95% CIs; diamonds = pooled over "
                 f"coverage at fixed $g$", fontsize=FS_TITLE)
    _save(fig, stem)


# ---------------------------------------------------------------- appendix
def figA1_gate1() -> None:
    d = load("gate1_bifurcation.json")
    bt = d["branch_table"]
    rel = d["relaxation_table"]
    ns = sorted({r["n"] for r in bt})
    if not HAVE_MPL:                                   # pragma: no cover
        ss = []
        pal = ["#440154", "#31688e", "#35b779", "#fde725"]
        for i, n in enumerate(ns):
            for br, dash in (("upper", None), ("lower", "5,4")):
                rs = sorted([r for r in bt if r["n"] == n and r["branch"] == br],
                            key=lambda r: r["a"])
                ss.append({"x": [r["a"] * 100 for r in rs],
                           "y": [r["z_mean"] for r in rs],
                           "label": f"N={n} {br}", "colour": pal[i % 4],
                           "dash": dash})
        svgw.lines(os.path.join(FIG, "figA1_gate1.svg"), ss, "A x100",
                   "mean repeat-tie degree z", "Gate 1: hysteresis")
        print("  wrote figures/figA1_gate1.svg (SVG fallback)")
        return
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(PRINT_W, 2.1))
    fig.subplots_adjust(wspace=0.42, left=0.10, right=0.99, top=0.80, bottom=0.20)
    cmap = plt.get_cmap("viridis")
    for i, n in enumerate(ns):
        col = cmap(i / max(1, len(ns) - 1))
        for br, ls, mk in (("upper", "-", "o"), ("lower", "--", "s")):
            rs = sorted([r for r in bt if r["n"] == n and r["branch"] == br],
                        key=lambda r: r["a"])
            ax.plot([r["a"] for r in rs], [r["z_mean"] for r in rs], ls,
                    marker=mk, ms=2.2, color=col, lw=1.0,
                    label=f"$N={n}$ {br}" if br == "upper" else None)
    ax.axhline(1.0, color="0.6", lw=0.7, ls=":")
    ax.set_xlabel("latent adaptive capacity $A$")
    ax.set_ylabel("mean repeat-tie degree $z$")
    ax.set_title("branch diagram\nsolid = from above, dashed = from below\n"
                 "(bistable window = hysteresis)", fontsize=FS_ANN)
    ax.legend(fontsize=FS_ANN)
    a = [r["a"] for r in rel]
    rt = [r["recovery_time_cycles"] for r in rel]
    ax2.plot(a, rt, "o-", color="crimson", lw=1.2, ms=2.5)
    ax2.set_yscale("log")
    ax2.invert_xaxis()
    ax2.set_xlabel("latent adaptive capacity $A$ (approaching the fold)")
    ax2.set_ylabel("post-shock recovery time (cycles)")
    ax2.set_title("critical slowing down measured\non the latent structure",
                  fontsize=FS_ANN)
    _save(fig, "figA1_gate1")


def main() -> None:
    os.makedirs(FIG, exist_ok=True)
    print(f"figures ({'matplotlib' if HAVE_MPL else 'pure-SVG fallback'}):")
    have = set(os.listdir(AGG))
    if "ea_headline.json" in have:
        for arm in ("mnar", "mcar"):
            fig1_headline(arm)
        fig2_confusion("mnar", "red")
        fig2_confusion("mnar", "white")
    if "eb_joint_sweep.json" in have:
        fig3_eb()
    if "gate1_bifurcation.json" in have:
        figA1_gate1()
    # ---- follow-up campaign 2026-08-05 -------------------------------------
    if "ea_headline_recal.json" in have:
        for arm in ("mnar", "mcar"):
            fig1_headline(arm, "ea_headline_recal.json", "_recal")
        fig2_confusion("mnar", "red", "ea_headline_recal.json", "_recal")
        fig2_confusion("mnar", "white", "ea_headline_recal.json", "_recal")
        fig4_isoline()
    if "calibration_null_large.json" in have:
        fig5_calibration()
    if "eb_joint_sweep_n500.json" in have:
        fig3_eb("eb_joint_sweep_n500.json", "fig3_eb_joint_n500")
        fig6_eb_forest()


if __name__ == "__main__":
    main()
