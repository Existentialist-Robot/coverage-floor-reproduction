"""
Minimal dependency-free SVG writer -- the fallback path for the figure scripts when
matplotlib is not installed.  Covers exactly what the three frozen figures need:
a heat-map grid with axis labels and a colour bar, poly-lines (for the iso-line and
its uncertainty band), and simple line/step plots.  Nothing here imports anything
outside the standard library.
"""

from __future__ import annotations

import html
from typing import Sequence

W, H = 900, 560
PAD_L, PAD_R, PAD_T, PAD_B = 90, 130, 60, 70


def _viridis(t: float) -> str:
    """A short piecewise-linear viridis approximation (no matplotlib)."""
    stops = [(0.0, (68, 1, 84)), (0.25, (59, 82, 139)), (0.5, (33, 145, 140)),
             (0.75, (94, 201, 98)), (1.0, (253, 231, 37))]
    t = min(1.0, max(0.0, t))
    for (t0, c0), (t1, c1) in zip(stops[:-1], stops[1:]):
        if t <= t1:
            f = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
            r, g, b = (int(round(a + (b_ - a) * f)) for a, b_ in zip(c0, c1))
            return f"rgb({r},{g},{b})"
    return "rgb(253,231,37)"


class Canvas:
    def __init__(self, width: int = W, height: int = H, title: str = ""):
        self.w, self.h = width, height
        self.parts: list[str] = []
        self.title = title

    def add(self, s: str) -> None:
        self.parts.append(s)

    def text(self, x, y, s, size=13, anchor="middle", weight="normal",
             fill="#222", rot=0):
        tr = f' transform="rotate({rot},{x},{y})"' if rot else ""
        self.add(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
                 f'font-family="Helvetica,Arial,sans-serif" text-anchor="{anchor}" '
                 f'font-weight="{weight}" fill="{fill}"{tr}>{html.escape(str(s))}</text>')

    def rect(self, x, y, w, h, fill, stroke="none", sw=1.0):
        self.add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
                 f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')

    def polyline(self, pts: Sequence[tuple[float, float]], stroke="#fff", sw=2.5,
                 dash: str | None = None):
        if len(pts) < 2:
            return
        d = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        da = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<polyline points="{d}" fill="none" stroke="{stroke}" '
                 f'stroke-width="{sw}"{da}/>')

    def line(self, x1, y1, x2, y2, stroke="#888", sw=1.0, dash=None):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                 f'stroke="{stroke}" stroke-width="{sw}"{da}/>')

    def save(self, path: str) -> None:
        body = "\n".join(self.parts)
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}" '
               f'height="{self.h}" viewBox="0 0 {self.w} {self.h}">'
               f'<rect width="100%" height="100%" fill="white"/>{body}</svg>')
        with open(path, "w", encoding="utf-8") as f:
            f.write(svg)


def heatmap_panels(path: str, panels: list[dict], xlabels: Sequence[str],
                   ylabels: Sequence[str], xlabel: str, ylabel: str,
                   cbar_label: str, title: str, vmin: float, vmax: float) -> None:
    """panels: [{"title":str, "grid":[[v]] rows=y, cols=x,
                 "isoline":[(xi,yi)...] in cell coords, "band":[...]}]"""
    n = len(panels)
    cw = (W - PAD_L - PAD_R) / n
    c = Canvas(W, H, title)
    c.text(W / 2, 26, title, size=15, weight="bold")
    for pi, p in enumerate(panels):
        x0 = PAD_L + pi * cw
        pw = cw - 40
        ph = H - PAD_T - PAD_B
        rows, cols = len(ylabels), len(xlabels)
        dx, dy = pw / cols, ph / rows
        c.text(x0 + pw / 2, PAD_T - 12, p["title"], size=13, weight="bold")
        for i in range(rows):
            for j in range(cols):
                v = p["grid"][i][j]
                if v is None:
                    c.rect(x0 + j * dx, PAD_T + i * dy, dx, dy, "#dddddd")
                    continue
                t = (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5
                c.rect(x0 + j * dx, PAD_T + i * dy, dx, dy, _viridis(t),
                       stroke="#ffffff", sw=0.6)
                c.text(x0 + (j + 0.5) * dx, PAD_T + (i + 0.6) * dy,
                       f"{v:+.0f}", size=10,
                       fill="#ffffff" if t < 0.6 else "#222222")
        for key, dash, sw in (("band", "6,4", 1.8), ("isoline", None, 2.8)):
            for seg in p.get(key) or []:
                pts = [(x0 + (xi + 0.5) * dx, PAD_T + (yi + 0.5) * dy)
                       for xi, yi in seg]
                c.polyline(pts, stroke="#ff2d2d", sw=sw, dash=dash)
        for j, lab in enumerate(xlabels):
            c.text(x0 + (j + 0.5) * dx, H - PAD_B + 18, lab, size=11)
        if pi == 0:
            for i, lab in enumerate(ylabels):
                c.text(x0 - 10, PAD_T + (i + 0.6) * dy, lab, size=11, anchor="end")
            c.text(28, H / 2, ylabel, size=13, rot=-90)
        c.text(x0 + pw / 2, H - PAD_B + 42, xlabel, size=13)
    bx = W - PAD_R + 40
    for k in range(60):
        t = k / 59.0
        c.rect(bx, H - PAD_B - t * (H - PAD_T - PAD_B) - 6,
               18, (H - PAD_T - PAD_B) / 60 + 1, _viridis(t))
    c.text(bx + 9, PAD_T - 12, cbar_label, size=11)
    c.text(bx + 26, PAD_T + 6, f"{vmax:+.0f}", size=10, anchor="start")
    c.text(bx + 26, H - PAD_B, f"{vmin:+.0f}", size=10, anchor="start")
    c.save(path)


def lines(path: str, series: list[dict], xlabel: str, ylabel: str, title: str,
          xlim=None, ylim=None) -> None:
    """series: [{"x":[...], "y":[...], "label":str, "colour":str, "dash":str|None}]"""
    c = Canvas(W, H, title)
    c.text(W / 2, 26, title, size=15, weight="bold")
    xs = [v for s in series for v in s["x"]]
    ys = [v for s in series for v in s["y"]]
    if not xs:
        c.save(path)
        return
    x0v, x1v = (xlim or (min(xs), max(xs)))
    y0v, y1v = (ylim or (min(ys), max(ys)))
    if x1v == x0v:
        x1v = x0v + 1
    if y1v == y0v:
        y1v = y0v + 1
    pw, ph = W - PAD_L - PAD_R, H - PAD_T - PAD_B

    def px(x):
        return PAD_L + (x - x0v) / (x1v - x0v) * pw

    def py(y):
        return PAD_T + ph - (y - y0v) / (y1v - y0v) * ph

    c.rect(PAD_L, PAD_T, pw, ph, "none", stroke="#bbb")
    for k in range(5):
        yv = y0v + (y1v - y0v) * k / 4
        c.line(PAD_L, py(yv), PAD_L + pw, py(yv), "#eee")
        c.text(PAD_L - 8, py(yv) + 4, f"{yv:.2f}", size=10, anchor="end")
    for k in range(6):
        xv = x0v + (x1v - x0v) * k / 5
        c.text(px(xv), H - PAD_B + 18, f"{xv:.0f}", size=10)
    for si, s in enumerate(series):
        c.polyline([(px(x), py(y)) for x, y in zip(s["x"], s["y"])],
                   stroke=s.get("colour", "#3b528b"), sw=2.2, dash=s.get("dash"))
        c.text(W - PAD_R + 12, PAD_T + 16 + si * 18, s.get("label", ""), size=11,
               anchor="start", fill=s.get("colour", "#3b528b"))
    c.text(W / 2, H - PAD_B + 42, xlabel, size=13)
    c.text(28, H / 2, ylabel, size=13, rot=-90)
    c.save(path)
