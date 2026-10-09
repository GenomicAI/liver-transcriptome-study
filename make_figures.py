#!/usr/bin/env python3
"""Draw the paper's three figures as vector PDFs from results/*.csv (written by liver_transcriptome_study.R).

    python3 make_figures.py        # from this folder; writes figures/fig_variance.pdf, fig_trends.pdf, fig_prediction.pdf

Needs Python 3.9+ and matplotlib (written and checked with 3.11); nothing else.  Nothing is typed by hand: every
plotted value is read from results/, and the numbers the paper's text quotes are asserted at the end, so a stale
results/ folder fails loudly.

Each figure is drawn at the size it is printed.  The paper's text block is 6.4 in wide (8.5 in page, 1.05 in
margins) and the PDFs go in with \\includegraphics[width=\\textwidth], that is at scale 1, so every label is exactly
10 pt on the page, the minimum the conference allows.  The font is STIX (Times-like, ships with matplotlib), embedded
as TrueType (pdf.fonttype 42), not as Type 3.  Colours are the Okabe-Ito palette; every series also has its own marker
and line style, so the figures still read in greyscale.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")                    # file output only (the macosx backend would resize the figure)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

W = 6.4                                  # inches: the paper's text width
FONT = 10                                # points: the conference minimum, used for every label
SKY, VERM, BLUE, ORANGE, GREEN, PURPLE = "#56B4E9", "#D55E00", "#0072B2", "#E69F00", "#009E73", "#CC79A7"
GREY = "#404040"
plt.rcParams.update({
    "font.family": "STIXGeneral", "font.size": FONT, "mathtext.fontset": "stix",
    "axes.titlesize": FONT, "axes.labelsize": FONT, "xtick.labelsize": FONT, "ytick.labelsize": FONT,
    "legend.fontsize": FONT,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 3, "ytick.major.size": 3, "xtick.major.pad": 2.5, "ytick.major.pad": 2.5,
    "axes.titlepad": 4, "axes.labelpad": 2.5,
    "grid.color": "#dfdfdf", "grid.linewidth": 0.5,
    "pdf.fonttype": 42,
})


def read(name):
    with open(os.path.join(RES, name), newline="") as fh:
        return list(csv.DictReader(fh))


# Short class names used on the axes (module keys are the R script's names).
SHORT = {
    "ecm": "Fibrosis / ECM",
    "ductular": "Ductular reaction",
    "inflam": "Inflammation",
    "lipid": "Lipid / steatosis",
    "synthetic": "Plasma proteins",
    "metabolic": "CYP metabolism",
    "background": "Background",
}
ORDER = ["ecm", "ductular", "inflam", "lipid", "synthetic", "metabolic", "background"]
STAGES = ["Healthy", "F0-F1", "F2", "F3", "F4"]
STAGE_TICKS = ["H", "F0–F1", "F2", "F3", "F4"]

mv = {r["module"]: r for r in read("module_variance.csv")}
disp = read("dispersion_by_expression.csv")
trend = {(r["module"], r["stage"]): r for r in read("module_stage_trend.csv")}
pred = {r["model"]: r for r in read("transcriptome_prediction.csv")}
bm = {r["module"]: r for r in read("transcriptome_prediction_by_module.csv")}
f = float


# ---- helpers ----------------------------------------------------------------------------------
def new_figure(height):
    return plt.figure(figsize=(W, height))


def axes_in(fig, left, bottom, width, height):
    """Axes placed by inches from the lower-left corner of the figure."""
    fw, fh = fig.get_size_inches()
    return fig.add_axes([left / fw, bottom / fh, width / fw, height / fh])


def style(ax, xgrid=False, ygrid=False):
    ax.set_axisbelow(True)
    if xgrid:
        ax.xaxis.grid(True)
    if ygrid:
        ax.yaxis.grid(True)


def class_axis(ax):
    """Seven gene classes down the vertical axis, first class at the top, names as tick labels."""
    ax.set_ylim(7.6, 0.4)
    ax.set_yticks(range(1, 8))
    ax.set_yticklabels([SHORT[m] for m in ORDER])
    ax.tick_params(axis="y", length=0, pad=4)


def legend_style():
    """Legends inside a plot get a faint white patch, so gridlines do not run through the words."""
    return dict(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9, fancybox=False,
                handletextpad=0.4, labelspacing=0.15, borderaxespad=0.3, borderpad=0.2)


def check_figure(fig, name):
    """Every visible label is at least FONT pt and lies inside the figure, so nothing is clipped at the edge."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = []
    for ax in fig.axes:
        texts += [ax.title, ax.xaxis.label, ax.yaxis.label, *ax.get_xticklabels(), *ax.get_yticklabels(), *ax.texts]
        if ax.get_legend():
            texts += ax.get_legend().get_texts()
    for lg in fig.legends:
        texts += lg.get_texts()
    shown = [t for t in texts if isinstance(t, Text) and t.get_visible() and t.get_text().strip()]
    assert shown, name
    right, top = fig.bbox.width, fig.bbox.height
    for t in shown:
        assert t.get_fontsize() >= FONT, "%s: %r is %.1f pt" % (name, t.get_text(), t.get_fontsize())
        bb = t.get_window_extent(renderer)
        assert bb.x0 >= 0 and bb.y0 >= 0 and bb.x1 <= right and bb.y1 <= top, \
            "%s: %r is cut off at the figure edge" % (name, t.get_text())
    return len(shown)


def save(fig, name):
    n = check_figure(fig, name)
    path = os.path.join(OUT, name + ".pdf")
    fig.savefig(path, format="pdf", metadata={"Title": name, "Creator": "make_figures.py (matplotlib)",
                                               "CreationDate": None})        # no date, so reruns give identical files
    plt.close(fig)
    print("wrote %-28s %.2f x %.2f in, %d labels, all >= %d pt and inside the figure"
          % (os.path.relpath(path, HERE), fig.get_figwidth(), fig.get_figheight(), n, FONT))


# ---- Figure 1: variance by gene class ---------------------------------------------------------
def fig_variance():
    fig = new_figure(2.2)
    # (a) median BCV by gene class, healthy vs diseased
    ax = axes_in(fig, 1.17, 0.47, 1.75, 1.5)
    hv = [f(mv[m]["bcv_healthy"]) for m in ORDER]
    dv = [f(mv[m]["bcv_diseased"]) for m in ORDER]
    h, off = 0.29, 0.17
    ax.barh([i - off for i in range(1, 8)], hv, height=h, color=SKY, linewidth=0, label="Healthy")
    ax.barh([i + off for i in range(1, 8)], dv, height=h, color=VERM, linewidth=0, label="Diseased")
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    class_axis(ax)
    style(ax, xgrid=True)
    ax.set_xlabel("BCV (median)")
    ax.set_title("(a) Gene class")
    ax.legend(loc="lower right", handlelength=0.9, handleheight=0.7, **legend_style())

    # (b) median BCV by mean expression level
    ax = axes_in(fig, 3.62, 0.47, 2.7, 1.5)
    x = range(1, len(disp) + 1)
    ax.plot(x, [f(r["bcv_healthy"]) for r in disp], color=SKY, marker="o", ms=4, lw=0.9)
    ax.plot(x, [f(r["bcv_diseased"]) for r in disp], color=VERM, marker="s", ms=3.6, lw=0.9)
    ax.set_xlim(0.6, 7.4)
    ax.set_xticks(list(x))
    ax.set_xticklabels([r["bin"].replace("-", "–") for r in disp])
    ax.set_ylim(0.2, 0.34)
    ax.set_yticks([0.2, 0.25, 0.3])
    ax.set_yticklabels(["0.20", "0.25", "0.30"])
    style(ax, ygrid=True)
    ax.set_xlabel("Mean expression ($\\mathrm{log}_2$ CPM)")
    ax.set_ylabel("BCV (median)")
    ax.set_title("(b) Expression level")
    save(fig, "fig_variance")


# ---- Figure 2: module score by fibrosis stage ---------------------------------------------------
def fig_trends():
    fig = new_figure(2.1)
    # (module key, colour, marker, legend text) for the two modules in each panel
    panels = [
        ("(a) Early, saturating", [("inflam", BLUE, "o", "Inflammation"), ("lipid", ORANGE, "s", "Lipid / steatosis")]),
        ("(b) Late, accelerating", [("ecm", VERM, "^", "Fibrosis / ECM"), ("ductular", PURPLE, "D", "Ductular reaction")]),
        ("(c) Steady or late loss", [("metabolic", SKY, "o", "CYP metabolism"), ("synthetic", GREEN, "s", "Plasma proteins")]),
    ]
    left, aw, gap = 0.57, 1.78, 0.2
    for k, (title, series) in enumerate(panels):
        ax = axes_in(fig, left + k * (aw + gap), 0.28, aw, 1.6)
        handles = []
        for j, (m, col, mk, lab) in enumerate(series):
            rows = [trend[(m, s)] for s in STAGES]
            med = [f(r["median"]) for r in rows]
            lo = [md - f(r["q1"]) for md, r in zip(med, rows)]
            hi = [f(r["q3"]) - md for md, r in zip(med, rows)]
            dx = (-0.06, 0.06)[j]                       # nudge the two modules apart so their error bars do not overlap
            ax.errorbar([i + dx for i in range(5)], med, yerr=[lo, hi], color=col, marker=mk, ms=3.8, lw=0.9,
                        elinewidth=0.6, capsize=1.6, capthick=0.6)
            ax.plot(range(5), [f(r["expected"]) for r in rows], color=col, lw=0.7, ls=(0, (4, 2.5)))
            handles.append(Line2D([], [], color=col, marker=mk, ms=3.8, lw=0.9, label=lab))
        ax.set_xlim(-0.35, 4.35)
        ax.set_xticks(range(5))
        ax.set_xticklabels(STAGE_TICKS)
        ax.set_ylim(-1.9, 3.6)
        ax.set_yticks([-1, 0, 1, 2, 3])
        if k:
            ax.tick_params(axis="y", labelleft=False)
        else:
            ax.set_ylabel("Module score ($\\mathrm{log}_2$)")
        style(ax, ygrid=True)
        ax.set_title(title)
        ax.legend(handles=handles, loc="upper left", handlelength=1.5, **legend_style())
    save(fig, "fig_trends")


# ---- Figure 3: whole-transcriptome prediction ---------------------------------------------------
def fig_prediction():
    fig = new_figure(2.62)
    # (a) per-liver r against the number of landmark genes
    ax = axes_in(fig, 0.57, 0.78, 2.1, 1.6)
    ks = (10, 25, 50, 100, 200, 400)
    ax.axhline(f(pred["Disease label"]["median_sample_r"]), color="gray", lw=1.0, ls=":")   # stage-only is 0.383, same line
    ax.text(0.985, 0.02, "labels only", transform=ax.transAxes, ha="right", va="bottom", color=GREY)
    ax.plot(ks, [f(pred["LM%d" % k]["median_sample_r"]) for k in ks], color=BLUE, marker="o", ms=4, lw=0.9,
            label="Landmarks")
    ax.plot(ks, [f(pred["LM%d+stage" % k]["median_sample_r"]) for k in ks], color=ORANGE, marker="s", ms=3.6, lw=0.9,
            ls="--", label="Landmarks + stage")
    ax.set_xscale("log")
    ax.set_xlim(7, 550)
    ax.set_xticks(ks)
    ax.set_xticklabels([str(k) for k in ks])
    ax.minorticks_off()
    ax.set_ylim(0.3, 0.86)
    ax.set_yticks([0.4, 0.5, 0.6, 0.7, 0.8])
    style(ax, ygrid=True)
    ax.set_xlabel("Landmark genes $K$")
    ax.set_ylabel("Median per-liver $r$")
    ax.set_title("(a) Landmark panels")
    ax.legend(loc="upper left", handlelength=1.5, **legend_style())

    # (b) median per-gene R^2 by gene class
    ax = axes_in(fig, 3.87, 0.78, 2.45, 1.6)
    y = range(1, 8)
    series = [("R2_label", ORANGE, "^", 5.5, "Disease label"), ("R2_stage", SKY, "s", 4.5, "Fibrosis stage"),
              ("R2_LM100", BLUE, "o", 4.5, "100 landmarks")]
    for z, (col_name, col, mk, ms, lab) in enumerate(series, start=2):
        ax.plot([f(bm[m][col_name]) for m in ORDER], y, ls="none", color=col, marker=mk, ms=ms, zorder=z)
    ax.plot([f(bm[m]["R2_oracle"]) for m in ORDER], y, ls="none", color="black", marker="|", ms=7, mew=1.2, zorder=5)
    ax.set_xlim(-0.05, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    class_axis(ax)
    style(ax, xgrid=True, ygrid=True)
    ax.set_xlabel("Median per-gene $R^2$")
    ax.set_title("(b) Gene class")
    handles = [Line2D([], [], ls="none", color=c, marker=m, ms=s, label=l) for _, c, m, s, l in series]
    handles.append(Line2D([], [], ls="none", color="black", marker="|", ms=7, mew=1.2, label="Ceiling"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, columnspacing=1.4, handletextpad=0.3,
               borderaxespad=0.15, numpoints=1)
    save(fig, "fig_prediction")


fig_variance()
fig_trends()
fig_prediction()

# ---- assertions: numbers quoted in the paper's text -------------------------------------------------
assert f(trend[("ductular", "F4")]["median"]) == 2.043          # "+2.04 only at F4"
assert f(trend[("ecm", "F4")]["median"]) == 2.359               # "+2.36 at F4"
assert f(trend[("lipid", "F2")]["median"]) == 1.398             # "+1.40"
assert f(trend[("lipid", "F4")]["median"]) == 0.899             # "+0.90 in cirrhosis"
assert f(trend[("inflam", "F0-F1")]["median"]) == 0.62          # "+0.62 log2 at F0-F1"
assert f(trend[("metabolic", "F4")]["median"]) == -1.284        # "-1.28 at F4"
assert f(trend[("synthetic", "F4")]["median"]) == -0.817        # "-0.82"
assert f(pred["LM100"]["median_sample_r"]) == 0.682             # "r = 0.68"
assert f(pred["LM10"]["median_sample_r"]) == 0.616              # "r = 0.62 with 10 landmarks"
assert f(pred["Disease label"]["median_sample_r"]) == 0.381     # "0.38"
assert f(pred["Fibrosis stage"]["median_sample_r"]) == 0.383    # "both 0.38"
assert f({r["bin"]: r for r in disp}[">10"]["bcv_diseased"]) == 0.303      # "increased to 0.30"
assert (mv["background"]["bcv_healthy"], mv["background"]["bcv_diseased"]) == ("0.249", "0.251")
vr = [f(mv[m]["var_ratio_D_over_H"]) for m in ORDER[:-1]]
assert (min(vr), max(vr)) == (2.31, 5.16)                       # "2.3--5.2 times more variable"
assert (f(bm["synthetic"]["R2_LM100"]), f(bm["ecm"]["R2_LM100"])) == (0.609, 0.815)      # "0.61--0.82"
assert (f(bm["ductular"]["R2_stage"]), f(bm["ductular"]["R2_label"])) == (0.591, 0.158)  # "59% vs 16%"
print("all numbers quoted in the text match results/")
