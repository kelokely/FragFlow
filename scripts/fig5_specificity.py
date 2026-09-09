#!/usr/bin/env python3
"""
Figure 5 -- specificity, measured where the premise locates the failure.

Firing rate of the frozen detection rule (n_obs >= 1 AND ddG <= tau) against tau, for the
depth-resolved and the scalar bulk-water reference computed on identical voxels, four populations,
three proteins.

The four populations are A, bulk water; B, far-field acyl core; C, non-site acyl core; and D, the
validated site. Population C carries the same definition in all three systems: not within 5 A of any
detected cavity. Population B is the population the depth reference is fit on, so its depth-resolved
curve is circular by construction and is marked as such on the panel rather than plotted silently.

On TREK-1 the true-positive population is HS6 (zeta +0.55) and HS8 (zeta -0.34). It is not the ML335
site: three independent geometries around ML335 -- the 23 ligand atoms, the ligand centroid, and an
independent landmark 0.07 A from it -- each return zero benzene-visited voxels, and the site is
recorded as closed to two chemistries (measured ddG +0.77 and +2.28). HS1, at zeta +1.81, lies
outside the core mask. Section 3.6 of the manuscript gives the reasoning; the panel label states it.

Reads   ../data/fp2_frozen.json         firing rates, parsed from the measurement
        ../data/site_depth_frame.json   reduced site depth in the membrane frame
Writes  fig5_specificity.pdf, fig5_specificity.png, fig5_specificity.csv (Data File S2)

Color carries population and line style carries reference, so identity is never color alone and the
reader's primary comparison, solid against dashed within one color, is the one the figure leads
with. Only C and D take categorical hues; A and B are context and are drawn in neutral gray. The y
axis is logarithmic because the claim is a ratio of about a thousandfold. A firing rate of 0.00 %
cannot be drawn on a log axis, so those points are plotted at a labeled floor with an open downward
triangle meaning that no voxel fired.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

TAU = np.array([-0.5, -1.0, -1.5, -2.0, -2.5, -3.0])
FLOOR = 0.004          # log-axis floor; below the 1-voxel resolution of every population here
TAU_FROZEN = -2.0      # the manifest's frozen threshold

# population -> (depth-resolved, scalar), percentages, in TAU order
DATA = {
    # n = 8: all eight production replicates, 6008 frames, 546 728 probe observations.
    "AAC": dict(
        site="cav39", zeta=("AAC", "cav39"), nvox=dict(A=33285, B=43019, C=3690, D=8),
        A=([4.63, 0.23, 0.00, 0.00, 0.00, 0.00], [4.57, 0.20, 0.00, 0.00, 0.00, 0.00]),
        B=([0.06, 0.00, 0.00, 0.00, 0.00, 0.00], [47.01, 37.07, 25.64, 13.92, 0.12, 0.00]),
        C=([5.75, 1.33, 0.24, 0.05, 0.00, 0.00], [74.36, 67.59, 60.19, 44.17, 9.78, 2.17]),
        D=([100.00, 75.00, 75.00, 37.50, 37.50, 0.00], [100.00, 100.00, 100.00, 75.00, 75.00, 62.50])),
    # n = 8: all eight production replicates, 4006 frames, 364 546 probe observations.
    # D reads 0 of 13 voxels at tau = -2.0 on this system: the pre-fixed threshold sits past the
    # site distribution here, which is why the figure reports the whole curve and not one point.
    "C5aR1": dict(
        site="9P2", zeta=("C5aR1", "cav3"), nvox=dict(A=32987, B=43384, C=2993, D=13),
        A=([8.79, 1.25, 0.02, 0.00, 0.00, 0.00], [8.79, 1.25, 0.01, 0.00, 0.00, 0.00]),
        B=([0.16, 0.00, 0.00, 0.00, 0.00, 0.00], [55.31, 45.97, 35.78, 22.72, 7.03, 0.00]),
        C=([8.32, 1.84, 0.47, 0.10, 0.00, 0.00], [84.70, 78.78, 71.13, 56.36, 28.30, 5.85]),
        D=([69.23, 46.15, 30.77, 0.00, 0.00, 0.00], [100.00, 100.00, 84.62, 84.62, 76.92, 69.23])),
    # n = 4: the four apo production replicates. The holo replicates are a different state and
    # are not pooled with them.
    "TREK-1": dict(
        site="HS6, HS8", zeta=("TREK-1", "HS6"), nvox=dict(A=22481, B=28104, C=5009, D=29),
        A=([3.08, 3.08, 1.52, 0.02, 0.00, 0.00], [3.08, 3.08, 0.16, 0.01, 0.00, 0.00]),
        B=([0.23, 0.01, 0.00, 0.00, 0.00, 0.00], [50.00, 50.00, 40.95, 34.95, 22.99, 4.49]),
        C=([10.12, 2.93, 1.10, 0.44, 0.24, 0.04], [74.89, 74.89, 68.34, 63.61, 53.84, 29.39]),
        D=([24.14, 20.69, 13.79, 6.90, 0.00, 0.00], [44.83, 44.83, 41.38, 41.38, 37.93, 31.03])),
}

POPS = [
    ("C", "#2a78d6", 1.9, "C  non-site acyl core", 3.2),
    ("D", "#eb6834", 1.9, "D  validated site",     3.3),
    ("A", "#6e6e68", 1.1, "A  bulk water",         2.0),
    ("B", "#9a9a94", 1.1, "B  far-field core ⚠ circular", 2.1),
]

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a85"
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 7.6, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": "#c9c9c4",
    "figure.facecolor": "white", "axes.facecolor": "white",
})

# ==================================================================================================
# DATA above is a transcription: 3 systems x 4 populations x (1 voxel count + 2 references x 6
# thresholds), 156 numbers, and it draws the figure that carries the paper's central claim. It is
# therefore checked against the measurement before anything is drawn. fp2_frozen.json was produced
# by parsing the analysis output rather than by reading numbers off it, so the chain runs
# measurement -> parse -> artifact -> assert -> figure with no hand transcription in it.
# ==================================================================================================
import json
import os

FROZEN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "data", "fp2_frozen.json")


def check_data_against_frozen():
    # `python -O` strips every assert, which would turn this entire check into a no-op and draw
    # the figure with whatever DATA happens to say. Refuse instead: an unverified figure is worse
    # than no figure.
    if not __debug__:
        raise RuntimeError(
            "this figure's data check is built on assertions and python -O strips them; "
            "run without -O so the numbers are verified before anything is drawn.")

    with open(FROZEN) as fh:
        ref = json.load(fh)
    cur = ref["current"]

    def same(a, b):
        return len(a) == len(b) and all(abs(float(x) - float(y)) < 5e-3 for x, y in zip(a, b))

    assert set(DATA) == set(cur), (
        f"SYSTEM SET DRIFT: DATA has {sorted(DATA)}, frozen file has {sorted(cur)}.")

    for name in sorted(DATA):
        d, c = DATA[name], cur[name]
        for pop in "ABCD":
            ok = c["pops"][pop]
            assert d["nvox"][pop] == ok["nvox"], (
                f"DATA DRIFT: {name} nvox[{pop}] is {d['nvox'][pop]}, {ok['source']} says "
                f"{ok['nvox']}.")
            for arm, mine in (("depth", d[pop][0]), ("scalar", d[pop][1])):
                assert len(mine) == len(TAU), f"{name}/{pop}/{arm}: {len(mine)} of {len(TAU)} taus"
                assert same(mine, ok[arm]), (
                    f"DATA DRIFT: {name} {pop} {arm} is {list(mine)}, {ok['source']} says "
                    f"{ok[arm]}. DATA is a transcription; the job output is the measurement. "
                    f"Reconcile before drawing.")

    print("  DATA verified against %s: %s" % (
        os.path.basename(FROZEN),
        ", ".join(f"{n} n={cur[n]['reps']}" for n in sorted(cur))), flush=True)


check_data_against_frozen()


# ==================================================================================================
# Site depth is read from an artifact rather than written as a literal, and from the same artifact
# every other depth in the package comes from. The check above guards the firing rates and does not
# look at the depths, so they get their own source and their own validation.
# ==================================================================================================
def resolve_zeta():
    with open(FROZEN.replace("fp2_frozen.json", "site_depth_frame.json")) as fh:
        sdf = json.load(fh)
    if not sdf.get("_validation_replicate1", {}).get("passed"):
        raise AssertionError("site_depth_frame.json did not pass its own replicate-1 "
                             "validation; refusing to label sites with depths from it")
    for name, d in DATA.items():
        key = d["zeta"]
        if not isinstance(key, tuple):
            raise AssertionError(f"{name}: zeta must name (system, cavity) in the artifact, "
                                 f"not carry a literal {key!r}")
        sysn, cav = key
        d["zeta"] = round(sdf[sysn][cav]["mean_of_replicates"], 2)
    print("  site depths from site_depth_frame.json: "
          + ", ".join(f"{n} {d['zeta']:+.2f}" for n, d in DATA.items()), flush=True)


resolve_zeta()

fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.15), sharey=True)


def fmt_x(v):
    if v == 0:
        return "0"
    return f"{v:,.0f}" if v >= 100 else f"{v:.1f}"


for ax, (name, d) in zip(axes, DATA.items()):
    ax.axvspan(TAU_FROZEN - 0.06, TAU_FROZEN + 0.06, color="#ecebe5", zorder=0, lw=0)
    ax.axhspan(0.0025, FLOOR * 1.3, color="#f5f4f0", zorder=0, lw=0)

    for key, color, lw, _lab, _z in POPS:
        for vals, style, dash in ((d[key][1], "scalar", (3.2, 1.7)), (d[key][0], "depth", None)):
            y = np.array(vals, float)
            drawn = np.where(y <= 0, FLOOR, y)
            ax.plot(TAU, drawn, color=color, lw=lw, ls="-",
                    dashes=dash if dash else (None, None),
                    solid_capstyle="round", zorder=_z,
                    marker="o" if style == "depth" else None,
                    ms=3.0, mfc="white", mew=1.0, mec=color,
                    alpha=1.0 if key in ("C", "D") else 0.85)
            zero = y <= 0
            if zero.any():
                ax.plot(TAU[zero], np.full(zero.sum(), FLOOR), ls="none", marker="v",
                        ms=3.6, mfc="white", mew=0.9, mec=color, zorder=_z + 0.1)

    cs, cd = d["C"][1][3], d["C"][0][3]
    ds, dd = d["D"][1][3], d["D"][0][3]
    bbox = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.88)
    ax.annotate(f"{cs:.1f}%", (TAU_FROZEN, cs), textcoords="offset points", xytext=(7, 7),
                fontsize=7.3, color=INK, fontweight="bold", ha="left", va="bottom",
                bbox=bbox, zorder=6)
    ax.annotate(f"{cd:.2f}%", (TAU_FROZEN, max(cd, FLOOR)), textcoords="offset points",
                xytext=(7, -8), fontsize=7.3, color=INK, fontweight="bold",
                ha="left", va="top", bbox=bbox, zorder=6)

    # Table 1 quotes this ratio from the EXACT voxel counts. Dividing the rounded rates gives a
    # different number -- in the submitted version 50.00/0.03 printed 1667x here against the
    # table's 50.00/(1/3686) = 1843x. Recover the count and divide the same way the table does.
    _nC = d["nvox"]["C"]
    _kC = round(cd * _nC / 100.0)
    disc_d = (dd / (100.0 * _kC / _nC)) if _kC else float("inf")
    zlab = (f"$\\zeta$ = {d['zeta']:+.2f}" if name != "TREK-1"
            else r"$\zeta$ +0.55, $-$0.34   (ML335: no benzene)")
    ax.set_title(f"{name}\n{d['site']}   {zlab}\n"
                 f"site ÷ non-site:  {fmt_x(disc_d)}× vs {ds / cs:.2f}×",
                 fontsize=8.3, color=INK, pad=6, linespacing=1.5)

    ax.set_yscale("log")
    ax.set_ylim(0.0025, 320)
    ax.set_xlim(-0.32, -3.18)
    ax.set_xticks(TAU)
    ax.set_xticklabels([f"{t:.1f}" for t in TAU], fontsize=7.0)
    ax.set_xlabel(r"threshold $\tau$  (kcal mol$^{-1}$)", fontsize=7.8, labelpad=2)
    ax.grid(axis="y", color="#e9e8e3", lw=0.6, zorder=-1)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

axes[0].set_ylabel("voxels firing  (%)", fontsize=7.8)
axes[0].set_yticks([0.01, 0.1, 1, 10, 100])
axes[0].set_yticklabels(["0.01", "0.1", "1", "10", "100"], fontsize=7.0)
axes[0].annotate("0", xy=(0, FLOOR), xycoords=("axes fraction", "data"),
                 xytext=(-5, 0), textcoords="offset points",
                 fontsize=7.0, color=INK3, ha="right", va="center", annotation_clip=False)

from matplotlib.patches import Patch
handles = [Line2D([], [], color=c, lw=lw + 0.3, label=lab) for _k, c, lw, lab, _z in POPS]
handles += [Line2D([], [], color=INK2, lw=1.7, marker="o", ms=3.2, mfc="white",
                   mec=INK2, label="depth-resolved reference"),
            Line2D([], [], color=INK2, lw=1.7, dashes=(3.2, 1.7), label="scalar bulk-water reference"),
            Patch(facecolor="#ecebe5", edgecolor="none", label=r"frozen rule  $\tau$ = $-$2.0"),
            Line2D([], [], color=INK3, lw=0, marker="v", ms=3.6, mfc="white", mec=INK3,
                   label="no voxel fired (drawn at 0)")]
leg = fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
                 fontsize=7.0, handlelength=2.3, columnspacing=1.9,
                 handletextpad=0.7, bbox_to_anchor=(0.5, 0.004))
for t in leg.get_texts():
    t.set_color(INK2)

fig.subplots_adjust(left=0.075, right=0.995, top=0.775, bottom=0.235, wspace=0.10)
for ext in ("pdf", "png"):
    fig.savefig(f"fig5_specificity.{ext}", dpi=400, bbox_inches="tight", facecolor="white")
print("wrote fig5_specificity.pdf / .png")

with open("fig5_specificity.csv", "w") as fh:
    fh.write("protein,site,zeta,population,n_voxels,reference," +
             ",".join(f"tau={t:+.1f}" for t in TAU) + "\n")
    for name, d in DATA.items():
        for key, _c, _lw, lab, _z in POPS:
            for arm, vals in (("depth-resolved", d[key][0]), ("scalar", d[key][1])):
                # site is quoted: TREK-1 names two sites, and the comma inside that field would
                # otherwise split the row into 13 fields where every other row has 12.
                fh.write(f"{name},\"{d['site']}\",{d['zeta']:+.2f},\"{lab}\",{d['nvox'][key]},{arm}," +
                         ",".join(f"{v:.2f}" for v in vals) + "\n")
print("wrote fig5_specificity.csv  (the table view)")

for name, d in DATA.items():
    cd, cs = d["C"][0][3], d["C"][1][3]
    dd, ds = d["D"][0][3], d["D"][1][3]
    # exact-count ratio, identical to the panel title and to Table 1 -- NOT dd/cd on rounded rates
    _nC = d["nvox"]["C"]; _kC = round(cd * _nC / 100.0)
    _disc = (dd / (100.0 * _kC / _nC)) if _kC else float("inf")
    print(f"  {name:7s} tau=-2.0  C {cs:6.2f}% -> {cd:5.2f}% ({_kC} voxels of {_nC})   "
          f"discrimination depth {_disc:8.1f}x  scalar {ds/cs:.2f}x")
