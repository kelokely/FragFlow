#!/usr/bin/env python3
"""
Figure 4 -- the reference, and why one transfer free energy is not enough.

(a) Far-field benzene free-energy profile against reduced depth zeta, for all eight TREK-1 POPC
    replicates drawn individually, so that the reproducible shape and the varying level are both
    visible. Drawing individual replicates rather than a mean and a band is deliberate: the
    between-replicate spread is the quantity the panel reports.
(b) The five-probe transfer-free-energy ladder, spanning 9.9 kcal/mol.
(c) Cosolvent occupancy in the bilayer band against the ladder, with each point's own
    no-preference baseline drawn beneath it.

Reads   ../data/perrep/{apo,holo}_rep{1..4}.npz   per-replicate far-field probe depths
        ../data/fig1c_boxes.json                  per-system box geometry for panel (c)
Writes  fig4_reference.pdf, fig4_reference.png, fig4_reference.csv (Data File S1)

Two categorical hues are used and no more: apo against holo in panel (a), converged against pending
in panels (b) and (c).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

RT = 0.0019872041 * 303.15
ZB = (1.32, 1.53)
W = 0.10
EDG = np.arange(0.0, 1.6001, W)
CEN = 0.5 * (EDG[:-1] + EDG[1:])
NB = 11

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8a85"

# --- the five-probe ladder: each probe measured by two independent umbrella directions -----------
LADDER = [("propane", -3.95, True), ("benzene", -3.34, True), ("phenol", +0.22, True),
          ("2-propanol", +2.41, False), ("acetamide", +5.91, False)]
# Cosolvent occupancy in the band |z - z_mid| < 11 A, for the three converged probes.
# Phenol's "3.2x" is the odds ratio against uniform, (f/unif)/((1-f)/(1-unif)) = 3.21. It is an
# odds ratio and not a ratio of occupancies, and it is written as one wherever it is quoted.
OCC = {"benzene": 75.00, "phenol": 40.15, "acetamide": 3.00}

# The no-preference line is one number per system, not one number overall. The 22.0 A band is a
# fixed width, so the fraction of the box it covers depends on that box: the three systems these
# points were measured in have mean equilibrated box heights of 111.17, 126.80 and 128.59 A, giving
# geometric baselines of 19.79 %, 17.35 % and 17.11 %. Benzene was measured in GPR40, phenol and
# acetamide in TREK-1 apo. The baselines are read from the artifact rather than written as literals,
# so the panel cannot drift back to a single line.
import json as _json
import os as _os
_BOXES = _json.load(open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                       "..", "data", "fig1c_boxes.json")))
NO_PREF = {k: _BOXES[k]["no_preference_pct"] for k in OCC}
NO_PREF_LO, NO_PREF_HI = _BOXES["_range_pct"]
assert abs(NO_PREF["benzene"] - NO_PREF["phenol"]) > 2.0, (
    "the three baselines have converged; check the artifact before drawing them as one line")


def _odds(pct):
    p = pct / 100.0
    return p / (1.0 - p)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 7.6, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": "#c9c9c4",
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def profile(z):
    cnt = np.histogram(z, bins=EDG)[0].astype(float)
    n_zb = float(((z >= ZB[0]) & (z <= ZB[1])).sum())
    r_zb = n_zb / (ZB[1] - ZB[0])
    with np.errstate(divide="ignore"):
        return -RT * np.log((cnt / W) / r_zb)


fig = plt.figure(figsize=(7.4, 2.55))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0], wspace=0.42,
                      left=0.075, right=0.985, top=0.86, bottom=0.30)

# ---- (a) the far-field profile, per replicate ---------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
Pm = np.eye(NB) - np.ones((NB, NB)) / NB
raw = {}
for arm, colour in (("apo", BLUE), ("holo", ORANGE)):
    for rep in (1, 2, 3, 4):
        f = os.path.join(DATA, f"perrep/{arm}_rep{rep}.npz")
        g = profile(np.load(f, allow_pickle=True)["zeta"])[:NB]
        raw.setdefault(arm, []).append(g)
        ax.plot(CEN[:NB], g, color=colour, lw=1.1, alpha=0.75, zorder=2)
allg = np.array(raw["apo"] + raw["holo"])
lvl = allg[:, :5].mean(1)
ax.axhline(0.0, color=INK3, lw=0.7, ls=(0, (4, 3)), zorder=1)
ax.annotate(r"phosphate plane  $\zeta$=1", xy=(1.0, -0.15), fontsize=6.2, color=INK3,
            ha="center", va="top", rotation=0)
ax.axvline(1.0, color="#d8d7d2", lw=0.9, zorder=0)
ax.set_xlabel(r"reduced depth  $\zeta$ = 2·depth/$\langle$P–P$\rangle$", fontsize=7.6, labelpad=2)
ax.set_ylabel(r"$\Delta G_{\rm ref}(\zeta)$  (kcal mol$^{-1}$)", fontsize=7.6)
ax.set_title("a   the reference, 8 replicates", fontsize=8.0, loc="left", pad=5, color=INK)
ax.set_xlim(0, 1.35)
ax.grid(axis="y", color="#eeede9", lw=0.6)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
sig_raw = allg.std(0, ddof=1).mean()
sig_prj = (allg @ Pm).std(0, ddof=1)[:5].mean()
ax.text(0.97, 0.06,
        f"level spread  {lvl.std(ddof=1):.2f}\nshape spread  {sig_prj:.3f}",
        transform=ax.transAxes, fontsize=6.6, color=INK2, linespacing=1.4,
        va="bottom", ha="right",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#e2e1dc", lw=0.6))
ax.legend(handles=[Line2D([], [], color=BLUE, lw=1.6, label="apo ×4"),
                   Line2D([], [], color=ORANGE, lw=1.6, label="holo ×4")],
          frameon=False, fontsize=6.6, loc="upper left", handlelength=1.6, borderpad=0.2)

# ---- (b) the ladder -----------------------------------------------------------------------------
ax = fig.add_subplot(gs[0, 1])
names = [n for n, _v, _c in LADDER]
vals = [v for _n, v, _c in LADDER]
cols = [BLUE if c else INK3 for _n, _v, c in LADDER]
y = np.arange(len(LADDER))[::-1]
ax.barh(y, vals, height=0.52, color=cols, zorder=2)
ax.axvline(0, color=INK2, lw=0.8, zorder=3)
for yy, v, n in zip(y, vals, names):
    if v < 0:                       # inside the bar, white on the fill — no tick-label collision
        ax.text(v + 0.22, yy, f"{v:+.2f}", va="center", ha="left",
                fontsize=6.8, color="white", fontweight="bold", zorder=4)
    else:
        ax.text(v + 0.30, yy, f"{v:+.2f}", va="center", ha="left", fontsize=6.8, color=INK)
ax.set_yticks(y)
ax.set_yticklabels(names, fontsize=7.0)
ax.set_xlabel(r"$\Delta G$ mid-plane  (kcal mol$^{-1}$)", fontsize=7.6, labelpad=2)
ax.set_title("b   probe ladder, 9.9 kcal mol$^{-1}$", fontsize=8.0, loc="left", pad=5, color=INK)
ax.set_xlim(-5.0, 7.9)
ax.grid(axis="x", color="#eeede9", lw=0.6)
ax.set_axisbelow(True)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)

# ---- (c) occupancy against the ladder ----------------------------------------------------------
ax = fig.add_subplot(gs[0, 2])
# Each point carries its own geometric baseline, because each was measured in a different box.
# The band spans the three; the tick under each point is that point's own.
ax.axhspan(NO_PREF_LO, NO_PREF_HI, color=INK3, alpha=0.13, lw=0, zorder=0)
# The three baselines differ by 2.7 points on a 0-92 axis, which is about 3 % of the panel
# height -- visible as three ticks at different heights, but not readable as values. The range is
# therefore written out, so the reader gets the number rather than squinting at the gap.
ax.annotate("no depth preference: %.1f-%.1f %%\n(each box's own geometry)"
            % (NO_PREF_LO, NO_PREF_HI), xy=(-4.3, NO_PREF_HI),
            xytext=(0, 4), textcoords="offset points",
            fontsize=6.2, color=INK3, ha="left", va="bottom", linespacing=1.25)
for n, occ in OCC.items():
    dg = dict((a, b) for a, b, _c in LADDER)[n]
    ax.plot([dg - 0.42, dg + 0.42], [NO_PREF[n]] * 2, color=INK3, lw=1.1,
            ls=(0, (3, 2)), zorder=2, solid_capstyle="butt")
    ax.plot([dg], [occ], marker="o", ms=6.5, mfc=ORANGE if n == "phenol" else BLUE,
            mec="white", mew=1.2, zorder=3, ls="none")
    ax.annotate(n, (dg, occ), textcoords="offset points", xytext=(0, 9),
                fontsize=6.8, color=INK, ha="center")
_PHEN_ODDS = _odds(OCC["phenol"]) / _odds(NO_PREF["phenol"])
ax.annotate(r"$\Delta G$=+0.22, yet" "\n" + r"odds %.1f× uniform" % _PHEN_ODDS,
            xy=(0.22, OCC["phenol"]), xytext=(2.2, 63), textcoords="data",
            fontsize=6.4, color=ORANGE, ha="left", va="center", linespacing=1.35,
            arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.7, shrinkA=2, shrinkB=5))
ax.set_xlabel(r"$\Delta G$ mid-plane  (kcal mol$^{-1}$)", fontsize=7.6, labelpad=2)
ax.set_ylabel("occupancy in bilayer band  (%)", fontsize=7.6)
ax.set_title("c   rank is not enough", fontsize=8.0, loc="left", pad=5, color=INK)
ax.set_xlim(-4.6, 7.2)
ax.set_ylim(-4, 92)
ax.grid(axis="y", color="#eeede9", lw=0.6)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

for ext in ("pdf", "png"):
    fig.savefig(f"fig4_reference.{ext}", dpi=400, bbox_inches="tight", facecolor="white")
print("wrote fig4_reference.pdf / .png")
print(f"  panel a: level spread {lvl.std(ddof=1):.3f}, shape spread {sig_prj:.4f} kcal/mol")
print("  panel c baselines: " + "  ".join(
    f"{n} {NO_PREF[n]:.2f}% ({_BOXES[n]['mean_box_z_A']:.1f} A, {_BOXES[n]['system']})"
    for n in OCC))
print(f"  panel c: phenol odds {_PHEN_ODDS:.2f}x against its own baseline")

with open("fig4_reference.csv", "w") as fh:
    fh.write("panel,series,x,y\n")
    for arm in ("apo", "holo"):
        for i, g in enumerate(raw[arm], 1):
            for z, v in zip(CEN[:NB], g):
                fh.write(f"a,{arm}_rep{i},{z:.2f},{v:.4f}\n")
    for n, v, c in LADDER:
        fh.write(f"b,{n},,{v:+.2f}\n")
    for n, occ in OCC.items():
        fh.write(f"c,{n},{dict((a, b) for a, b, _ in LADDER)[n]:+.2f},{occ:.2f}\n")
    for n in OCC:
        fh.write(f"c,{n}_no_preference,{dict((a, b) for a, b, _ in LADDER)[n]:+.2f},"
                 f"{NO_PREF[n]:.2f}\n")
print("wrote fig4_reference.csv")
