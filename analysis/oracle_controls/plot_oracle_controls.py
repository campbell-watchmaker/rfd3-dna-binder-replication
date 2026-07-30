#!/usr/bin/env python3
"""Figure for the ΔminPAE oracle control experiment.

Panels
  A, B   minPAE heatmap (protein x DNA target), one per oracle. Sequential
         single-hue blue ramp, INVERTED: low minPAE (a confident interface) is
         dark and salient; high minPAE (no interface) recedes toward the surface,
         which is the semantically right direction since absence-of-signal should
         recede. Each oracle keeps its own colorbar because their PAE scales are
         differently calibrated -- forcing a shared scale would imply a
         comparability that does not exist.
  C      ΔminPAE per specific TF, grouped by oracle. Positive = the TF's own site
         is read more confidently than the best off-target.
  D      minPAE spread across the panel for every protein, coloured by class.
         Spread needs no designated on-target, so it is the metric that lets
         non-specific binders and non-binders be compared on the same axis.

Colour: categorical slots 1-3 (blue/orange/aqua) for the three protein classes --
the first three slots are the set validated for all-pairs CVD separation.
Identity is never colour-alone: every class is also directly labelled.

Usage:
    python plot_oracle_controls.py --summary rf3=results/oracle_controls/rf3_summary.csv \
                                   --summary protenix=results/oracle_controls/protenix_summary.csv \
                                   --folds   rf3=results/oracle_controls/rf3_folds.csv \
                                   --folds   protenix=results/oracle_controls/protenix_folds.csv \
                                   --out figures/oracle_controls.png
"""
from __future__ import annotations
import argparse
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

# --- design-system parameters (see dataviz references/palette.md) -------------
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
TEXT_MUTED = "#8a8880"
GRID = "#e4e3df"
# sequential blue ramp, steps 100 -> 700
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
# categorical slots 1-3 (validated all-pairs in both modes)
CLASS_COLOR = {
    "specific": "#2a78d6",
    "nonspecific_binder": "#eb6834",
    "nonbinder": "#1baf7a",
}
CLASS_LABEL = {
    "specific": "sequence-specific TF",
    "nonspecific_binder": "non-specific duplex binder",
    "nonbinder": "non-binder",
}
CLASS_ORDER = ["specific", "nonspecific_binder", "nonbinder"]

# minPAE low = confident, so invert the ramp: dark = low = salient
CMAP = LinearSegmentedColormap.from_list("minpae", list(reversed(BLUE_RAMP)))

DNA_ORDER = ["zif268_site", "lambda_OL1", "ebox", "hd_taatta", "tata",
             "prnp", "scramble", "polygc"]
DNA_LABEL = {
    "zif268_site": "Zif268\nGCGTGGGCGT", "lambda_OL1": "λ OL1\n17 bp",
    "ebox": "E-box\nCACGTG", "hd_taatta": "HD\nTAATTA", "tata": "TATA\nTATAAAA",
    "prnp": "PRNP\n(our target)", "scramble": "scramble", "polygc": "poly-GC",
}


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", action="append", default=[], metavar="ORACLE=PATH")
    ap.add_argument("--folds", action="append", default=[], metavar="ORACLE=PATH")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    summaries, folds = {}, {}
    for spec in args.summary:
        o, _, p = spec.partition("=")
        if os.path.exists(p):
            summaries[o] = read_csv(p)
    for spec in args.folds:
        o, _, p = spec.partition("=")
        if os.path.exists(p):
            folds[o] = read_csv(p)
    if not folds:
        print("no per-fold CSVs found -- nothing to plot")
        return 1

    oracles = [o for o in ("rf3", "protenix") if o in folds]

    # protein order: by class, then by name, from whichever summary we have
    prot_class = {}
    for rows in list(summaries.values()) + list(folds.values()):
        for r in rows:
            prot_class[r["protein"]] = r["klass"]
    proteins = sorted(prot_class, key=lambda p: (CLASS_ORDER.index(prot_class[p]), p))

    on_target_of = {}
    for rows in summaries.values():
        for r in rows:
            if r.get("on_target"):
                on_target_of[r["protein"]] = r["on_target"]

    fig = plt.figure(figsize=(15.5, 11), facecolor=SURFACE)
    gs = fig.add_gridspec(2, len(oracles), height_ratios=[1.5, 1.0],
                          hspace=0.30, wspace=0.24,
                          left=0.115, right=0.975, top=0.90, bottom=0.085)

    # ---------------- Panels A/B: heatmaps ----------------
    for k, oracle in enumerate(oracles):
        ax = fig.add_subplot(gs[0, k])
        ax.set_facecolor(SURFACE)
        vals = {(r["protein"], r["dna_id"]): fnum(r["min_pae"]) for r in folds[oracle]}
        dna_present = [d for d in DNA_ORDER if any((p, d) in vals for p in proteins)]
        M = np.full((len(proteins), len(dna_present)), np.nan)
        for i, p in enumerate(proteins):
            for j, d in enumerate(dna_present):
                v = vals.get((p, d))
                if v is not None:
                    M[i, j] = v

        im = ax.imshow(M, cmap=CMAP, aspect="auto", interpolation="nearest")
        ax.set_xticks(range(len(dna_present)))
        ax.set_xticklabels([DNA_LABEL.get(d, d) for d in dna_present],
                           fontsize=7.5, color=TEXT_SECONDARY)
        ax.set_yticks(range(len(proteins)))
        ax.set_yticklabels(proteins, fontsize=9, color=TEXT_PRIMARY)
        # colour the y tick labels by class -> identity is redundant with position
        for lbl, p in zip(ax.get_yticklabels(), proteins):
            lbl.set_color(CLASS_COLOR[prot_class[p]])
        ax.set_title(f"{'AB'[k]}   minPAE — {oracle}", fontsize=11.5, weight="bold",
                     color=TEXT_PRIMARY, loc="left", pad=10)

        # direct labels: the number in every cell (8x8 is small enough to read)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isnan(M[i, j]):
                    ax.text(j, i, "·", ha="center", va="center",
                            fontsize=9, color=TEXT_MUTED)
                    continue
                rng = np.nanmax(M) - np.nanmin(M)
                dark = rng > 0 and (M[i, j] - np.nanmin(M)) / rng < 0.45
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7.4,
                        color="#ffffff" if dark else TEXT_SECONDARY)

        # ring the on-target cell for each specific TF (2px surface ring)
        for i, p in enumerate(proteins):
            on = on_target_of.get(p)
            if on in dna_present:
                j = dna_present.index(on)
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor=TEXT_PRIMARY, lw=2.0, zorder=5))
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
        cb = fig.colorbar(im, ax=ax, fraction=0.036, pad=0.02)
        cb.set_label("minPAE (Å) — lower = more confident interface",
                     fontsize=7.8, color=TEXT_SECONDARY)
        cb.ax.tick_params(labelsize=7, colors=TEXT_SECONDARY, length=0)
        cb.outline.set_visible(False)

    # ---------------- Panel C: ΔminPAE for specific TFs ----------------
    axc = fig.add_subplot(gs[1, 0])
    axc.set_facecolor(SURFACE)
    spec_prots = [p for p in proteins if prot_class[p] == "specific"]
    n_o = len(oracles)
    width = 0.8 / max(n_o, 1)
    any_delta = False
    for k, oracle in enumerate(oracles):
        d = {r["protein"]: fnum(r.get("delta_min_pae")) for r in summaries.get(oracle, [])}
        xs, ys = [], []
        for i, p in enumerate(spec_prots):
            v = d.get(p)
            if v is not None:
                xs.append(i + (k - (n_o - 1) / 2) * width)
                ys.append(v)
                any_delta = True
        if xs:
            axc.bar(xs, ys, width=width * 0.92,
                    color=["#2a78d6", "#eb6834"][k % 2],
                    label=oracle, zorder=3, linewidth=0)
    axc.axhline(0, color=TEXT_SECONDARY, lw=1.0, zorder=4)
    axc.set_xticks(range(len(spec_prots)))
    axc.set_xticklabels(spec_prots, fontsize=8.5, color=TEXT_PRIMARY)
    axc.set_ylabel("ΔminPAE (Å)", fontsize=9, color=TEXT_SECONDARY)
    axc.set_title("C   ΔminPAE per specific TF   (>0 = own site read best)",
                  fontsize=11.5, weight="bold", color=TEXT_PRIMARY, loc="left", pad=10)
    if any_delta:
        axc.legend(frameon=False, fontsize=8.5, labelcolor=TEXT_SECONDARY, ncol=n_o)
    axc.grid(axis="y", color=GRID, lw=0.8, zorder=0)
    axc.set_axisbelow(True)
    for s in ("top", "right", "bottom", "left"):
        axc.spines[s].set_visible(False)
    axc.tick_params(length=0, labelsize=8, colors=TEXT_SECONDARY)

    # ---------------- Panel D: best minPAE per protein, by oracle ----------------
    # This is the decisive panel: it shows whether a given oracle's PAE puts
    # non-binders anywhere near binders. A useful oracle must leave a gap.
    if len(oracles) > 1:
        axd = fig.add_subplot(gs[1, 1])
        axd.set_facecolor(SURFACE)
        ORACLE_COLOR = {"rf3": "#2a78d6", "protenix": "#eb6834"}
        ypos = {p: len(proteins) - 1 - i for i, p in enumerate(proteins)}
        for oracle in oracles:
            rows = {r["protein"]: r for r in summaries.get(oracle, [])}
            xs = [fnum(rows[p]["min_pae_min"]) for p in proteins if p in rows]
            ys = [ypos[p] for p in proteins if p in rows]
            axd.plot(xs, ys, marker="o", ls="", markersize=9,
                     color=ORACLE_COLOR.get(oracle, "#2a78d6"),
                     markeredgecolor=SURFACE, markeredgewidth=2.0,
                     label=oracle, zorder=4)
        # class bands behind the dots so the three arms are readable at a glance
        for p in proteins:
            axd.axhspan(ypos[p] - 0.5, ypos[p] + 0.5,
                        color=CLASS_COLOR[prot_class[p]], alpha=0.07, zorder=0)
        axd.set_yticks([ypos[p] for p in proteins])
        axd.set_yticklabels(proteins, fontsize=8.5)
        for lbl, p in zip(axd.get_yticklabels(), proteins):
            lbl.set_color(CLASS_COLOR[prot_class[p]])
        axd.set_xlabel("best minPAE achieved anywhere on the panel (Å)",
                       fontsize=9, color=TEXT_SECONDARY)
        axd.set_title("D   can the oracle tell a non-binder from a binder?",
                      fontsize=11.5, weight="bold", color=TEXT_PRIMARY, loc="left", pad=10)
        axd.grid(axis="x", color=GRID, lw=0.8, zorder=1)
        axd.set_axisbelow(True)
        for s in ("top", "right", "bottom", "left"):
            axd.spines[s].set_visible(False)
        axd.tick_params(length=0, labelsize=8, colors=TEXT_SECONDARY)
        handles = [plt.Line2D([], [], marker="o", ls="", markersize=9,
                              markeredgecolor=SURFACE, markeredgewidth=2.0,
                              color=ORACLE_COLOR[o], label=o) for o in oracles]
        handles += [plt.Line2D([], [], marker="s", ls="", markersize=7,
                               color=CLASS_COLOR[c], label=CLASS_LABEL[c])
                    for c in CLASS_ORDER if c in set(prot_class.values())]
        axd.legend(handles=handles, frameon=False, fontsize=7.8,
                   labelcolor=TEXT_SECONDARY, loc="upper right")

    fig.suptitle("Does ΔminPAE separate specific DNA binding from non-specific and non-binding?",
                 fontsize=13.5, weight="bold", color=TEXT_PRIMARY, x=0.115, ha="left", y=0.965)
    fig.text(0.115, 0.928,
             "8 natural controls × 8 DNA targets, folded MSA-free against a shared 24-bp panel, "
             "on two open AF3-class oracles.  Black ring = that TF's cognate site.",
             fontsize=9, color=TEXT_SECONDARY, ha="left")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out, dpi=200, facecolor=SURFACE)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
