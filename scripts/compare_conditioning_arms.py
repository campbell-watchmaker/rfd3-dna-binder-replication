#!/usr/bin/env python3
"""Compare binder-block outcomes across H-bond conditioning settings.

The rfd3na spec conditions diffusion on a chosen subset of major-groove H-bond
donor/acceptor atoms (scripts/make_rfd3na_specs.py --hbond-central-bases). Too few and
the binder has no reason to sit in the major groove; too many and diffusion is
over-constrained. This prints the side-by-side needed to pick the setting.

METRICS, and why each is here
-----------------------------
protein-only Ca-RMSD   does the designed sequence FOLD as intended, ignoring where it
                       sits? Isolates foldability from placement.
DNA-aligned Ca-RMSD    the paper's self-consistency gate: superpose on the DNA, then
                       measure protein displacement. Conflates fold and placement, so
                       it is only interpretable next to the protein-only number.
major-groove H-bonds   does the binder READ BASES, or just grip the phosphate
                       backbone? This is the metric the conditioning setting is
                       supposed to move. Counted on the REFOLD (designed sidechains):
                       counting it on the pre-design backbone is meaningless, because
                       rfd3na emits placeholder sidechains (heavily Asn/Ala/Thr/Gln)
                       that are not the designed sequence.
ipTM                   the oracle's own interface confidence.

A CAUTION ON READING THIS
-------------------------
Each arm is an independent draw of a handful of backbones, and rfd3na exposes no seed
control, so the arms are not paired. Only a LARGE difference is resolvable; a shifted
median with overlapping ranges is noise. The script prints n and the full ranges so
that is visible rather than hidden behind a median.

Usage:
    python compare_conditioning_arms.py \
        --arm "6 atoms=results/binder_block/all_designs.csv" \
        --arm "8 atoms=results/binder_block/all_designs_hb4.csv"
"""
from __future__ import annotations
import argparse
import csv
import statistics as st
import sys


def load(path):
    rows = list(csv.DictReader(open(path)))
    out = []
    for r in rows:
        try:
            out.append({
                "dna_rmsd": float(r["dna_aligned_ca_rmsd"]) if r["dna_aligned_ca_rmsd"] else None,
                "prot_rmsd": float(r["protein_only_ca_rmsd"]) if r.get("protein_only_ca_rmsd") else None,
                "iptm": float(r["iptm"]) if r["iptm"] else None,
                "hb": int(r["protein_dna_hbonds"]) if r["protein_dna_hbonds"] else None,
                "mg": int(r["major_groove_hbonds"]) if r["major_groove_hbonds"] else None,
            })
        except (TypeError, ValueError):
            continue
    return out


def summarise(vals):
    v = [x for x in vals if x is not None]
    if not v:
        return "        n/a"
    return f"{min(v):6.2f} {st.median(v):7.2f} {max(v):7.2f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", action="append", required=True, metavar="LABEL=CSV")
    ap.add_argument("--rmsd-gate", type=float, default=8.0,
                    help="the paper's pre-resample gate (default 8.0 A)")
    args = ap.parse_args()

    arms = []
    for spec in args.arm:
        label, _, path = spec.partition("=")
        arms.append((label, load(path)))

    print(f"{'metric':26s} " + " ".join(f"{'   min  median     max':>22s}" for _ in arms))
    print(f"{'':26s} " + " ".join(f"{l:>22s}" for l, _ in arms))
    print("-" * (26 + 23 * len(arms)))
    for key, name in (("prot_rmsd", "protein-only Ca-RMSD"),
                      ("dna_rmsd", "DNA-aligned Ca-RMSD"),
                      ("mg", "major-groove H-bonds"),
                      ("hb", "total protein-DNA H-b"),
                      ("iptm", "ipTM")):
        print(f"{name:26s} " + " ".join(summarise([r[key] for r in rows]) for _, rows in arms))

    print()
    for label, rows in arms:
        n = len(rows)
        passers = sum(1 for r in rows if r["dna_rmsd"] is not None
                      and r["dna_rmsd"] < args.rmsd_gate)
        folds = sum(1 for r in rows if r["prot_rmsd"] is not None and r["prot_rmsd"] < 3.0)
        reads = sum(1 for r in rows if r["mg"] is not None and r["mg"] >= 3)
        print(f"{label}: n={n}")
        print(f"   pass {args.rmsd_gate:.0f} A DNA-aligned gate : {passers}/{n}")
        print(f"   folds correctly (<3 A prot)  : {folds}/{n}")
        print(f"   reads bases (>=3 major-groove): {reads}/{n}")
    print()
    print("Independent draws, no shared seed -- treat a shifted median with overlapping")
    print("ranges as noise, not signal.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
