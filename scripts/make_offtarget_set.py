#!/usr/bin/env python3
"""Build the off-target DNA set for the specificity block.

The specificity block ranks a binder by ΔminPAE = min over off-targets of
minPAE(off) − minPAE(on). That needs an off-target panel. For the PRNP-site the
paper evaluates specificity two ways, and we build both:

  1. **Single-base-substitution variants** of the on-target site. The paper
     characterizes the PRNP binder DBS5 as "specific over 35/40 single-base
     variants", so the SBS panel is the fine-grained specificity test: every
     position × every alternative base (3 per position → 3·L variants for an
     L-bp site).
  2. **Unrelated decoy sites** — the other Table 1 targets, as gross off-targets
     a good binder should reject outright.

Each off-target is emitted as a duplex (both strands) so it can be folded with
the on-target-designed protein in the templated all-by-all fold.

Usage:
    python make_offtarget_set.py --on-target TGAGGAGAGGAG \
        --out specs/specificity_block/offtargets.json
"""
from __future__ import annotations
import argparse
import json

_COMP = str.maketrans("ACGT", "TGCA")
BASES = "ACGT"

# Corrected 2026-07-31 after a close read of the paper's Methods. The two panels
# below serve DIFFERENT purposes and must not be conflated:
#
#   ranking  on-target + the other Table 1 targets. This is what the paper's
#            ΔminPAE all-by-all actually folds against ("on-target + the six core
#            targets, + ten additional targets if the design's on-target was in
#            the additional set").
#   sbs      the single-base-variant sweep. The paper's "specific over 35/40
#            single-base variants" is WET-LAB characterisation of one binder
#            (DBS5) after ranking -- it is NOT the ΔminPAE ranking panel.
#
# The earlier version of this script put all 3*L single-base variants into one
# panel with the decoys, which (a) inflated the all-by-all ~4x in GPU cost and
# (b) silently changed the metric: ΔminPAE is a MINIMUM over off-targets, so
# including sequences one base from the on-target makes it a near-worst-case
# statistic rather than the discrimination-against-unrelated-sites statistic the
# paper reports.
#
# UNRESOLVED: which of the Table 1 targets are the paper's "six core" vs "ten
# additional" is not established from the accessible text, so `ranking` uses all
# other Table 1 targets we have (a superset of the core six).
PANEL_NOTE = (
    "panel=ranking reproduces the paper's ΔminPAE all-by-all (on-target + other "
    "Table 1 targets). panel=sbs is the single-base-variant sweep, which the paper "
    "used for wet-lab characterisation of an already-selected binder, NOT for "
    "ΔminPAE ranking. Do not rank on the sbs panel."
)


def revcomp(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]


# Other Table 1 targets (Sehgal et al. 2026, Table 1) used as unrelated decoys.
# Sequences transcribed from Table 1 of the Sehgal et al. 2026 paper. VERIFY each
# against the published Table 1 before a production run -- they were read from the
# paper text and have not been cross-checked against a second source.
TABLE1_DECOYS = {
    "Oct4gRNA2": "GGGCTTGCGA",
    "TBP": "CGTATAAACG",
    "CAG": "CAGCAGCAGCAG",
    "HSTelo": "AGGGTTAGGGTT",
    "NFkB": "GGGGATTCCCCC",
    "HD": "GCTTAATTAGCG",
    "P53": "AGACATGTCT",
    "Tbox": "AGGTGTGAAG",
    "FKH": "GCGTAAACAA",
}


def single_base_variants(seq: str):
    seq = seq.upper()
    out = []
    for i, wt in enumerate(seq):
        for b in BASES:
            if b == wt:
                continue
            var = seq[:i] + b + seq[i + 1:]
            out.append((f"sbs_{i+1}{wt}>{b}", var))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--on-target", required=True, help="on-target sense strand 5'->3'")
    ap.add_argument("--out", required=True)
    ap.add_argument("--panel", choices=["ranking", "sbs", "both"], default="ranking",
                    help="ranking (default): on-target + other Table 1 targets -- the panel "
                         "the paper's ΔminPAE all-by-all actually uses. sbs: on-target + every "
                         "single-base variant, for characterising ONE already-ranked design. "
                         "both: the union (NOT what the paper ranks on; see the note below).")
    args = ap.parse_args()

    on = args.on_target.upper()
    want_decoys = args.panel in ("ranking", "both")
    want_sbs = args.panel in ("sbs", "both")

    entries = []
    # on-target itself (reference point for ΔminPAE)
    entries.append({"id": "on_target", "kind": "on_target", "sense": on, "antisense": revcomp(on)})
    if want_decoys:
        # other Table 1 targets. Only equal-length ones can share a template; keep all, tag length.
        for name, seq in TABLE1_DECOYS.items():
            entries.append({
                "id": f"decoy_{name}", "kind": "decoy", "sense": seq,
                "antisense": revcomp(seq), "same_length_as_on": len(seq) == len(on),
            })
    if want_sbs:
        for name, var in single_base_variants(on):
            entries.append({"id": name, "kind": "sbs", "sense": var, "antisense": revcomp(var)})

    bundle = {
        "on_target": on,
        "length_bp": len(on),
        "panel": args.panel,
        "n_sbs": sum(1 for e in entries if e["kind"] == "sbs"),
        "n_decoys": sum(1 for e in entries if e["kind"] == "decoy"),
        "offtargets": entries,
        "_note": PANEL_NOTE,
    }
    with open(args.out, "w") as f:
        json.dump(bundle, f, indent=2)
    print(f"on-target {on} ({len(on)} bp)   panel={args.panel}")
    print(f"{bundle['n_decoys']} Table 1 off-targets, {bundle['n_sbs']} single-base variants")
    print(f"{len(entries)} total folds per design")
    if args.panel == "both":
        print("\nWARNING: --panel both is NOT the paper's ranking panel. ΔminPAE taken over "
              "single-base variants is a far harsher statistic than over unrelated sites, "
              "and it inflates the all-by-all ~4x. Use --panel ranking to rank.")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
