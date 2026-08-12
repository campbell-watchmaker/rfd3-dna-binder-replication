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

PADDING TO A COMMON LENGTH (added 2026-08-06, on by default)
------------------------------------------------------------
4 of the 9 Table 1 decoys are 10 bp against a 12-bp on-target. minPAE is a MINIMUM
over protein x DNA token pairs, so a shorter duplex simply offers fewer pairs to
minimise over and is systematically disadvantaged as an off-target -- which biases
ΔminPAE upward for reasons that have nothing to do with specificity. Every target is
therefore centred in a verified-neutral flank at one fixed length, reusing
`build_duplex()` and `verify_panel()` from analysis/oracle_controls/control_panel.py.
This is the same confound the oracle control panel was padded to 24 bp to remove, and
padding there did not degrade discrimination (rf3 argmin 4/5 on padded duplexes).

`verify_panel()` additionally checks that no padded decoy has picked up another
target's motif in its flank -- contamination there would compress ΔminPAE for the
wrong reason. A previous `same_length_as_on` field recorded the problem and was never
read by anything.

This is a DEVIATION from the paper, which appears to fold Table 1 sites at native
length. Recorded as such in docs/replication_log.md. `--no-pad` restores native
lengths.

Usage:
    python make_offtarget_set.py --on-target TGAGGAGAGGAG \
        --out specs/specificity_block/offtargets.json
"""
from __future__ import annotations
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "analysis", "oracle_controls"))

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


def pad_entries(entries, fixed_bp, on_target_id="on_target"):
    """Centre every entry's motif in the neutral flank at `fixed_bp`, verified clean.

    Mutates each entry to carry the padded duplex and keeps the bare motif alongside,
    so a downstream consumer can still report which site was folded. Returns
    (entries, problems) -- `problems` is verify_panel()'s report; a non-empty list means
    a padded duplex carries a motif it should not and the panel must not be used.
    """
    from control_panel import build_duplex, verify_panel

    panel, motifs = {}, {}
    for e in entries:
        sense, anti, lpad, rpad = build_duplex(e["sense"], fixed_bp=fixed_bp)
        e["motif"] = e["sense"]
        e["motif_len"] = len(e["motif"])
        e["sense"], e["antisense"] = sense, anti
        e["left_pad"], e["right_pad"] = lpad, rpad
        e["padded_to_bp"] = fixed_bp
        panel[e["id"]] = {"sense": sense}
        motifs[e["id"]] = e["motif"]
    return entries, verify_panel(panel, motifs=motifs)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--on-target", required=True, help="on-target sense strand 5'->3'")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pad-to-bp", type=int, default=None,
                    help="pad every duplex to this length (default: the longest target "
                         "in the panel, so nothing is truncated and padding is minimal)")
    ap.add_argument("--no-pad", action="store_true",
                    help="emit native lengths. minPAE then favours the on-target purely "
                         "because shorter decoys offer fewer token pairs -- use only to "
                         "reproduce the unpadded comparison.")
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
        for name, seq in TABLE1_DECOYS.items():
            entries.append({
                "id": f"decoy_{name}", "kind": "decoy", "sense": seq,
                "antisense": revcomp(seq),
            })
    if want_sbs:
        for name, var in single_base_variants(on):
            entries.append({"id": name, "kind": "sbs", "sense": var, "antisense": revcomp(var)})

    fixed_bp, problems = None, []
    if not args.no_pad:
        fixed_bp = args.pad_to_bp or max(len(e["sense"]) for e in entries)
        entries, problems = pad_entries(entries, fixed_bp)

    bundle = {
        "on_target": on,
        "length_bp": len(on),
        "padded_to_bp": fixed_bp,
        "panel": args.panel,
        "n_sbs": sum(1 for e in entries if e["kind"] == "sbs"),
        "n_decoys": sum(1 for e in entries if e["kind"] == "decoy"),
        "offtargets": entries,
        "_note": PANEL_NOTE,
    }
    print(f"on-target {on} ({len(on)} bp)   panel={args.panel}")
    print(f"{bundle['n_decoys']} Table 1 off-targets, {bundle['n_sbs']} single-base variants")
    print(f"{len(entries)} total folds per design")
    if fixed_bp:
        native = sorted({e["motif_len"] for e in entries})
        print(f"padded every duplex to {fixed_bp} bp (native lengths present: {native})")
        print("  -> removes the token-count bias that favours the on-target when decoys "
              "are shorter")
    else:
        print("NOT padded: shorter decoys offer fewer protein x DNA token pairs, so "
              "minPAE favours the on-target for reasons unrelated to specificity")
    if problems:
        print(f"\nREFUSING TO WRITE -- {len(problems)} padded duplex(es) carry a motif they "
              "should not; ΔminPAE would be compressed for the wrong reason:")
        for p in problems:
            print("  ! " + p)
        return 1

    with open(args.out, "w") as f:
        json.dump(bundle, f, indent=2)
    if args.panel == "both":
        print("\nWARNING: --panel both is NOT the paper's ranking panel. ΔminPAE taken over "
              "single-base variants is a far harsher statistic than over unrelated sites, "
              "and it inflates the all-by-all ~4x. Use --panel ranking to rank.")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
