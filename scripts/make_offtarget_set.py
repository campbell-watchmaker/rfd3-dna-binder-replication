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
    ap = argparse.ArgumentParser()
    ap.add_argument("--on-target", required=True, help="on-target sense strand 5'->3'")
    ap.add_argument("--out", required=True)
    ap.add_argument("--include-decoys", action="store_true", default=True)
    args = ap.parse_args()

    on = args.on_target.upper()
    entries = []
    # on-target itself (reference point for ΔminPAE)
    entries.append({"id": "on_target", "kind": "on_target", "sense": on, "antisense": revcomp(on)})
    # single-base substitutions
    for name, var in single_base_variants(on):
        entries.append({"id": name, "kind": "sbs", "sense": var, "antisense": revcomp(var)})
    # unrelated decoys (only those of equal length can share a template; keep all, tag length)
    if args.include_decoys:
        for name, seq in TABLE1_DECOYS.items():
            entries.append({
                "id": f"decoy_{name}", "kind": "decoy", "sense": seq,
                "antisense": revcomp(seq), "same_length_as_on": len(seq) == len(on),
            })

    bundle = {
        "on_target": on,
        "length_bp": len(on),
        "n_sbs": sum(1 for e in entries if e["kind"] == "sbs"),
        "n_decoys": sum(1 for e in entries if e["kind"] == "decoy"),
        "offtargets": entries,
    }
    with open(args.out, "w") as f:
        json.dump(bundle, f, indent=2)
    print(f"on-target {on} ({len(on)} bp)")
    print(f"{bundle['n_sbs']} single-base variants (3 x {len(on)} = {3*len(on)} expected)")
    print(f"{bundle['n_decoys']} unrelated decoys")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
