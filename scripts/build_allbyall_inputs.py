#!/usr/bin/env python3
"""Build templated all-by-all fold inputs for the specificity block.

For each specificity-resampled design sequence, emit one complex-JSON fold input
per DNA target in the off-target set (on-target + single-base variants + decoys).
Each fold pairs the SAME protein sequence with a DIFFERENT DNA duplex, so folding
them all and comparing minPAE reveals which sites the protein reads confidently.

Also emits a `folds_manifest.json` skeleton listing every (design, dna) pair with
its `kind`, ready for compute_delta_minpae.py once the PAE paths are filled in
after the folds return.

Usage:
    python build_allbyall_inputs.py \
        --design-fasta survivors.fasta \
        --offtargets specs/specificity_block/offtargets.json \
        --out-dir specs/specificity_block/allbyall_inputs \
        [--skip-wt]
"""
from __future__ import annotations
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_fold_inputs import read_fasta  # reuse the FASTA reader


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design-fasta", required=True)
    ap.add_argument("--offtargets", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--skip-wt", action="store_true")
    args = ap.parse_args()

    off = json.load(open(args.offtargets))
    targets = off["offtargets"]
    designs = read_fasta(args.design_fasta)
    if args.skip_wt and designs:
        designs = designs[1:]

    os.makedirs(args.out_dir, exist_ok=True)
    manifest = []
    n = 0
    for dname, seq in designs:
        for t in targets:
            fold_id = f"{dname}__{t['id']}"
            cj = {
                "id": fold_id,
                "chains": [
                    {"id": "A", "type": "protein", "sequence": seq},
                    {"id": "B", "type": "dna", "sequence": t["sense"]},
                    {"id": "C", "type": "dna", "sequence": t["antisense"]},
                ],
            }
            with open(os.path.join(args.out_dir, f"{fold_id}.json"), "w") as f:
                json.dump(cj, f, indent=2)
            manifest.append({
                "design_id": dname, "dna_id": t["id"], "kind": t["kind"],
                "fold_input": f"{fold_id}.json",
                "pae_path": "FILL_AFTER_FOLD", "oracle": "protenix",
                "protein_chain": "A", "dna_chains": ["B", "C"],
            })
            n += 1

    with open(os.path.join(args.out_dir, "folds_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"{len(designs)} designs x {len(targets)} targets = {n} fold inputs -> {args.out_dir}")
    print(f"folds_manifest.json written ({len(manifest)} entries; fill pae_path after folds return)")


if __name__ == "__main__":
    main()
