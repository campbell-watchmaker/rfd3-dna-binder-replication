#!/usr/bin/env python3
"""Build per-design complex-JSON fold inputs from a LigandMPNN FASTA + the target DNA.

Each designed protein is folded WITH the target DNA duplex so the oracle predicts
the protein-DNA complex (and yields interface ipTM). Both DNA strands are always
written explicitly -- esmfold2 does not auto-generate the complementary strand,
and protenix/openfold3 accept the same shape, so one file serves all three oracles.

Usage:
    python build_fold_inputs.py --fasta designs.fasta --dna TGAGGAGAGGAG \
        --out-dir specs/binder_block/fold_inputs [--skip-wt]
"""
from __future__ import annotations
import argparse
import json
import os

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]


def read_fasta(path):
    recs = []
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--dna", required=True, help="target DNA sense strand 5'->3'")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--skip-wt", action="store_true",
                    help="skip the first FASTA record (LigandMPNN emits the WT input as record 0)")
    args = ap.parse_args()

    dna_sense = args.dna.upper()
    dna_anti = revcomp(dna_sense)
    os.makedirs(args.out_dir, exist_ok=True)

    recs = read_fasta(args.fasta)
    if args.skip_wt and recs:
        recs = recs[1:]

    n = 0
    for k, (name, seq) in enumerate(recs):
        cj = {
            "id": name,
            "chains": [
                {"id": "A", "type": "protein", "sequence": seq},
                {"id": "B", "type": "dna", "sequence": dna_sense},
                {"id": "C", "type": "dna", "sequence": dna_anti},
            ],
        }
        with open(os.path.join(args.out_dir, f"{name}.json"), "w") as f:
            json.dump(cj, f, indent=2)
        n += 1
    print(f"wrote {n} complex-JSON fold inputs to {args.out_dir}")
    print(f"DNA sense 5'->3': {dna_sense}  |  antisense: {dna_anti}")


if __name__ == "__main__":
    main()
