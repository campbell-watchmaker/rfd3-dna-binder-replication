#!/usr/bin/env python3
"""Compute rfd3na conditioning inputs (ori tokens + major-groove H-bond atoms)
from a folded B-DNA duplex.

Run this on the DNA duplex CIF produced by the GPU-side fold (protenix/openfold3),
e.g.:

    python compute_conditioning.py \
        --duplex ../results/prnp_duplex/prnp_duplex_model.cif \
        --out ../targets/prnp/prnp_conditioning.json \
        --window 6 --offset 3.0

Emits a JSON bundle with the ori-token placements and the candidate major-groove
donor/acceptor atoms, plus a chemical validation block (purine-N7-on-major-side
fraction and mean helical twist) so a reviewer can confirm the geometry is sane.
"""
from __future__ import annotations
import argparse, json, sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from target_prep import (
    load_duplex, base_pairs, ori_tokens, hbond_candidates, validate_major_groove,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--duplex", required=True, help="folded B-DNA duplex (.cif/.pdb)")
    ap.add_argument("--out", required=True, help="output conditioning JSON")
    ap.add_argument("--window", type=int, default=6, help="bp per ori token")
    ap.add_argument("--offset", type=float, default=3.0, help="Å into major groove")
    args = ap.parse_args()

    dna = load_duplex(args.duplex)
    pairs = base_pairs(dna)
    frac_n7, twist = validate_major_groove(dna, pairs)
    if frac_n7 < 0.9:
        print(f"WARNING: only {frac_n7:.0%} of purine N7 on major side — "
              "check duplex geometry / strand pairing.", file=sys.stderr)

    bundle = {
        "duplex_file": os.path.basename(args.duplex),
        "n_base_pairs": len(pairs),
        "params": {"window_bp": args.window, "major_groove_offset_A": args.offset},
        "validation": {
            "purine_N7_on_major_frac": round(frac_n7, 3),
            "mean_interior_twist_deg": round(twist, 1),
        },
        "ori_tokens": ori_tokens(dna, window=args.window, offset=args.offset),
        "hbond_candidates": hbond_candidates(dna),
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(bundle, f, indent=2)
    print(f"wrote {args.out}: {len(bundle['ori_tokens'])} ori tokens, "
          f"{len(bundle['hbond_candidates'])} H-bond candidates; "
          f"N7-major {frac_n7:.0%}, twist {twist:.1f}°")


if __name__ == "__main__":
    main()
