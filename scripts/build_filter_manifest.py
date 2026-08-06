#!/usr/bin/env python3
"""Build the filter_binder_block.py manifest from a completed refold batch.

The binder-block self-consistency metric is a DNA-aligned protein Cα-RMSD: superpose
the refold onto the *design* using the DNA atoms only, then measure how far the
protein moved. Low RMSD means the designed sequence actually folds back into the
backbone it was designed for, in the same place relative to the duplex.

That needs three things paired up per design, which live in three different places:

  design_path   the RELAXED rfd3na backbone -- the structure LigandMPNN designed onto,
                not the raw diffusion output. (Using the unrelaxed one would measure
                the relax as if it were refold error.)
  refold_path   the rf3 prediction of the designed sequence + the same duplex
  iptm          from the rf3 summary confidence, for the ipTM gate

Usage:
    python build_filter_manifest.py \
        --refold-manifest folds/refold/folds_manifest.json \
        --relaxed-dir     relaxed \
        --raw-dir         raw/stage6 \
        --out             filter_manifest.json
"""
from __future__ import annotations
import argparse
import glob
import json
import os
import sys


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refold-manifest", required=True)
    ap.add_argument("--relaxed-dir", required=True,
                    help="directory of relaxed rfd3na backbones (<backbone>.pdb)")
    ap.add_argument("--raw-dir", required=True,
                    help="downloaded rf3 refold results, one subdir per fold_id")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    recs = json.load(open(args.refold_manifest))
    rows, missing = [], []
    for r in recs:
        design = os.path.join(args.relaxed_dir, r["backbone"] + ".pdb")
        # rf3 writes <name>_model.cif; prefer the top-level (best-ranked) copy
        cands = sorted(
            glob.glob(os.path.join(args.raw_dir, r["fold_id"], "**", "*_model.cif"),
                      recursive=True),
            key=lambda p: (p.count(os.sep), len(p)))
        if not os.path.isfile(design):
            missing.append(f"{r['fold_id']}: no relaxed design at {design}")
            continue
        if not cands:
            missing.append(f"{r['fold_id']}: no refold CIF under {args.raw_dir}")
            continue

        iptm = None
        summ = sorted(glob.glob(os.path.join(args.raw_dir, r["fold_id"], "**",
                                             "*_summary_confidences.json"),
                                recursive=True),
                      key=lambda p: (p.count(os.sep), len(p)))
        if summ:
            try:
                iptm = json.load(open(summ[0])).get("iptm")
            except (OSError, ValueError):
                pass

        # The per-token PAE, for the specificity block's on-target minPAE entry gate.
        # `*_summary_confidences.json` holds only scalars; the matrix is in the sibling
        # `*_confidences.json`, so the glob must EXCLUDE the summary file or it matches
        # both and may pick the one with no `pae` key.
        pae_path = None
        conf = [p for p in sorted(glob.glob(os.path.join(args.raw_dir, r["fold_id"], "**",
                                                         "*_confidences.json"),
                                            recursive=True),
                                  key=lambda p: (p.count(os.sep), len(p)))
                if not p.endswith("_summary_confidences.json")]
        if conf:
            pae_path = conf[0]

        rows.append({
            "design_id": r["fold_id"], "oracle": "rf3",
            "design_path": design, "refold_path": cands[0],
            "iptm": iptm, "pae_path": pae_path,
            "protein_chain": r.get("protein_chain", "A"),
            "dna_chains": r.get("dna_chains", ["B", "C"]),
            "protein_len": r.get("protein_len"),
            "backbone": r["backbone"], "seq_id": r["seq_id"],
            "overall_confidence": r.get("overall_confidence"),
            "ligand_confidence": r.get("ligand_confidence"),
        })

    with open(args.out, "w") as f:
        json.dump(rows, f, indent=2)
    print(f"{len(rows)} design/refold pairs -> {args.out}")
    n_iptm = sum(1 for r in rows if r["iptm"] is not None)
    n_pae = sum(1 for r in rows if r["pae_path"] is not None)
    print(f"  ipTM recovered for {n_iptm}/{len(rows)}")
    print(f"  PAE matrix found for {n_pae}/{len(rows)} (needed for the specificity-block "
          "minPAE entry gate)")
    if missing:
        print(f"  MISSING ({len(missing)}):")
        for m in missing[:10]:
            print("    ! " + m)
    return 0


if __name__ == "__main__":
    sys.exit(main())
