#!/usr/bin/env python3
"""Extract predicted protein-only templates for the templated rf3 arm.

WHY A PREDICTED TEMPLATE RATHER THAN THE CRYSTAL CHAIN
-----------------------------------------------------
The first attempt templated from the deposited PDB chain. That failed a
token-count assertion: 1AAY chain A models only 85 of Zif268's 90 residues
(disordered termini are absent from crystals), so the templated fold carried a
different molecule from the sequence-based one and the two were not comparable.

A *predicted* structure has no unmodelled residues -- every residue in the input
sequence is placed -- so the length mismatch disappears by construction.

It is also the faithful analogue of the pipeline. Sehgal et al. template the
specificity block's all-by-all with "the most recent AF3 prediction before the
all-by-all folding", i.e. the protein chain from the design's own prior
prediction, never an experimental structure. Using a prediction here therefore
(a) fixes the length mismatch, (b) matches what the pipeline will actually do to
designs, and (c) removes the objection that crystal coordinates leak
experimental information into a benchmark.

WHICH PREDICTION IS THE TEMPLATE
--------------------------------
The pipeline's template is the protein chain from the fold that immediately
preceded the all-by-all -- a protein+on-target-DNA prediction. So:

  specific TFs        protein chain from that TF's own ON-TARGET untemplated fold
  Sac7d / non-binders no on-target exists, so the neutral `scramble` fold is used
                      as a fixed reference

The template is CONSTANT across all 8 DNA targets for a given protein, so it
cannot bias which target wins the argmin -- which is the property that makes this
safe for the specificity comparison. Only the DNA differs between the 8 folds.

Note the protein conformation is on-target-conditioned, exactly as in the paper.
That is an asymmetry the paper accepts; it is shared equally by every off-target
fold, so it shifts absolute minPAE without favouring any one target.

Only the protein chains are written. DNA and ligands are dropped, so no
protein-DNA docking geometry reaches the template -- verified by asserting zero
nucleotide residues in the output.

TWO MODES (the design mode added 2026-08-06)
--------------------------------------------
    --controls-panel   the 5+3 natural-TF control panel. Labels and expected residue
                       counts come from curated_controls.json; the source fold is the
                       protein's on-target (or neutral `scramble`) row.
    --design-fasta     the real pipeline: one template per DESIGN, expected residue count
                       from the FASTA sequence itself and expect_chains=1 (a design is a
                       single chain, unlike the panel's 2-copy dimers). The source fold is
                       that design's own Stage-6 refold -- which is exactly the paper's
                       "most recent AF3 prediction before the all-by-all folding".

Both modes share `extract()` verbatim, so the design path keeps the two gates that were
paid for in the control arm: zero surviving nucleotide residues, and an exact
residue-count match against the input sequence. Either one failing is what stops a
DNA-carrying or truncated template silently poisoning the all-by-all.

Usage:
    # control panel
    python make_predicted_templates.py --controls-panel \
        --raw-dir <raw/rf3> --out-dir <templates_predicted>

    # designs, for the specificity block
    python make_predicted_templates.py --design-fasta survivors.fasta \
        --raw-dir <raw/stage6> --out-dir <templates_predicted>
"""
from __future__ import annotations
import argparse
import glob
import json
import os
import sys

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from control_panel import load_controls, ON_TARGET  # noqa: E402

NEUTRAL_REFERENCE = "scramble"  # template source for proteins with no cognate site


def source_fold(label: str) -> str:
    """Which existing untemplated fold supplies this protein's template."""
    return f"{label}__{ON_TARGET.get(label, NEUTRAL_REFERENCE)}"


def extract(cif_path: str, out_path: str, expect_residues: int, expect_chains: int):
    arr = pdbx.get_structure(pdbx.CIFFile.read(cif_path), model=1)
    prot = arr[struc.filter_amino_acids(arr)]
    if prot.array_length() == 0:
        raise ValueError(f"{cif_path}: no amino-acid atoms found")

    n_res = struc.get_residue_count(prot)
    chains = sorted(set(prot.chain_id.tolist()))

    # A template CIF containing DNA would hand the model the docking geometry and
    # make the experiment circular, so this is a hard gate, not a warning.
    leftover_nuc = struc.get_residue_count(prot[struc.filter_nucleotides(prot)]) \
        if prot.array_length() else 0
    if leftover_nuc:
        raise ValueError(f"{out_path}: {leftover_nuc} nucleotide residues survived the filter")

    if n_res != expect_residues:
        raise ValueError(
            f"{out_path}: extracted {n_res} residues, expected {expect_residues} "
            "(protein_len x copies). A predicted structure should place every input "
            "residue -- do not proceed with a mismatched template.")
    if len(chains) != expect_chains:
        raise ValueError(
            f"{out_path}: got protein chains {chains}, expected {expect_chains} chain(s)")

    f = pdbx.CIFFile()
    pdbx.set_structure(f, prot)
    f.write(out_path)
    return n_res, chains


def _control_jobs(controls_path):
    """[(label, source_fold, expect_residues, expect_chains, is_on_target)] for the panel."""
    here = os.path.dirname(os.path.abspath(__file__))
    controls = load_controls(controls_path or os.path.join(here, "curated_controls.json"))
    return [(c["label"], source_fold(c["label"]),
             c["protein_length"] * c["copies"], c["copies"],
             c["label"] in ON_TARGET)
            for c in controls]


def _design_jobs(fasta_path):
    """Same tuple shape, for designs read from a survivors FASTA.

    The source fold is the design's own Stage-6 refold directory, named by the FASTA
    record -- which is the fold_id build_fold_inputs.py emitted. expect_chains is 1: a
    design is a single protein chain, so anything else means the refold CIF is not what
    this label claims.
    """
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "scripts"))
    from build_fold_inputs import read_fasta
    return [(name, name, len(seq), 1, True) for name, seq in read_fasta(fasta_path)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", required=True,
                    help="downloaded untemplated rf3 results, one subdir per fold")
    ap.add_argument("--out-dir", required=True)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--controls-panel", action="store_true",
                      help="the natural-TF control panel (labels from curated_controls.json)")
    mode.add_argument("--design-fasta",
                      help="one template per design; residue count from the sequence, "
                           "expect_chains=1")
    ap.add_argument("--controls", default=None, help="--controls-panel only")
    args = ap.parse_args()

    jobs = _design_jobs(args.design_fasta) if args.design_fasta \
        else _control_jobs(args.controls)
    if not jobs:
        print("no templates to build (empty controls list / FASTA)", file=sys.stderr)
        return 1
    os.makedirs(args.out_dir, exist_ok=True)

    report, problems = [], []
    for label, fold, exp_res, exp_chains, on_target in jobs:
        hits = sorted(glob.glob(os.path.join(args.raw_dir, fold, "**", "*_model.cif"),
                                recursive=True),
                      key=lambda p: (p.count(os.sep), len(p)))
        if not hits:
            problems.append(f"{label}: no predicted CIF under {args.raw_dir}/{fold}")
            continue
        out = os.path.join(args.out_dir, f"{label}_template.cif")
        try:
            n_res, chains = extract(hits[0], out, expect_residues=exp_res,
                                    expect_chains=exp_chains)
        except Exception as e:
            problems.append(f"{label}: {e}")
            continue
        report.append({
            "label": label, "template_cif": out, "source_fold": fold,
            "source_is_on_target": on_target,
            "chain_ids": chains, "n_residues": n_res,
            "expected_residues": exp_res,
            "protein_copies": exp_chains, "contains_nucleotide": False,
        })

    with open(os.path.join(args.out_dir, "templates_manifest.json"), "w") as f:
        json.dump(report, f, indent=2)

    print(f"predicted templates written: {len(report)}/{len(jobs)} -> {args.out_dir}\n")
    w = max([11] + [len(r["label"]) for r in report])
    print(f"  {'protein':{w}s} {'res':>5s} {'exp':>5s} {'chains':>8s}  source fold")
    print("  " + "-" * (w + 51))
    for r in report:
        print(f"  {r['label']:{w}s} {r['n_residues']:5d} {r['expected_residues']:5d} "
              f"{','.join(r['chain_ids']):>8s}  {r['source_fold']}"
              f"{'' if r['source_is_on_target'] else '  (neutral ref)'}")
    if problems:
        print(f"\nPROBLEMS ({len(problems)}):")
        for p in problems:
            print("  ! " + p)
        return 1
    print("\nall templates length-matched to their input sequence; no nucleotides present")
    return 0


if __name__ == "__main__":
    sys.exit(main())
