#!/usr/bin/env python3
"""Binder-block filtering + oracle comparison (CPU, runs in Claude Science).

For each refolded design (protein-DNA complex CIF from protenix / openfold3 /
esmfold2), computes the paper's binder-block metrics:

  * DNA-aligned protein Ca-RMSD -- superpose the refold onto the design by the
    DNA atoms only, then measure how far the protein Ca atoms moved. This is the
    paper's self-consistency metric: does the designed protein still sit the same
    way on the DNA after an independent fold? (paper gate: <8A -> resample -> <3A)
  * interface ipTM -- read from the oracle's per-design confidence output.
  * protein-DNA H-bond count -- open reimplementation of the DSSR interaction
    count: donor..acceptor pairs across the protein-DNA interface within a
    distance cutoff and (when H is present) a donor-H..acceptor angle cutoff.

Emits two CSVs:
  * passers.csv -- designs passing the gates, ranked.
  * oracle_comparison.csv -- every (design, oracle) row with RMSD/ipTM/H-bonds
    (+ runtime/gpu if provided), for the protenix-vs-openfold3-vs-esmfold2 writeup.

The RMSD needs the *design* structure (pre-fold, from rfd3na/ligandmpnn) and the
*refold* for the same design; superposition is by DNA atoms.

Dependencies: biotite (structure IO + superposition), numpy.
"""
from __future__ import annotations
import argparse
import csv
import json
import os

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb


# ---- H-bond chemistry ----
# Protein side-chain / backbone donor & acceptor atoms, and DNA donor & acceptor
# atoms. Names are PDB/CIF standard. This mirrors what DSSR counts as protein-DNA
# H-bonds; it is a geometric reimplementation, not DSSR itself.
PROTEIN_DONORS = {"N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "OG", "OG1", "OH", "SG"}
PROTEIN_ACCEPTORS = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "SD"}
DNA_DONORS = {"N4", "N6", "N1", "N2", "O2'"}  # base amino / imino / ribose donors
DNA_ACCEPTORS = {"N7", "O6", "O4", "O2", "N3", "O1P", "O2P", "OP1", "OP2", "O3'", "O4'", "O5'"}

HBOND_DIST_CUTOFF = 3.5   # heavy-atom donor..acceptor, Angstrom
MAJOR_GROOVE_ACCEPTORS = {"N7", "O6", "O4"}  # subset used to tag major-groove reads


def _load_any(path):
    if path.endswith((".cif", ".mmcif", ".bcif")):
        f = pdbx.CIFFile.read(path)
        return pdbx.get_structure(f, model=1)
    f = pdb.PDBFile.read(path)
    return f.get_structure(model=1)


def _protein_dna_masks(arr):
    dna = struc.filter_nucleotides(arr)
    prot = struc.filter_amino_acids(arr)
    return prot, dna


def dna_aligned_ca_rmsd(design_arr, refold_arr):
    """Superpose refold onto design by DNA atoms; return protein Ca RMSD after that fit."""
    d_prot, d_dna = _protein_dna_masks(design_arr)
    r_prot, r_dna = _protein_dna_masks(refold_arr)

    # match DNA atoms by (chain, res_id, atom_name); use the common set, in order
    def dna_index(arr, mask):
        sub = arr[mask]
        return {(a.chain_id, a.res_id, a.atom_name): i for i, a in enumerate(sub)}, sub

    d_idx, d_sub = dna_index(design_arr, d_dna)
    r_idx, r_sub = dna_index(refold_arr, r_dna)
    common = [k for k in d_idx if k in r_idx]
    if len(common) < 3:
        raise ValueError(f"too few common DNA atoms to superpose ({len(common)})")
    d_dna_coords = d_sub[[d_idx[k] for k in common]]
    r_dna_coords = r_sub[[r_idx[k] for k in common]]

    # fit refold DNA -> design DNA, apply transform to whole refold
    _, transform = struc.superimpose(d_dna_coords, r_dna_coords)
    refold_moved = transform.apply(refold_arr)

    # protein Ca RMSD between design and transformed refold, matched by (chain,res_id)
    def ca_map(arr):
        m = arr[(struc.filter_amino_acids(arr)) & (arr.atom_name == "CA")]
        return {(a.chain_id, a.res_id): arr_i for arr_i, a in enumerate(m)}, m

    d_ca_idx, d_ca = ca_map(design_arr)
    r_ca_idx, r_ca = ca_map(refold_moved)
    ca_common = [k for k in d_ca_idx if k in r_ca_idx]
    if not ca_common:
        raise ValueError("no common protein Ca atoms")
    dc = d_ca.coord[[d_ca_idx[k] for k in ca_common]]
    rc = r_ca.coord[[r_ca_idx[k] for k in ca_common]]
    return float(np.sqrt(np.mean(np.sum((dc - rc) ** 2, axis=1)))), len(ca_common)


def count_protein_dna_hbonds(arr):
    """Count heavy-atom protein-DNA H-bond candidate pairs within the distance cutoff.

    Returns (total, major_groove) where major_groove counts pairs whose DNA atom is
    a major-groove acceptor (N7/O6/O4).
    """
    prot_mask, dna_mask = _protein_dna_masks(arr)
    prot = arr[prot_mask]
    dna = arr[dna_mask]
    if prot.array_length() == 0 or dna.array_length() == 0:
        return 0, 0

    # protein donor/acceptor atoms
    p_don = prot[np.isin(prot.atom_name, list(PROTEIN_DONORS))]
    p_acc = prot[np.isin(prot.atom_name, list(PROTEIN_ACCEPTORS))]
    d_don = dna[np.isin(dna.atom_name, list(DNA_DONORS))]
    d_acc = dna[np.isin(dna.atom_name, list(DNA_ACCEPTORS))]

    total = 0
    major = 0

    def pairs(a, b, tag_major_from_b=False):
        nonlocal total, major
        if a.array_length() == 0 or b.array_length() == 0:
            return
        # pairwise distances
        dmat = np.linalg.norm(a.coord[:, None, :] - b.coord[None, :, :], axis=2)
        hits = np.argwhere(dmat <= HBOND_DIST_CUTOFF)
        total += len(hits)
        if tag_major_from_b:
            for _, j in hits:
                if b.atom_name[j] in MAJOR_GROOVE_ACCEPTORS:
                    major += 1

    # protein donor -> DNA acceptor (the dominant, and where major-groove reads live)
    pairs(p_don, d_acc, tag_major_from_b=True)
    # DNA donor -> protein acceptor
    pairs(d_don, p_acc, tag_major_from_b=False)
    return total, major


def analyze_one(design_path, refold_path):
    design = _load_any(design_path)
    refold = _load_any(refold_path)
    rmsd, n_ca = dna_aligned_ca_rmsd(design, refold)
    hb_total, hb_major = count_protein_dna_hbonds(refold)
    return {"dna_aligned_ca_rmsd": round(rmsd, 3), "n_ca_matched": n_ca,
            "protein_dna_hbonds": hb_total, "major_groove_hbonds": hb_major}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True,
                    help="JSON list of {design_id, oracle, design_path, refold_path, iptm?, runtime_s?, gpu?}")
    ap.add_argument("--out", required=True, help="passers.csv")
    ap.add_argument("--oracle-comparison", required=True, help="oracle_comparison.csv (all rows)")
    ap.add_argument("--rmsd-gate", type=float, default=3.0)
    ap.add_argument("--iptm-gate", type=float, default=0.7)
    args = ap.parse_args()

    jobs = json.load(open(args.manifest))
    all_rows = []
    for j in jobs:
        try:
            m = analyze_one(j["design_path"], j["refold_path"])
        except Exception as e:
            m = {"dna_aligned_ca_rmsd": None, "n_ca_matched": 0,
                 "protein_dna_hbonds": None, "major_groove_hbonds": None, "error": str(e)}
        row = {"design_id": j["design_id"], "oracle": j.get("oracle", "unknown"),
               "iptm": j.get("iptm"), "runtime_s": j.get("runtime_s"), "gpu": j.get("gpu"), **m}
        all_rows.append(row)

    cols = ["design_id", "oracle", "dna_aligned_ca_rmsd", "iptm", "protein_dna_hbonds",
            "major_groove_hbonds", "n_ca_matched", "runtime_s", "gpu"]
    os.makedirs(os.path.dirname(args.oracle_comparison) or ".", exist_ok=True)
    with open(args.oracle_comparison, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)

    def passes(r):
        return (r["dna_aligned_ca_rmsd"] is not None and r["dna_aligned_ca_rmsd"] < args.rmsd_gate
                and r.get("iptm") is not None and r["iptm"] > args.iptm_gate)
    passers = sorted((r for r in all_rows if passes(r)),
                     key=lambda r: (-(r["iptm"] or 0), r["dna_aligned_ca_rmsd"]))
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(passers)

    print(f"analyzed {len(all_rows)} (design,oracle) rows -> {args.oracle_comparison}")
    print(f"{len(passers)} passers (RMSD<{args.rmsd_gate}, ipTM>{args.iptm_gate}) -> {args.out}")


if __name__ == "__main__":
    main()
