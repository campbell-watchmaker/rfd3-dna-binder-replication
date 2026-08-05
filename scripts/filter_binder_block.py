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


def _dna_strand_sequences(arr):
    """{chain_id: one-letter base sequence} for each DNA chain, 5'->3' by res_id."""
    _, dna_mask = _protein_dna_masks(arr)
    dna = arr[dna_mask]
    out = {}
    for ch in sorted(set(dna.chain_id.tolist())):
        c = dna[dna.chain_id == ch]
        out[ch] = "".join(c.res_name[c.res_id == r][0][-1] for r in np.unique(c.res_id))
    return out


def _match_dna_chains(design_arr, refold_arr):
    """Map refold DNA chain id -> design DNA chain id, BY BASE SEQUENCE.

    Chain letters are NOT comparable between these two structures. rfd3na writes the
    duplex first and the designed protein last (DNA = A,B; protein = C), while the
    refold spec declares the protein first (protein = A; DNA = B,C). Matching DNA on
    the chain letter therefore pairs the design's ANTISENSE strand with the refold's
    SENSE strand -- and because the sugar-phosphate backbone atom names are identical
    in every nucleotide, ~130 atoms "match" and the superposition silently succeeds on
    the wrong strand, returning a confident, meaningless RMSD.

    So pair the strands by sequence instead, which is unambiguous for a
    non-palindromic duplex.
    """
    d_seqs = _dna_strand_sequences(design_arr)
    r_seqs = _dna_strand_sequences(refold_arr)
    mapping, used = {}, set()
    for r_ch, r_seq in r_seqs.items():
        for d_ch, d_seq in d_seqs.items():
            if d_ch not in used and d_seq == r_seq:
                mapping[r_ch] = d_ch
                used.add(d_ch)
                break
    if len(mapping) != len(r_seqs):
        raise ValueError(
            f"could not pair DNA strands by sequence: design {d_seqs} vs refold {r_seqs}")
    return mapping


def dna_aligned_ca_rmsd(design_arr, refold_arr):
    """Superpose refold onto design by DNA atoms; return protein Ca RMSD after that fit.

    Both correspondences are established WITHOUT trusting chain letters (see
    _match_dna_chains): DNA strands are paired by base sequence, and protein Ca atoms
    are matched in sequential order along the single designed chain. The refold is a
    prediction of the same sequence, so residue i corresponds to residue i.
    """
    chain_map = _match_dna_chains(design_arr, refold_arr)

    def dna_index(arr, remap=None):
        _, mask = _protein_dna_masks(arr)
        sub = arr[mask]
        idx = {}
        for i, a in enumerate(sub):
            ch = remap.get(a.chain_id, a.chain_id) if remap else a.chain_id
            idx[(ch, a.res_id, a.atom_name)] = i
        return idx, sub

    d_idx, d_sub = dna_index(design_arr)
    r_idx, r_sub = dna_index(refold_arr, remap=chain_map)
    common = sorted(k for k in d_idx if k in r_idx)
    if len(common) < 3:
        raise ValueError(f"too few common DNA atoms to superpose ({len(common)})")
    d_dna_coords = d_sub[[d_idx[k] for k in common]]
    r_dna_coords = r_sub[[r_idx[k] for k in common]]

    # fit refold DNA -> design DNA, apply transform to whole refold
    _, transform = struc.superimpose(d_dna_coords, r_dna_coords)
    refold_moved = transform.apply(refold_arr)

    # protein Ca, matched in sequential order (chain letters differ; see above)
    def ca_ordered(arr):
        m = arr[struc.filter_amino_acids(arr) & (arr.atom_name == "CA")]
        order = np.lexsort((m.res_id, m.chain_id))
        return m[order]

    d_ca = ca_ordered(design_arr)
    r_ca = ca_ordered(refold_moved)
    n = min(d_ca.array_length(), r_ca.array_length())
    if n == 0:
        raise ValueError("no protein Ca atoms")
    if d_ca.array_length() != r_ca.array_length():
        # a length mismatch means these are not the same design; refuse rather than
        # silently compare a truncated prefix
        raise ValueError(
            f"protein length mismatch: design has {d_ca.array_length()} Ca, "
            f"refold has {r_ca.array_length()}")
    dc, rc = d_ca.coord, r_ca.coord
    return float(np.sqrt(np.mean(np.sum((dc - rc) ** 2, axis=1)))), n


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
    ap.add_argument("--stage", choices=["pre_resample", "post_resample"],
                    default="post_resample",
                    help="Which binder-block checkpoint this is. The paper gates TWICE: "
                         "pre_resample keeps DNA-aligned RMSD < 8 A to decide what goes into "
                         "LigandMPNN resampling (no ipTM gate at that point); post_resample "
                         "applies RMSD < 3 A + ipTM > 0.7 to the resampled folds.")
    ap.add_argument("--rmsd-gate", type=float, default=None,
                    help="override the stage default (8.0 pre_resample, 3.0 post_resample)")
    ap.add_argument("--iptm-gate", type=float, default=0.7,
                    help="ignored at --stage pre_resample")
    args = ap.parse_args()

    # Paper sequence (Methods, "Binder block"): fold -> RMSD < 8 A -> LigandMPNN
    # resample -> fold -> RMSD < 3 A, ipTM > 0.7, high H-bond counts. The 8 A gate
    # was missing here, which would have sent every diffused backbone into
    # resampling instead of only the self-consistent ones.
    if args.rmsd_gate is None:
        args.rmsd_gate = 8.0 if args.stage == "pre_resample" else 3.0
    use_iptm = args.stage == "post_resample"

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
        if r["dna_aligned_ca_rmsd"] is None or r["dna_aligned_ca_rmsd"] >= args.rmsd_gate:
            return False
        if not use_iptm:
            return True
        return r.get("iptm") is not None and r["iptm"] > args.iptm_gate
    passers = sorted((r for r in all_rows if passes(r)),
                     key=lambda r: (-(r["iptm"] or 0), r["dna_aligned_ca_rmsd"]))
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(passers)

    print(f"analyzed {len(all_rows)} (design,oracle) rows -> {args.oracle_comparison}")
    gate = f"RMSD<{args.rmsd_gate}" + (f", ipTM>{args.iptm_gate}" if use_iptm else " (no ipTM gate)")
    print(f"stage={args.stage}: {len(passers)} passers ({gate}) -> {args.out}")
    if args.stage == "pre_resample":
        print("  -> feed these to LigandMPNN resampling, then re-run with "
              "--stage post_resample on the resampled folds")


if __name__ == "__main__":
    main()
