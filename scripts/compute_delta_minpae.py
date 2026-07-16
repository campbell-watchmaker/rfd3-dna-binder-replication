#!/usr/bin/env python3
"""Compute minPAE and ΔminPAE for the specificity block (CPU, runs in Claude Science).

The specificity block folds each on-target-designed protein against the on-target
DNA and every off-target DNA (templated all-by-all), then ranks designs by how
much better they read the on-target than the best off-target:

    minPAE(complex) = min over protein-residue i, DNA-residue j of PAE(i, j)
    ΔminPAE(design) = min over off-targets of minPAE(off) − minPAE(on)

A large positive ΔminPAE means every off-target is predicted with a WORSE
(higher) best protein-DNA PAE than the on-target -- i.e. the binder is confidently
paired only with its intended site. This is the paper's specificity metric.

PAE is a per-residue-pair matrix produced by AF3-class folders (protenix /
openfold3). **esmfold2 does not emit a PAE matrix, so the specificity block uses
protenix / openfold3 only** -- unlike the binder block's three-oracle comparison.

Input: a manifest listing, per (design, dna_target) fold, the path to the fold's
PAE output and enough chain/length info to slice the protein-vs-DNA submatrix.
The PAE file is read as:
  * a .json/.npz/.npy with a square PAE array (protenix/openfold3 confidence
    output), OR
  * a .json with {"pae": [[...]], "token_chain_ids": [...]} (chain labels per token).

If token chain ids are present we use them to pick protein vs DNA tokens; otherwise
the manifest must give protein_len and dna token ranges explicitly.

Usage:
    python compute_delta_minpae.py --manifest folds_manifest.json \
        --out results/specificity_block/delta_minpae.csv
"""
from __future__ import annotations
import argparse
import csv
import json
import os

import numpy as np


def load_pae(path):
    """Return (pae 2D array, token_chain_ids or None)."""
    if path.endswith(".npy"):
        return np.load(path), None
    if path.endswith(".npz"):
        z = np.load(path)
        key = "pae" if "pae" in z else z.files[0]
        return z[key], None
    # json
    d = json.load(open(path))
    if isinstance(d, list):
        return np.asarray(d, dtype=float), None
    pae = None
    for k in ("pae", "predicted_aligned_error", "pae_matrix"):
        if k in d:
            pae = np.asarray(d[k], dtype=float)
            break
    if pae is None:
        raise ValueError(f"no PAE array found in {path}")
    chains = d.get("token_chain_ids") or d.get("chain_ids") or d.get("asym_id")
    return pae, (list(chains) if chains is not None else None)


def protein_dna_token_masks(n_tokens, chains, protein_chain, dna_chains, protein_len, dna_ranges):
    """Boolean masks over tokens for protein vs DNA. Prefer chain labels; else use explicit ranges."""
    if chains is not None:
        chains = np.asarray(chains)
        prot = chains == protein_chain
        dna = np.isin(chains, list(dna_chains))
        return prot, dna
    # fallback: explicit index ranges (0-based, [lo,hi) ) from the manifest
    prot = np.zeros(n_tokens, bool)
    dna = np.zeros(n_tokens, bool)
    plo, phi = protein_len
    prot[plo:phi] = True
    for lo, hi in dna_ranges:
        dna[lo:hi] = True
    return prot, dna


def min_pae(pae, prot_mask, dna_mask):
    """min over protein-DNA residue pairs of PAE(i,j), using both PAE orientations."""
    block1 = pae[np.ix_(prot_mask, dna_mask)]
    block2 = pae[np.ix_(dna_mask, prot_mask)]
    vals = []
    if block1.size:
        vals.append(block1.min())
    if block2.size:
        vals.append(block2.min())
    if not vals:
        raise ValueError("empty protein-DNA PAE block")
    return float(min(vals))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True,
                    help="JSON list of {design_id, dna_id, kind(on_target|sbs|decoy), pae_path, "
                         "oracle, protein_chain?, dna_chains?, protein_len?, dna_ranges?}")
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-complex-out", default=None, help="optional CSV of every minPAE")
    args = ap.parse_args()

    jobs = json.load(open(args.manifest))

    # minPAE per (design, dna, oracle)
    rows = []
    for j in jobs:
        pae, chains = load_pae(j["pae_path"])
        prot_mask, dna_mask = protein_dna_token_masks(
            pae.shape[0], chains,
            j.get("protein_chain", "A"), j.get("dna_chains", ["B", "C"]),
            tuple(j["protein_len"]) if j.get("protein_len") else (0, 0),
            [tuple(r) for r in j.get("dna_ranges", [])],
        )
        mp = min_pae(pae, prot_mask, dna_mask)
        rows.append({"design_id": j["design_id"], "dna_id": j["dna_id"],
                     "kind": j["kind"], "oracle": j.get("oracle", "unknown"), "min_pae": round(mp, 4)})

    if args.per_complex_out:
        with open(args.per_complex_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["design_id", "dna_id", "kind", "oracle", "min_pae"])
            w.writeheader(); w.writerows(rows)

    # ΔminPAE per (design, oracle)
    from collections import defaultdict
    by_design = defaultdict(dict)   # (design,oracle) -> {dna_id: (kind, minpae)}
    for r in rows:
        by_design[(r["design_id"], r["oracle"])][r["dna_id"]] = (r["kind"], r["min_pae"])

    out_rows = []
    for (design, oracle), d in by_design.items():
        on = [v for v in d.values() if v[0] == "on_target"]
        offs = [v[1] for v in d.values() if v[0] in ("sbs", "decoy")]
        if not on or not offs:
            continue
        on_mp = on[0][1]
        best_off = min(offs)                     # most competitive off-target
        delta = best_off - on_mp                 # >0 = on-target preferred
        # also report separated sbs / decoy worst cases
        sbs = [v[1] for v in d.values() if v[0] == "sbs"]
        decoy = [v[1] for v in d.values() if v[0] == "decoy"]
        out_rows.append({
            "design_id": design, "oracle": oracle,
            "on_target_min_pae": round(on_mp, 4),
            "best_offtarget_min_pae": round(best_off, 4),
            "delta_min_pae": round(delta, 4),
            "min_sbs_min_pae": round(min(sbs), 4) if sbs else None,
            "min_decoy_min_pae": round(min(decoy), 4) if decoy else None,
            "n_offtargets": len(offs),
        })

    out_rows.sort(key=lambda r: -r["delta_min_pae"])   # rank by specificity, descending
    cols = ["design_id", "oracle", "delta_min_pae", "on_target_min_pae",
            "best_offtarget_min_pae", "min_sbs_min_pae", "min_decoy_min_pae", "n_offtargets"]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader(); w.writerows(out_rows)

    print(f"{len(rows)} (design,dna,oracle) complexes -> minPAE")
    print(f"{len(out_rows)} (design,oracle) ranked by ΔminPAE -> {args.out}")
    if out_rows:
        top = out_rows[0]
        print(f"top: {top['design_id']} ({top['oracle']}) ΔminPAE={top['delta_min_pae']}")


if __name__ == "__main__":
    main()
