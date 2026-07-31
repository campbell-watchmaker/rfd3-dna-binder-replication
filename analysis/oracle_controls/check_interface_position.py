#!/usr/bin/env python3
"""Where on the duplex does minPAE actually land?

THE HOLE THIS CLOSES
--------------------
minPAE is a minimum over ALL protein-token x DNA-token pairs. It never checks
*where* on the duplex the best-paired position sits. If a protein's minimum falls
on flanking DNA rather than inside the motif window, then minPAE is not measuring
motif recognition -- and because every target in the panel shares the SAME neutral
flank, flank-reading would also compress differences between targets and deflate
ΔminPAE toward zero.

So: for every fold, find the DNA token attaining minPAE, map it to a position on
its strand, and ask whether that position is inside the motif.

BASELINE
--------
A motif of length m centred in a duplex of length L occupies m/L of the positions,
so a positionally-random argmin lands in-motif m/L of the time. That fraction is
reported per fold, and the aggregate in-motif rate is compared against the
mean baseline -- an in-motif rate at or below baseline means the metric is
picking up no positional signal.

Motif window, per strand
------------------------
sense:     [len(lpad), len(lpad) + m)
antisense: antisense = revcomp(sense), so antisense position i pairs with sense
           position L-1-i. The motif therefore occupies
           [L - (len(lpad) + m), L - len(lpad))

Usage:
    python check_interface_position.py --manifest folds/rf3/folds_manifest.json \
        --out-csv results/oracle_controls/rf3_interface_positions.csv
"""
from __future__ import annotations
import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
from control_panel import build_panel, FIXED_BP  # noqa: E402
from compute_delta_minpae import _norm_chain  # noqa: E402


def strand_token_indices(rec, d):
    """Return (protein_mask, [(strand_name, [token indices 5'->3'])]) for one fold."""
    with open(rec["pae_path"]) as f:
        blob = json.load(f)

    if rec["oracle"] == "rf3":
        pae = np.asarray(blob["pae"], dtype=float)
        chains = np.asarray([_norm_chain(c) for c in blob["token_chain_ids"]])
        want_p = {_norm_chain(c) for c in rec["protein_chains"]}
        prot = np.isin(chains, list(want_p))
        strands = []
        for name, cid in zip(("sense", "antisense"), rec["dna_chains"]):
            idx = np.where(chains == _norm_chain(cid))[0]
            strands.append((name, idx))
        return pae, prot, strands

    # protenix: positional, asym ids ascending = entity order
    key = next(k for k in ("token_pair_pae", "pae") if k in blob)
    pae = np.asarray(blob[key], dtype=float)
    asym = np.asarray(blob["token_asym_id"])
    order = sorted(set(asym.tolist()))
    nc = rec["protein_copies"]
    prot = np.isin(asym, order[:nc])
    strands = []
    for name, a in zip(("sense", "antisense"), order[nc:nc + 2]):
        strands.append((name, np.where(asym == a)[0]))
    return pae, prot, strands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()

    manifest = [m for m in json.load(open(args.manifest))
                if m.get("pae_path") and os.path.exists(m["pae_path"])]
    panel = build_panel()

    rows, errors = [], []
    for rec in manifest:
        d = panel[rec["dna_id"]]
        m, L = d["motif_len"], len(d["sense"])
        lp = len(d["left_pad"])
        windows = {
            "sense": (lp, lp + m),
            "antisense": (L - (lp + m), L - lp),
        }
        try:
            pae, prot, strands = strand_token_indices(rec, d)
            best = None
            for name, idx in strands:
                if len(idx) != L:
                    raise ValueError(f"strand {name} has {len(idx)} tokens, expected {L}")
                # min over protein tokens for each position on this strand
                sub = np.minimum(pae[np.ix_(prot, idx)].min(axis=0),
                                 pae[np.ix_(idx, prot)].min(axis=1))
                pos = int(np.argmin(sub))
                val = float(sub[pos])
                if best is None or val < best[2]:
                    best = (name, pos, val)
            name, pos, val = best
            lo, hi = windows[name]
            rows.append({
                "fold_id": rec["fold_id"], "protein": rec["protein"],
                "klass": rec["klass"], "dna_id": rec["dna_id"],
                "is_on_target": rec["is_on_target"], "oracle": rec["oracle"],
                "min_pae": round(val, 4), "strand": name, "position": pos,
                "motif_start": lo, "motif_end": hi - 1, "motif_len": m,
                "in_motif": lo <= pos < hi,
                "baseline_frac": round(m / L, 4),
            })
        except Exception as e:
            errors.append(f"{rec['fold_id']}: {e}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)
    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    def rate(rs):
        return (sum(r["in_motif"] for r in rs) / len(rs)) if rs else float("nan")

    def base(rs):
        return (sum(r["baseline_frac"] for r in rs) / len(rs)) if rs else float("nan")

    oracle = rows[0]["oracle"]
    print(f"=== {oracle}: where does minPAE land? ({len(rows)} folds, {len(errors)} errors) ===")
    for e in errors:
        print("  ! " + e)
    print(f"\nALL folds          in-motif {rate(rows):5.1%}   "
          f"positional baseline {base(rows):5.1%}")
    for k in ("specific", "nonspecific_binder", "nonbinder"):
        rs = [r for r in rows if r["klass"] == k]
        if rs:
            print(f"  {k:19s} in-motif {rate(rs):5.1%}   baseline {base(rs):5.1%}")

    ons = [r for r in rows if r["is_on_target"]]
    print(f"\nON-TARGET folds only (the critical set, n={len(ons)}):")
    hdr = f"  {'protein':11s} {'dna':12s} {'minPAE':>7s} {'strand':>10s} {'pos':>4s} {'motif':>9s} {'in?':>4s}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in sorted(ons, key=lambda r: r["min_pae"]):
        print(f"  {r['protein']:11s} {r['dna_id']:12s} {r['min_pae']:7.2f} "
              f"{r['strand']:>10s} {r['position']:4d} "
              f"{str(r['motif_start']) + '-' + str(r['motif_end']):>9s} "
              f"{'YES' if r['in_motif'] else 'no':>4s}")
    print(f"\n  on-target in-motif rate: {rate(ons):.1%} (baseline {base(ons):.1%})")
    print(f"\nwrote {args.out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
