#!/usr/bin/env python3
"""Validate the specificity block's ranking code on the natural-TF control panel.

WHY THIS EXISTS
---------------
Before the specificity block is trusted on de novo designs -- whose correct answer
nobody knows -- run its actual ranking code on a panel whose answer is already
measured. The control panel folded 8 proteins x 8 DNA targets on rf3 and
`analysis/oracle_controls/RESULTS.md` records the outcome: argmin lands on the
protein's own cognate site for **4 of 5** sequence-specific TFs, and ΔminPAE > 0 for
4 of 5.

That result was produced by `compute_control_metrics.py`, which is control-panel
specific. The specificity block ranks with `scripts/compute_delta_minpae.py` instead.
Two implementations of the same statistic is exactly where a silent inconsistency
lives, so this script reshapes the control panel's manifest into the specificity
block's manifest schema, runs the SPECIFICITY BLOCK'S OWN script over it, and checks
the answer against the recorded one.

If it reproduces, the block's plumbing is sound and any later disagreement on designs
is about the designs. If it does not, the bug is ours.

Costs nothing: it reads PAE files already on disk. Re-download them free with
`collect_control_results.py` if a scratchpad has been cleaned.

The shape translation, which is the substance of this script:

    control panel                    specificity block
    -----------------------------    ---------------------------------
    protein   (Zif268, ...)      ->  design_id
    dna_id    (zif268_site, ...) ->  dna_id
    is_on_target True/False      ->  kind: "on_target" / "decoy"
    protein_chains[0]            ->  protein_chain
    protein_len x protein_copies ->  protein_len (a residue COUNT; see
                                     _protein_len_range in the consumer)

Each of the 5 specific TFs becomes one "design" with 1 on-target and 7 off-targets --
the other four TFs' cognate sites plus scramble/polygc/prnp. The 3 non-specific
controls have no cognate site and are therefore unrankable by construction; they are
excluded, and the consumer is expected to report any that slip through rather than
drop them silently.

Usage:
    python validate_specificity_block.py \
        --control-manifest <rf3_manifest.json, pae_path filled> \
        --out-dir <scratch/validate>
"""
from __future__ import annotations
import argparse
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..")
sys.path.insert(0, HERE)

# The recorded rf3 answer this validation must reproduce
# (analysis/oracle_controls/RESULTS.md, 64 folds, 2026-07-30).
EXPECTED_ARGMIN_CORRECT = 4
EXPECTED_DELTA_POSITIVE = 4
EXPECTED_N_SPECIFIC = 5


def to_specificity_manifest(control_manifest, oracle="rf3"):
    """Reshape control-panel records into the specificity block's manifest schema."""
    from control_panel import ON_TARGET

    out, skipped = [], []
    for rec in control_manifest:
        if rec.get("oracle") != oracle:
            continue
        protein = rec["protein"]
        if protein not in ON_TARGET:
            # no cognate site, so no ΔminPAE reference point -- unrankable by
            # construction, not by omission
            skipped.append((protein, rec["dna_id"], "no on-target defined"))
            continue
        if not rec.get("pae_path") or not os.path.exists(rec["pae_path"]):
            skipped.append((protein, rec["dna_id"], "no PAE on disk"))
            continue
        out.append({
            "design_id": protein,
            "dna_id": rec["dna_id"],
            "kind": "on_target" if rec["is_on_target"] else "decoy",
            "oracle": oracle,
            "pae_path": rec["pae_path"],
            "protein_chain": rec["protein_chains"][0],
            "dna_chains": rec["dna_chains"],
            "protein_len": rec["protein_len"],
            "protein_copies": rec.get("protein_copies", 1),
        })
    return out, skipped


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--control-manifest", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--oracle", default="rf3")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    control = json.load(open(args.control_manifest))
    jobs, skipped = to_specificity_manifest(control, args.oracle)
    if not jobs:
        print(f"no usable {args.oracle} folds with a PAE on disk. Re-download free with "
              "collect_control_results.py (WITHOUT --skip-status, which refuses to fetch "
              "anything not already present).", file=sys.stderr)
        return 1

    mpath = os.path.join(args.out_dir, "specificity_shaped_manifest.json")
    json.dump(jobs, open(mpath, "w"), indent=2)
    n_designs = len({j["design_id"] for j in jobs})
    print(f"reshaped {len(jobs)} control folds -> {n_designs} 'designs' "
          f"({len(skipped)} folds excluded: no cognate site / no PAE)")

    out_csv = os.path.join(args.out_dir, "delta_minpae_controls.csv")
    per_complex = os.path.join(args.out_dir, "min_pae_controls.csv")
    r = subprocess.run([sys.executable,
                        os.path.join(REPO, "scripts", "compute_delta_minpae.py"),
                        "--manifest", mpath, "--out", out_csv,
                        "--per-complex-out", per_complex],
                       capture_output=True, text=True)
    print("\n--- scripts/compute_delta_minpae.py ---")
    print(r.stdout.rstrip())
    if r.returncode != 0:
        print(r.stderr[-2000:], file=sys.stderr)
        return 1

    rows = list(csv.DictReader(open(out_csv)))
    per = list(csv.DictReader(open(per_complex)))

    # argmin: does the panel-wide lowest minPAE for a protein land on its own site?
    from control_panel import ON_TARGET
    by_design = {}
    for p in per:
        by_design.setdefault(p["design_id"], []).append(p)
    argmin_ok, argmin_detail = 0, []
    for design, recs in sorted(by_design.items()):
        best = min(recs, key=lambda p: float(p["min_pae"]))
        ok = best["dna_id"] == ON_TARGET[design]
        argmin_ok += ok
        argmin_detail.append((design, best["dna_id"], ON_TARGET[design], ok,
                              float(best["min_pae"])))

    n_pos = sum(1 for r_ in rows if float(r_["delta_min_pae"]) > 0)

    print(f"\n{'protein':11s} {'ΔminPAE':>9s} {'on':>7s} {'best off':>9s}  "
          f"argmin -> (expected)")
    print("-" * 68)
    delta_by = {r_["design_id"]: r_ for r_ in rows}
    for design, got, want, ok, mp in argmin_detail:
        d = delta_by.get(design, {})
        print(f"{design:11s} {d.get('delta_min_pae', '?'):>9s} "
              f"{d.get('on_target_min_pae', '?'):>7s} "
              f"{d.get('best_offtarget_min_pae', '?'):>9s}  "
              f"{'OK ' if ok else 'MISS'} {got} -> ({want})")

    print(f"\nargmin on own site : {argmin_ok}/{len(argmin_detail)}  "
          f"(recorded: {EXPECTED_ARGMIN_CORRECT}/{EXPECTED_N_SPECIFIC})")
    print(f"ΔminPAE > 0        : {n_pos}/{len(rows)}  "
          f"(recorded: {EXPECTED_DELTA_POSITIVE}/{EXPECTED_N_SPECIFIC})")

    problems = []
    if len(rows) != EXPECTED_N_SPECIFIC:
        problems.append(f"ranked {len(rows)} designs, expected {EXPECTED_N_SPECIFIC} "
                        "specific TFs -- folds are missing from the panel")
    if argmin_ok != EXPECTED_ARGMIN_CORRECT:
        problems.append(f"argmin {argmin_ok}/{len(argmin_detail)} != recorded "
                        f"{EXPECTED_ARGMIN_CORRECT}/{EXPECTED_N_SPECIFIC}")
    if n_pos != EXPECTED_DELTA_POSITIVE:
        problems.append(f"ΔminPAE>0 count {n_pos} != recorded {EXPECTED_DELTA_POSITIVE}")

    if problems:
        print("\nDISAGREES with the recorded control-panel result:")
        for p in problems:
            print("  ! " + p)
        print("\nThe two paths compute the same statistic, so a disagreement is a bug in "
              "one of them -- not a property of the designs. Do not run the specificity "
              "block on designs until this reconciles.")
        return 1

    print("\nPASS: the specificity block's own ranking code reproduces the control "
          "panel's recorded answer. Its plumbing is sound; later disagreement on designs "
          "is about the designs.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
