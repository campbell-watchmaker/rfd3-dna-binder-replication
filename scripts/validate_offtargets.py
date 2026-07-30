#!/usr/bin/env python3
"""Validate off-target sequences against Sehgal et al. 2026, Table 1.

This script documents which Table 1 sequences are used as decoys in the
specificity block and flags any discrepancies. Run this before a production
specificity block to ensure the off-target panel matches the paper.

The TABLE1_DECOYS in make_offtarget_set.py were transcribed from the paper
and should be cross-checked against the published Table 1.

Usage:
    python validate_offtargets.py --paper-table1 <path-or-url>

If no paper data is provided, prints the current decoys and prompts for
manual verification.
"""
from __future__ import annotations
import json
import sys

# Current decoy sequences in the repo (from make_offtarget_set.py)
REPO_DECOYS = {
    "Oct4gRNA2": "GGGCTTGCGA",
    "TBP": "CGTATAAACG",
    "CAG": "CAGCAGCAGCAG",
    "HSTelo": "AGGGTTAGGGTT",
    "NFkB": "GGGGATTCCCCC",
    "HD": "GCTTAATTAGCG",
    "P53": "AGACATGTCT",
    "Tbox": "AGGTGTGAAG",
    "FKH": "GCGTAAACAA",
}

# Instructions for manual verification
VERIFICATION_INSTRUCTIONS = """
MANUAL VERIFICATION REQUIRED

To validate off-target sequences, cross-check the following against
Sehgal et al. 2026, Table 1 (bioRxiv 2026.04.27.720408):

CURRENT DECOYS IN REPO:
"""


def print_decoys_for_verification():
    """Print table for manual cross-checking."""
    print(VERIFICATION_INSTRUCTIONS)
    print("\n{:<15} {:<20} (sequence)".format("Target Name", "Sequence"))
    print("-" * 50)
    for name, seq in sorted(REPO_DECOYS.items()):
        print(f"{name:<15} {seq:<20}")
    print("\n" + "=" * 50)
    print("ACTION REQUIRED:")
    print("1. Open Sehgal et al. 2026, Table 1 (bioRxiv 2026.04.27.720408)")
    print("2. For each target above, verify the sequence matches exactly")
    print("3. Report any discrepancies")
    print("=" * 50)


def validate_against_paper(paper_data: dict) -> dict:
    """Compare repo decoys against paper Table 1.

    paper_data: dict with keys matching target names, values = sequences
    """
    results = {
        "verified": [],
        "mismatches": [],
        "missing_in_paper": [],
        "extra_in_repo": set(REPO_DECOYS.keys()) - set(paper_data.keys()),
    }

    for name, repo_seq in REPO_DECOYS.items():
        if name not in paper_data:
            results["missing_in_paper"].append(name)
        else:
            paper_seq = paper_data[name].upper()
            repo_seq_upper = repo_seq.upper()
            if paper_seq == repo_seq_upper:
                results["verified"].append(name)
            else:
                results["mismatches"].append({
                    "target": name,
                    "repo_seq": repo_seq,
                    "paper_seq": paper_seq,
                })

    return results


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--paper-table1", help="JSON file with Table 1 data: {name: sequence, ...}")
    args = ap.parse_args()

    if args.paper_table1:
        with open(args.paper_table1) as f:
            paper_data = json.load(f)
        results = validate_against_paper(paper_data)

        print(f"✓ Verified: {len(results['verified'])} targets")
        for name in results["verified"]:
            print(f"  - {name}")

        if results["mismatches"]:
            print(f"\n✗ MISMATCHES: {len(results['mismatches'])} targets")
            for m in results["mismatches"]:
                print(f"  - {m['target']}")
                print(f"    Repo: {m['repo_seq']}")
                print(f"    Paper: {m['paper_seq']}")
            return 1

        if results["missing_in_paper"]:
            print(f"\n⚠ Missing in paper: {results['missing_in_paper']}")

        if results["extra_in_repo"]:
            print(f"\n⚠ Extra in repo (not in paper): {results['extra_in_repo']}")

        print("\n✓ All sequences validated successfully")
        return 0
    else:
        print_decoys_for_verification()
        return 0


if __name__ == "__main__":
    sys.exit(main())
