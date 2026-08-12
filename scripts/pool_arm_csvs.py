#!/usr/bin/env python3
"""Concatenate per-arm filter CSVs so a factorial's MAIN EFFECTS can be compared.

Phase 2.5 is a 2x2: CFG (on/off) x H-bond conditioning (fixed/sampled), 20 designs per
cell. `compare_conditioning_arms.py` compares whatever labelled CSVs it is handed, so the
four cells go in directly -- but a main effect needs the two cells sharing a level pooled
into one CSV first (n=40 rather than n=20). That is all this does.

Rows are tagged with their source arm in an `arm` column, so a pooled CSV stays
decomposable and nothing is silently merged beyond recovery.

Why bother pooling at all: at n=20 per cell only a large effect is resolvable, and rfd3na
exposes no seed control so the cells are independent draws rather than paired. Pooling to
n=40 per level is the whole reason the factorial was run crossed instead of as two
separate A/B pairs.

Usage:
    python pool_arm_csvs.py --out pooled_cfg_on.csv \
        fixed_cfgon=<csv> sampled_cfgon=<csv>
"""
from __future__ import annotations
import argparse
import csv
import os
import sys


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("inputs", nargs="+", metavar="LABEL=CSV")
    args = ap.parse_args()

    rows, fields = [], []
    for spec in args.inputs:
        if "=" not in spec:
            ap.error(f"expected LABEL=CSV, got {spec!r}")
        label, path = spec.split("=", 1)
        if not os.path.isfile(path):
            print(f"missing {path}", file=sys.stderr)
            return 1
        with open(path) as f:
            r = csv.DictReader(f)
            for name in r.fieldnames or []:
                if name not in fields:
                    fields.append(name)
            n = 0
            for row in r:
                row["arm"] = label
                rows.append(row)
                n += 1
        print(f"  {label}: {n} row(s) from {os.path.basename(path)}")

    if "arm" not in fields:
        fields.append("arm")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} row(s) -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
