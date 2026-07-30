#!/usr/bin/env python3
"""Prepare + submit the control-panel folds via pecli, under a hard spend cap.

`pecli prepare` is free and validates the input; `pecli submit` spends money. This
driver does prepare -> parse the run id -> submit, one fold at a time, and refuses
to start a fold that would take projected spend past --max-spend. It records every
run id into the manifest so results can be collected later.

Per-oracle flags are not interchangeable:
  rf3       --diffusion-batch-size 1 --seed 42     (one sample, reproducible)
  protenix  --need-atom-confidence true --sample 1 --seeds 42
            The PAE flag is REQUIRED: without it protenix writes only a summary
            confidence and no per-token PAE, so ΔminPAE cannot be computed.

Usage:
    python submit_control_folds.py --manifest folds/rf3/folds_manifest.json \
        --unit-cost 0.06 --max-spend 6.00
    python submit_control_folds.py --manifest folds/protenix/folds_manifest.json \
        --unit-cost 0.13 --max-spend 9.00 --dry-run
"""
from __future__ import annotations
import argparse
import json
import os
import re
import subprocess
import sys

RUN_ID_RE = re.compile(r"(?:Prepared|prepared)\s+\S+\s+(?:run|pipeline)\s+([0-9a-f]{6})")

ORACLE_FLAGS = {
    "rf3": ["--diffusion-batch-size", "1", "--seed", "42"],
    "protenix": ["--need-atom-confidence", "true", "--sample", "1", "--seeds", "42"],
}


def run(cmd, timeout=900):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--unit-cost", type=float, required=True,
                    help="realised $/run for this oracle (see `pecli cost`)")
    ap.add_argument("--max-spend", type=float, required=True,
                    help="hard cap; a fold that would exceed it is not submitted")
    ap.add_argument("--group", default="oracle-controls")
    ap.add_argument("--limit", type=int, default=None, help="only the first N pending folds")
    ap.add_argument("--dry-run", action="store_true", help="prepare only, never submit")
    args = ap.parse_args()

    mpath = os.path.abspath(args.manifest)
    manifest = json.load(open(mpath))
    fold_dir = os.path.dirname(mpath)
    oracle = manifest[0]["oracle"]
    flags = ORACLE_FLAGS[oracle]

    pending = [m for m in manifest if not m.get("run_id")]
    if args.limit:
        pending = pending[:args.limit]

    projected = len(pending) * args.unit_cost
    print(f"oracle={oracle}  pending={len(pending)}  unit=${args.unit_cost:.3f}  "
          f"projected=${projected:.2f}  cap=${args.max_spend:.2f}")
    if projected > args.max_spend:
        n_ok = int(args.max_spend // args.unit_cost)
        print(f"  projected spend exceeds cap -- will submit at most {n_ok} folds")
        pending = pending[:n_ok]

    spent = 0.0
    submitted, failed = 0, []
    for i, rec in enumerate(pending, 1):
        if spent + args.unit_cost > args.max_spend:
            print(f"  cap reached at ${spent:.2f}; stopping with "
                  f"{len(pending) - i + 1} folds unsubmitted")
            break
        inp = os.path.join(fold_dir, rec["fold_input"])
        rc, out = run(["pecli", "prepare", oracle, "--input", inp, *flags,
                       "--description", f"oracle-controls {rec['fold_id']} ({oracle}, MSA-free)"])
        m = RUN_ID_RE.search(out)
        if rc != 0 or not m:
            failed.append((rec["fold_id"], "prepare: " + out.strip()[-300:]))
            continue
        run_id = m.group(1)
        if "MSA pipeline" in out or "auto-routed" in out:
            # a pipeline means the input lacked a precomputed MSA -- that is a paid
            # extra step and breaks the MSA-free design, so refuse it.
            failed.append((rec["fold_id"],
                           f"prepare auto-routed to an MSA pipeline ({run_id}); "
                           "input is missing its precomputed MSA -- not submitting"))
            continue
        if args.dry_run:
            print(f"  [{i}/{len(pending)}] {rec['fold_id']}: prepared {run_id} (dry-run)")
            rec["prepared_id"] = run_id
            continue
        rc2, out2 = run(["pecli", "submit", run_id, "-y", "--group", args.group])
        if rc2 != 0:
            failed.append((rec["fold_id"], "submit: " + out2.strip()[-300:]))
            continue
        rec["run_id"] = run_id
        spent += args.unit_cost
        submitted += 1
        print(f"  [{i}/{len(pending)}] {rec['fold_id']}: submitted {run_id}  "
              f"(spent ≈ ${spent:.2f})")
        # persist after every submit so an interruption never loses run ids
        with open(mpath, "w") as f:
            json.dump(manifest, f, indent=2)

    with open(mpath, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\nsubmitted {submitted}  failed {len(failed)}  spend ≈ ${spent:.2f}")
    for fid, err in failed:
        print(f"  ! {fid}: {err}")
    still = len([m for m in manifest if not m.get("run_id")])
    print(f"manifest updated: {mpath}  ({still} folds still unsubmitted)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
