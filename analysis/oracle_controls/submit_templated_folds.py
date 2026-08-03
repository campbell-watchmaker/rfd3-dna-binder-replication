#!/usr/bin/env python3
"""Prepare + stage + submit the TEMPLATED rf3 arm, under a hard spend cap.

A templated rf3 fold cannot be prepared directly. pecli decides "bare fold vs paid
msa->fold pipeline" at PREPARE time and its only signal is an in-input MSA carrier
hung on a *sequence* component. A templated input has no protein `seq` component to
hang one on (`msa_path` is a SequenceComponent field -- passing it to a CIF `path`
component is a TypeError in atomworks, and hanging it on a DNA component makes rf3
raise "Unsupported chain type for MSAs: polydeoxyribonucleotide"). So a templated
input prepared directly ALWAYS auto-routes to the paid pipeline (~$2.46/fold).
This is a second, nastier instance of pecli issue #174.

pecli's documented prepare -> review/edit -> submit flow gets there anyway:

  1. prepare the UNTEMPLATED carrier-bearing input  -> a BARE run
  2. swap in the TEMPLATED json + drop the template CIF in as an `aux_files`
     companion under templates/  (stage_templated_run.py does this and re-runs
     submit-time validation locally, free)
  3. submit -- which re-checks only the input extension and config.json, and does
     NOT re-run the auto-route decision, so the run stays bare and MSA-free

Templates come from make_predicted_templates.py: the protein chain of that
protein's own prior untemplated prediction, never a crystal chain. See that script
for why.

Usage:
    python submit_templated_folds.py \
        --templated-dir   folds/rf3_templated \
        --untemplated-dir folds/rf3 \
        --template-dir    templates_predicted \
        --unit-cost 0.022 --max-spend 2.50
"""
from __future__ import annotations
import argparse
import json
import os
import re
import subprocess
import sys

RUN_ID_RE = re.compile(r"(?:Prepared|prepared)\s+\S+\s+(?:run|pipeline)\s+([0-9a-f]{6})")
STAGING_RE = re.compile(r"(/\S*?/\.pecli/staging/[^\s/]+)/?")
HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd, timeout=900):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--templated-dir", required=True,
                    help="templated fold inputs + folds_manifest.json")
    ap.add_argument("--untemplated-dir", required=True,
                    help="the carrier-bearing untemplated inputs, used only to get a BARE prepare")
    ap.add_argument("--template-dir", required=True, help="predicted protein-only CIFs")
    ap.add_argument("--unit-cost", type=float, required=True)
    ap.add_argument("--max-spend", type=float, required=True)
    ap.add_argument("--group", default="oracle-controls-templated")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--proteins", default=None,
                    help="comma-separated control labels to restrict to, e.g. "
                         "'LambdaRep,Engrailed'. Keeps each protein's WHOLE 8-target row, "
                         "which is required for ΔminPAE and argmin to be computable -- "
                         "on-target folds alone would only show the absolute-minPAE move and "
                         "could not detect templating trading discrimination for confidence.")
    args = ap.parse_args()

    mpath = os.path.join(os.path.abspath(args.templated_dir), "folds_manifest.json")
    manifest = json.load(open(mpath))
    pending = [m for m in manifest if not m.get("run_id")]
    if args.proteins:
        want = {p.strip() for p in args.proteins.split(",") if p.strip()}
        unknown = want - {m["protein"] for m in manifest}
        if unknown:
            print(f"unknown protein label(s): {sorted(unknown)}")
            return 1
        pending = [m for m in pending if m["protein"] in want]
        print(f"restricted to {sorted(want)}: {len(pending)} folds")
    if args.limit:
        pending = pending[:args.limit]

    projected = len(pending) * args.unit_cost
    print(f"templated rf3  pending={len(pending)}  unit=${args.unit_cost:.3f}  "
          f"projected=${projected:.2f}  cap=${args.max_spend:.2f}")
    if projected > args.max_spend:
        n_ok = int(args.max_spend // args.unit_cost)
        print(f"  projected exceeds cap -- submitting at most {n_ok}")
        pending = pending[:n_ok]

    spent, submitted, failed = 0.0, 0, []
    for i, rec in enumerate(pending, 1):
        if spent + args.unit_cost > args.max_spend:
            print(f"  cap reached at ${spent:.2f}; {len(pending) - i + 1} folds unsubmitted")
            break
        fid = rec["fold_id"]
        untmpl = os.path.join(os.path.abspath(args.untemplated_dir), f"{fid}.json")
        tmpl = os.path.join(os.path.abspath(args.templated_dir), rec["fold_input"])
        cif = os.path.join(os.path.abspath(args.template_dir),
                           f"{rec['protein']}_template.cif")
        for p in (untmpl, tmpl, cif):
            if not os.path.isfile(p):
                failed.append((fid, f"missing {p}"))
                break
        else:
            # 1. prepare the UNTEMPLATED input -> guaranteed bare run
            rc, out = run(["pecli", "prepare", "rf3", "--input", untmpl,
                           "--diffusion-batch-size", "1", "--seed", "42",
                           "--description",
                           f"oracle-controls TEMPLATED {fid} (rf3, MSA-free, protein templated)"])
            m, s = RUN_ID_RE.search(out), STAGING_RE.search(out)
            if rc != 0 or not m or not s:
                failed.append((fid, "prepare: " + out.strip()[-250:]))
                continue
            if "MSA pipeline" in out or "auto-routed" in out:
                failed.append((fid, "prepare auto-routed to a PAID MSA pipeline -- "
                                    "the untemplated input lost its MSA carrier"))
                continue
            run_id, sdir = m.group(1), s.group(1)

            # 2. swap in the templated spec + stage the CIF companion (free, validated)
            rc2, out2 = run([sys.executable,
                             os.path.join(HERE, "stage_templated_run.py"), sdir, tmpl, cif])
            if rc2 != 0 or "PASS:" not in out2:
                failed.append((fid, "stage: " + out2.strip()[-250:]))
                continue

            # 3. submit -- does not re-run the auto-route decision
            rc3, out3 = run(["pecli", "submit", run_id, "-y", "--group", args.group])
            if rc3 != 0:
                failed.append((fid, "submit: " + out3.strip()[-250:]))
                continue
            rec["run_id"] = run_id
            spent += args.unit_cost
            submitted += 1
            print(f"  [{i}/{len(pending)}] {fid}: {run_id}  (spent ~${spent:.2f})")
            with open(mpath, "w") as f:
                json.dump(manifest, f, indent=2)

    with open(mpath, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nsubmitted {submitted}  failed {len(failed)}  spend ~${spent:.2f}")
    for fid, err in failed:
        print(f"  ! {fid}: {err}")
    print(f"manifest: {mpath} "
          f"({sum(1 for m in manifest if not m.get('run_id'))} still unsubmitted)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
