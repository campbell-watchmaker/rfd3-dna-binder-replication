#!/usr/bin/env python3
"""Prepare + stage + submit one rfd3na arm of the Phase-2.5 A/B comparison, spend-capped.

An rfd3na spec references its target duplex by CONTAINER path (`/workspace/…cif`), and
`pecli prepare` stages only config.json and the spec itself -- so a run submitted without
`stage_rfd3na_target.py` starts and immediately fails on a CIF that was never uploaded.
This drives the whole prepare -> stage -> submit sequence for every spec in an arm, so
the staging step cannot be forgotten (it has been, twice).

Every submission is counted against `--max-spend` BEFORE it is sent, and the cap is
checked per run rather than only up front, so a mid-batch failure cannot overshoot.

Usage:
    python submit_arm_diffusion.py --spec-dir <dir of *.json specs> \
        --config <sampler config.json> --duplex-cif targets/prnp/prnp_duplex.cif \
        --group p25-fixed-cfgon --designs-per-run 10 \
        --unit-cost 0.012 --max-spend 0.30 [--dry-run]
"""
from __future__ import annotations
import argparse
import glob
import json
import os
import re
import subprocess
import sys

RUN_ID_RE = re.compile(r"(?:Prepared|prepared)\s+\S+\s+(?:run|pipeline)\s+([0-9a-f]{6})")
STAGING_RE = re.compile(r"(/\S*?/\.pecli/staging/[^\s/]+)/?")
HERE = os.path.dirname(os.path.abspath(__file__))

# rfd3na's sampler knobs, as pecli actually exposes them (pecli/tools/rfd3na.py).
# `pecli prepare` takes ONE flag per config field -- there is no --config and no
# --design-inputs, though specs/binder_block/PIPELINE.md claimed both until 2026-08-06.
# Config keys are the leaf names; pecli's parser maps hyphens back to underscores.
KNOWN_FIELDS = {
    "diffusion_batch_size", "n_batches", "kind", "num_timesteps", "step_scale",
    "noise_scale", "gamma_0", "use_classifier_free_guidance", "cfg_scale",
    "ckpt_path", "read_sequence_from_sequence_head", "prevalidate_inputs",
    "dump_trajectories", "low_memory_mode", "gpu",
}
# Keys our sampler_config.json carries that rfd3na does NOT expose. `cfg_features` is
# the notable one: it is written into the config with a long rationale, but pecli has no
# such field, so it has never reached the sampler on any run. Skipped explicitly and
# loudly rather than passed and rejected.
UNSUPPORTED = {"cfg_features"}


def config_to_flags(cfg: dict):
    """Turn a sampler-config dict into pecli prepare flags. Returns (flags, skipped)."""
    flags, skipped = [], []
    for k, v in cfg.items():
        if k.startswith("_"):
            continue
        if k in UNSUPPORTED or k not in KNOWN_FIELDS:
            skipped.append(k)
            continue
        flag = "--" + k.replace("_", "-")
        # bools are `--flag value`; a bare flag means True, so pass the value
        # explicitly to keep `false` expressible
        flags += [flag, ("true" if v else "false") if isinstance(v, bool) else str(v)]
    return flags, skipped


def run(cmd, timeout=1800):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spec-dir", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--duplex-cif", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--designs-per-run", type=int, required=True,
                    help="diffusion_batch_size * n_batches, for the cost projection")
    ap.add_argument("--unit-cost", type=float, default=0.012, help="$ per DESIGN")
    ap.add_argument("--max-spend", type=float, required=True)
    ap.add_argument("--description", default=None)
    ap.add_argument("--dry-run", action="store_true",
                    help="prepare and stage only -- free, and validates every input shape")
    ap.add_argument("--out", default=None, help="JSON manifest of submitted runs")
    args = ap.parse_args()

    specs = sorted(p for p in glob.glob(os.path.join(args.spec_dir, "*.json"))
                   if os.path.basename(p) != "manifest.json")
    if not specs:
        print(f"no specs in {args.spec_dir}", file=sys.stderr)
        return 1
    cif = os.path.abspath(args.duplex_cif)
    if not os.path.isfile(cif):
        print(f"no duplex CIF at {cif}", file=sys.stderr)
        return 1

    flags, skipped = config_to_flags(json.load(open(args.config)))
    print(f"sampler flags: {' '.join(flags)}")
    if skipped:
        print(f"  SKIPPED (rfd3na exposes no such field): {sorted(skipped)}")

    per_run = args.designs_per_run * args.unit_cost
    projected = len(specs) * per_run
    print(f"{args.group}: {len(specs)} spec(s) x {args.designs_per_run} designs "
          f"= {len(specs) * args.designs_per_run} designs")
    print(f"  ${per_run:.3f}/run  projected ${projected:.2f}  cap ${args.max_spend:.2f}")
    if projected > args.max_spend and not args.dry_run:
        print(f"  REFUSING: projected ${projected:.2f} exceeds cap ${args.max_spend:.2f}",
              file=sys.stderr)
        return 1

    desc = args.description or f"phase2.5 arm {args.group}"
    spent, submitted, failed = 0.0, [], []
    for i, spec in enumerate(specs, 1):
        name = os.path.splitext(os.path.basename(spec))[0]
        if not args.dry_run and spent + per_run > args.max_spend:
            print(f"  cap reached at ${spent:.2f}; {len(specs) - i + 1} spec(s) unsubmitted")
            break

        rc, out = run(["pecli", "prepare", "rfd3na", "--input", spec, *flags,
                       "--description", f"{desc} :: {name}"])
        m, s = RUN_ID_RE.search(out), STAGING_RE.search(out)
        if rc != 0 or not m or not s:
            failed.append((name, "prepare: " + out.strip()[-300:]))
            continue
        run_id, sdir = m.group(1), s.group(1)

        # the step whose omission silently produces a run that dies on a missing CIF
        rc2, out2 = run([sys.executable, os.path.join(HERE, "stage_rfd3na_target.py"),
                         sdir, cif])
        if rc2 != 0:
            failed.append((name, "stage: " + out2.strip()[-300:]))
            continue

        if args.dry_run:
            print(f"  [{i}/{len(specs)}] {name}: prepared {run_id}, staged OK (NOT submitted)")
            submitted.append({"spec": name, "run_id": run_id, "submitted": False})
            continue

        rc3, out3 = run(["pecli", "submit", run_id, "-y", "--group", args.group])
        if rc3 != 0:
            failed.append((name, "submit: " + out3.strip()[-300:]))
            continue
        spent += per_run
        submitted.append({"spec": name, "run_id": run_id, "submitted": True,
                          "group": args.group, "designs": args.designs_per_run})
        print(f"  [{i}/{len(specs)}] {name}: {run_id}  (spent ~${spent:.2f})")

    if args.out:
        json.dump(submitted, open(args.out, "w"), indent=2)
    verb = "prepared+staged" if args.dry_run else "submitted"
    print(f"\n{verb} {len(submitted)}  failed {len(failed)}  spend ~${spent:.2f}")
    for name, err in failed:
        print(f"  ! {name}: {err}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
