#!/usr/bin/env python3
"""Drive one Phase-2.5 arm from diffusion output to a filter CSV.

Stages, each idempotent so the driver can be re-run after a partial failure:

    fetch   pull the arm's rfd3na designs (free; group index, not `pecli results`
            -- see analysis/oracle_controls/fetch_group_results.py for why)
    relax   OpenMM, DNA positionally restrained, protein free (free, CPU)
    mpnn    LigandMPNN inverse folding, protein chain only  ($0.00 measured)
    refold  rf3, protein + the same duplex, MSA-free, NO template ($0.020/design)
    filter  DNA-aligned Ca-RMSD + protein-only RMSD + H-bonds + minPAE -> CSV

CHAIN LAYOUT, the thing that has bitten twice
---------------------------------------------
rfd3na emits the fixed duplex first, so its output (and the relaxed PDB) is
DNA A + DNA B + protein C. The rf3 refold is the OPPOSITE: protein A + DNA B,C. The
LigandMPNN config must therefore name chain C, and the RMSD step pairs DNA strands by
base sequence rather than by chain id. Do not "tidy" either to match the other.

Usage:
    python run_arm_downstream.py --group p25-fixed_cfgon --work-dir <dir> \
        --dna TGAGGAGAGGAG --max-spend 0.60 [--stages fetch,relax,mpnn]
"""
from __future__ import annotations
import argparse
import glob
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..")
RUN_ID_RE = re.compile(r"(?:Prepared|prepared)\s+\S+\s+(?:run|pipeline)\s+([0-9a-f]{6})")
PECLI_PY = "/Library/Frameworks/Python.framework/Versions/3.12/bin/python3.12"
ALL_STAGES = ("fetch", "relax", "mpnn", "refold", "filter")


def sh(cmd, timeout=3600):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def stage_fetch(a, W):
    out = os.path.join(W, "raw_diffusion")
    rc, o = sh([PECLI_PY if os.path.exists(PECLI_PY) else sys.executable,
                os.path.join(REPO, "analysis", "oracle_controls",
                             "fetch_group_results.py"),
                "--group", a.group, "--out-dir", out])
    print(o.strip()[-1500:])
    pdbs = glob.glob(os.path.join(out, "**", "*.pdb"), recursive=True)
    print(f"fetch: {len(pdbs)} design PDB(s)")
    return len(pdbs) > 0


def stage_relax(a, W):
    src = glob.glob(os.path.join(W, "raw_diffusion", "**", "*.pdb"), recursive=True)
    dest = os.path.join(W, "relaxed")
    os.makedirs(dest, exist_ok=True)
    ok, failed = 0, []
    for p in sorted(src):
        name = os.path.splitext(os.path.basename(p))[0]
        out = os.path.join(dest, f"{name}.pdb")
        if os.path.isfile(out):
            ok += 1
            continue
        rc, o = sh([sys.executable, os.path.join(HERE, "relax_openmm.py"),
                    "--complex", p, "--out", out, "--dna-chains", "A,B"], timeout=1800)
        if rc != 0 or not os.path.isfile(out):
            failed.append((name, o.strip()[-200:]))
        else:
            ok += 1
    print(f"relax: {ok} relaxed, {len(failed)} failed")
    for n, e in failed[:5]:
        print(f"  ! {n}: {e}")
    return ok > 0


def stage_mpnn(a, W):
    pdbs = sorted(glob.glob(os.path.join(W, "relaxed", "*.pdb")))
    submitted = _load(W, "mpnn_runs.json")
    done = {r["backbone"] for r in submitted}
    for p in pdbs:
        name = os.path.splitext(os.path.basename(p))[0]
        if name in done:
            continue
        rc, o = sh(["pecli", "prepare", "ligandmpnn", "--input", p,
                    "--chains-to-design", "C",          # the DESIGNED chain; see docstring
                    "--temperature", "0.1",
                    "--batch-size", str(a.seqs_per_backbone), "--number-of-batches", "1",
                    "--description", f"phase2.5 {a.group} mpnn :: {name}"])
        m = RUN_ID_RE.search(o)
        if rc != 0 or not m:
            print(f"  ! {name}: prepare: {o.strip()[-200:]}")
            continue
        rc2, o2 = sh(["pecli", "submit", m.group(1), "-y", "--group", f"{a.group}-mpnn"])
        if rc2 != 0:
            print(f"  ! {name}: submit: {o2.strip()[-200:]}")
            continue
        submitted.append({"backbone": name, "run_id": m.group(1)})
        # Written after EVERY submission, not once at the end: the record of what has
        # already been paid for must survive an interruption, or a re-run resubmits
        # everything it cannot remember.
        _save(W, "mpnn_runs.json", submitted)
    _save(W, "mpnn_runs.json", submitted)
    print(f"mpnn: {len(submitted)} run(s) submitted (LigandMPNN measured at $0.00)")
    return bool(submitted)


def stage_refold(a, W):
    raw = os.path.join(W, "raw_mpnn")
    rc, o = sh([PECLI_PY if os.path.exists(PECLI_PY) else sys.executable,
                os.path.join(REPO, "analysis", "oracle_controls",
                             "fetch_group_results.py"),
                "--group", f"{a.group}-mpnn", "--out-dir", raw])
    print(o.strip()[-600:])
    folds = os.path.join(W, "folds_refold")
    rc, o = sh([sys.executable, os.path.join(HERE, "build_fold_inputs.py"),
                "--ligandmpnn-dir", raw, "--dna", a.dna, "--out-dir", folds])
    print(o.strip()[-800:])
    if rc != 0:
        return False
    man = json.load(open(os.path.join(folds, "folds_manifest.json")))
    projected = 0.020 * len(man)
    print(f"refold: {len(man)} fold(s), projected ${projected:.2f}, cap ${a.max_spend:.2f}")
    if projected > a.max_spend:
        print("  REFUSING: projected refold spend exceeds --max-spend", file=sys.stderr)
        return False
    spent = 0.0
    for rec in man:
        if rec.get("run_id"):
            continue
        if spent + 0.020 > a.max_spend:
            print(f"  cap reached at ${spent:.2f}")
            break
        rc, o = sh(["pecli", "prepare", "rf3", "--input",
                    os.path.join(folds, rec["fold_input"]),
                    "--diffusion-batch-size", "1", "--seed", "42",
                    "--description", f"phase2.5 {a.group} refold :: {rec['fold_id']}"])
        m = RUN_ID_RE.search(o)
        if rc != 0 or not m:
            print(f"  ! {rec['fold_id']}: prepare: {o.strip()[-200:]}")
            continue
        if "MSA pipeline" in o or "auto-routed" in o:
            print(f"  ! {rec['fold_id']}: auto-routed to a PAID MSA pipeline; skipping")
            continue
        rc2, o2 = sh(["pecli", "submit", m.group(1), "-y",
                      "--group", f"{a.group}-refold"])
        if rc2 != 0:
            print(f"  ! {rec['fold_id']}: submit: {o2.strip()[-200:]}")
            continue
        rec["run_id"] = m.group(1)
        spent += 0.020
        json.dump(man, open(os.path.join(folds, "folds_manifest.json"), "w"), indent=2)
    json.dump(man, open(os.path.join(folds, "folds_manifest.json"), "w"), indent=2)
    print(f"refold: submitted, spend ~${spent:.2f}")
    return True


def stage_filter(a, W):
    raw = os.path.join(W, "raw_refold")
    rc, o = sh([PECLI_PY if os.path.exists(PECLI_PY) else sys.executable,
                os.path.join(REPO, "analysis", "oracle_controls",
                             "fetch_group_results.py"),
                "--group", f"{a.group}-refold", "--out-dir", raw])
    print(o.strip()[-400:])
    fm = os.path.join(W, "filter_manifest.json")
    rc, o = sh([sys.executable, os.path.join(HERE, "build_filter_manifest.py"),
                "--refold-manifest", os.path.join(W, "folds_refold",
                                                  "folds_manifest.json"),
                "--relaxed-dir", os.path.join(W, "relaxed"),
                "--raw-dir", raw, "--out", fm])
    print(o.strip()[-800:])
    if rc != 0:
        return False
    out_csv = os.path.join(W, "all_designs.csv")
    rc, o = sh([sys.executable, os.path.join(HERE, "filter_binder_block.py"),
                "--manifest", fm, "--stage", "pre_resample",
                "--out", os.path.join(W, "passers.csv"),
                "--oracle-comparison", out_csv])
    print(o.strip()[-1000:])
    return rc == 0


def _load(W, name):
    p = os.path.join(W, name)
    return json.load(open(p)) if os.path.isfile(p) else []


def _save(W, name, obj):
    json.dump(obj, open(os.path.join(W, name), "w"), indent=2)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--dna", required=True)
    ap.add_argument("--max-spend", type=float, required=True,
                    help="cap on the REFOLD stage, the only paid downstream step")
    ap.add_argument("--seqs-per-backbone", type=int, default=1,
                    help="1 at this scale: refolding is 89%% of a backbone's cost, and "
                         "the smoke test measured backbone quality dominating sequence "
                         "choice. The paper uses 5.")
    ap.add_argument("--stages", default=",".join(ALL_STAGES))
    args = ap.parse_args()

    W = os.path.abspath(args.work_dir)
    os.makedirs(W, exist_ok=True)
    stages = [s.strip() for s in args.stages.split(",") if s.strip()]
    bad = set(stages) - set(ALL_STAGES)
    if bad:
        ap.error(f"unknown stage(s) {sorted(bad)}; choose from {ALL_STAGES}")

    fns = {"fetch": stage_fetch, "relax": stage_relax, "mpnn": stage_mpnn,
           "refold": stage_refold, "filter": stage_filter}
    for s in stages:
        print(f"\n===== {args.group} :: {s} =====")
        if not fns[s](args, W):
            print(f"stage {s} did not complete; stopping", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
