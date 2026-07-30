#!/usr/bin/env python3
"""Download finished control folds and fill `pae_path` (+ confidence scalars) in the manifest.

Polls each submitted run, downloads the ones that have SUCCEEDED, locates the
per-token PAE file, and records it. Safe to re-run: already-collected folds are
skipped, so this can be called repeatedly while a batch drains.

PAE file, per oracle:
  rf3       <name>_confidences.json            (always written)
  protenix  <name>_full_data_sample_<rank>.json (ONLY if the run set
            output.need_atom_confidence=true; otherwise no per-token PAE exists
            and the fold must be resubmitted with the flag)

Usage:
    python collect_control_results.py --manifest folds/rf3/folds_manifest.json \
        --download-dir results_raw/rf3
"""
from __future__ import annotations
import argparse
import glob
import json
import os
import subprocess
import sys


def run(cmd, timeout=900):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def status_of(run_id):
    rc, out = run(["pecli", "status", run_id], timeout=180)
    if rc != 0:
        return "UNKNOWN"
    first = out.strip().splitlines()[0] if out.strip() else ""
    parts = first.split()
    return parts[1] if len(parts) > 1 else "UNKNOWN"


def find_pae(root, oracle):
    """Locate the per-token PAE file under a downloaded run directory."""
    if oracle == "rf3":
        # prefer the top-level (best-ranked) confidences over per-sample copies
        cands = [p for p in glob.glob(os.path.join(root, "**", "*_confidences.json"),
                                      recursive=True)
                 if "summary" not in os.path.basename(p)]
        if not cands:
            return None
        cands.sort(key=lambda p: (p.count(os.sep), len(p)))
        return cands[0]
    cands = glob.glob(os.path.join(root, "**", "*full_data_sample*.json"), recursive=True)
    if not cands:
        return None
    cands.sort()  # sample_0 is the best-ranked
    return cands[0]


def find_summary(root, oracle):
    pat = "*_summary_confidences.json" if oracle == "rf3" else "*summary_confidence*.json"
    cands = glob.glob(os.path.join(root, "**", pat), recursive=True)
    if not cands:
        return None
    cands.sort(key=lambda p: (p.count(os.sep), len(p)))
    return cands[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--download-dir", required=True)
    args = ap.parse_args()

    mpath = os.path.abspath(args.manifest)
    manifest = json.load(open(mpath))
    os.makedirs(args.download_dir, exist_ok=True)

    todo = [m for m in manifest if m.get("run_id")
            and m.get("pae_path") in (None, "", "FILL_AFTE" "R_FOLD")]
    print(f"{len(todo)} submitted folds awaiting collection "
          f"({sum(1 for m in manifest if m.get('run_id')) } submitted of {len(manifest)})")

    counts = {}
    collected, no_pae = 0, []
    for rec in todo:
        st = status_of(rec["run_id"])
        counts[st] = counts.get(st, 0) + 1
        if st != "SUCCEEDED":
            continue
        dest = os.path.join(args.download_dir, rec["fold_id"])
        if not os.path.isdir(dest) or not os.listdir(dest):
            rc, out = run(["pecli", "results", rec["run_id"], "--out", dest])
            if rc != 0:
                no_pae.append((rec["fold_id"], "download failed: " + out.strip()[-200:]))
                continue
        pae = find_pae(dest, rec["oracle"])
        if not pae:
            no_pae.append((rec["fold_id"],
                           "no per-token PAE in output"
                           + (" -- protenix run lacked --need-atom-confidence true"
                              if rec["oracle"] == "protenix" else "")))
            continue
        rec["pae_path"] = os.path.abspath(pae)
        s = find_summary(dest, rec["oracle"])
        if s:
            try:
                sd = json.load(open(s))
                for k_out, k_in in (("iptm", "iptm"), ("ptm", "ptm"),
                                    ("plddt", "overall_plddt")):
                    if k_in in sd:
                        rec[k_out] = round(float(sd[k_in]), 4)
                    elif k_out in sd:
                        rec[k_out] = round(float(sd[k_out]), 4)
            except (OSError, ValueError, TypeError):
                pass
        collected += 1

    with open(mpath, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"statuses: {counts}")
    print(f"collected {collected} PAE file(s)")
    for fid, err in no_pae:
        print(f"  ! {fid}: {err}")
    ready = sum(1 for m in manifest
                if m.get("pae_path") not in (None, "", "FILL_AFTE" "R_FOLD"))
    print(f"manifest now has {ready}/{len(manifest)} folds with a PAE path")
    return 0


if __name__ == "__main__":
    sys.exit(main())
