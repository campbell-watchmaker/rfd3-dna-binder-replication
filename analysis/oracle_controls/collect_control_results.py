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
  openfold3 <name>_seed_<s>_sample_<n>_confidences.json (written whenever the
            `pae_enabled` preset is on, which is the default)
  esmfold2  TWO files, not one: <id>_pae.npy (float16 [N,N] matrix, 0-32 A) plus
            <id>_pae_tokens.json (the per-token chain map). Both only exist if the
            run set output.emit_pae=true (--emit-pae true); esmfold2 emitted
            scalars only before pecli 41d504a. The .npy goes in `pae_path`, the
            sidecar in `pae_tokens_path`.

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
    """Locate the per-token PAE file under a downloaded run directory.

    Returns a path, or None. For esmfold2 the caller must ALSO call
    find_pae_tokens(): its PAE is two files, not one.
    """
    if oracle == "esmfold2":
        # <id>_pae.npy (float16 [N,N], 0-32 A) -- only present when the run was
        # submitted with --emit-pae true (pecli #175); before that patch esmfold2
        # emitted scalars only and minPAE was uncomputable.
        cands = glob.glob(os.path.join(root, "**", "*_pae.npy"), recursive=True)
        if not cands:
            return None
        cands.sort(key=lambda p: (p.count(os.sep), len(p)))
        return cands[0]
    if oracle == "openfold3":
        # <fold>/seed_<s>/<fold>_seed_<s>_sample_<n>_confidences.json holds
        # plddt/pde/pae; the sibling *_confidences_aggregated.json is the summary
        # (scalars only) and must not be picked as the PAE.
        cands = [p for p in glob.glob(os.path.join(root, "**", "*_confidences.json"),
                                      recursive=True)
                 if "aggregated" not in os.path.basename(p)]
        if not cands:
            return None
        cands.sort()  # sample_1 is the first (and, at 1 sample/seed, only) one
        return cands[0]
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


def find_pae_tokens(root):
    """esmfold2 only: <id>_pae_tokens.json, the per-token chain map that goes with
    the .npy matrix (it carries token_chain_ids + a `mapping` provenance string)."""
    cands = glob.glob(os.path.join(root, "**", "*_pae_tokens.json"), recursive=True)
    if not cands:
        return None
    cands.sort(key=lambda p: (p.count(os.sep), len(p)))
    return cands[0]


def find_summary(root, oracle):
    if oracle == "esmfold2":
        # complex_confidence.json carries ptm/iptm/pair_chains_iptm and, with
        # emit_pae, mean_pae (+ pair_chains_mean_pae when the token map is exact).
        pat = "complex_confidence.json"
    elif oracle == "openfold3":
        pat = "*_confidences_aggregated.json"
    else:
        pat = "*_summary_confidences.json" if oracle == "rf3" else "*summary_confidence*.json"
    cands = glob.glob(os.path.join(root, "**", pat), recursive=True)
    if not cands:
        return None
    cands.sort(key=lambda p: (p.count(os.sep), len(p)))
    return cands[0]


def _merge_write(mpath, updated):
    """Write our fields back WITHOUT clobbering concurrent edits.

    Learned the hard way: this collector and submit_control_folds.py both write
    the same manifest. Reading it, spending minutes downloading, then writing the
    whole thing back silently discarded 29 run_ids the submitter had added in the
    meantime -- the jobs were live and billed but no longer tracked. So re-read
    from disk at write time and only overwrite the specific keys this tool owns.
    """
    on_disk = json.load(open(mpath))
    by_id = {r["fold_id"]: r for r in on_disk}
    OWNED = ("pae_path", "pae_tokens_path", "iptm", "ptm", "plddt", "mean_pae")
    for rec in updated:
        tgt = by_id.get(rec["fold_id"])
        if tgt is None:
            continue
        for k in OWNED:
            if rec.get(k) is not None:
                tgt[k] = rec[k]
    with open(mpath, "w") as f:
        json.dump(on_disk, f, indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--download-dir", required=True)
    ap.add_argument("--skip-status", action="store_true",
                    help="don't call `pecli status` per fold; assume anything already "
                         "downloaded is complete. One status call per fold dominates "
                         "runtime on a large batch, and is redundant once the batch is "
                         "known to have succeeded.")
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
        dest = os.path.join(args.download_dir, rec["fold_id"])
        already = os.path.isdir(dest) and os.listdir(dest)
        if args.skip_status:
            if not already:
                counts["NOT_DOWNLOADED"] = counts.get("NOT_DOWNLOADED", 0) + 1
                continue
            counts["ASSUMED_OK"] = counts.get("ASSUMED_OK", 0) + 1
        else:
            st = status_of(rec["run_id"])
            counts[st] = counts.get(st, 0) + 1
            if st != "SUCCEEDED":
                continue
        if not already:
            rc, out = run(["pecli", "results", rec["run_id"], "--out", dest])
            if rc != 0:
                no_pae.append((rec["fold_id"], "download failed: " + out.strip()[-200:]))
                continue
        pae = find_pae(dest, rec["oracle"])
        if not pae:
            hint = ""
            if rec["oracle"] == "protenix":
                hint = " -- protenix run lacked --need-atom-confidence true"
            elif rec["oracle"] == "esmfold2":
                hint = " -- esmfold2 run lacked --emit-pae true"
            no_pae.append((rec["fold_id"], "no per-token PAE in output" + hint))
            continue
        rec["pae_path"] = os.path.abspath(pae)
        if rec["oracle"] == "esmfold2":
            # The matrix alone is unsliceable: the chain of each token lives in the
            # sidecar. Record it, but do NOT make it fatal here -- the metrics step
            # decides whether the map is usable (it can fall back positionally for a
            # monoatomic-ion composition) and is the right place to refuse.
            tok = find_pae_tokens(dest)
            if tok:
                rec["pae_tokens_path"] = os.path.abspath(tok)
            else:
                no_pae.append((rec["fold_id"],
                               "found _pae.npy but no _pae_tokens.json sidecar"))
        s = find_summary(dest, rec["oracle"])
        if s:
            try:
                sd = json.load(open(s))
                # each oracle names the same scalar differently: rf3/protenix
                # write overall_plddt, openfold3 writes avg_plddt
                for k_out, keys in (("iptm", ("iptm",)), ("ptm", ("ptm",)),
                                    ("plddt", ("overall_plddt", "avg_plddt", "plddt")),
                                    ("mean_pae", ("mean_pae",))):
                    for k_in in keys:
                        if k_in in sd:
                            rec[k_out] = round(float(sd[k_in]), 4)
                            break
            except (OSError, ValueError, TypeError):
                pass
        collected += 1

    _merge_write(mpath, manifest)

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
