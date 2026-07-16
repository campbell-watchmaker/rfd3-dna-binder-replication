#!/usr/bin/env python3
"""Generate rfd3na design-spec JSON(s) for the binder block from a conditioning bundle.

Turns the output of `compute_conditioning.py` (ori tokens + major-groove H-bond
candidate atoms, computed on the folded target duplex) into rfd3na input specs in
the exact schema the RFdiffusion3-NA checkpoint parses, per the upstream foundry
reference (rosettacommons.github.io/foundry/models/rfd3/input.html and the NA
binder tutorial).

Key schema facts this encodes (verified against the foundry docs, not assumed):
  * `ori_token` is a SINGLE [x,y,z] per spec -- it overrides the COM placement of
    the diffused protein. The paper places one ori per 6-bp stretch and runs
    "~5100 scaffolds per ori", i.e. a SEPARATE diffusion run per ori placement.
    So this generator emits one spec PER ori token; sweep over them at submit time.
  * H-bond conditioning uses two InputSelection dicts, `select_hbond_donor` and
    `select_hbond_acceptor`, keyed by DNA residue id ("A6", "B3-4") with
    comma-joined atom-name strings as values ("N7,O6"). Requires HBPLUS installed
    on the GPU side.
  * The DNA is held fixed via `select_fixed_atoms: {"<dna resid range>": "ALL"}`.
  * `contig` lists the fixed DNA chains plus the designed protein length range
    using the InputSelection mini-language.
  * `is_non_loopy: true` biases toward fewer loops (paper setting).

Sampler knobs (num_timesteps, step_scale/noise, gamma_0, CFG) are NOT written here
-- they are pecli `rfd3na` submit-time config (see sampler_config.json), the same
spec/config split pecli uses for every diffusion tool. This file carries the
biology (the design layout + conditioning); the config carries the sampler.

Usage:
    python make_rfd3na_specs.py \
        --conditioning targets/prnp/conditioning.json \
        --duplex-cif   targets/prnp/prnp_duplex.cif \
        --protein-len  120-150 \
        --design-name  prnp_binder \
        --out-dir      specs/binder_block/rfd3na_specs

Emits one `<design-name>_ori<k>.json` per ori token, plus a `manifest.json`
listing them for the submit driver to sweep.
"""
from __future__ import annotations
import argparse
import json
import os


def _resid_key(chain: str, res_id: int) -> str:
    """rfd3na residue id: chain letter immediately followed by number, e.g. 'A6'."""
    return f"{chain}{res_id}"


def _group_hbond(candidates, role):
    """Collapse [{chain,res_id,atom,role}, ...] into {resid: 'atom,atom'} for one role."""
    out: dict[str, list[str]] = {}
    for c in candidates:
        if c["role"] != role:
            continue
        key = _resid_key(c["chain"], c["res_id"])
        out.setdefault(key, [])
        if c["atom"] not in out[key]:
            out[key].append(c["atom"])
    return {k: ",".join(v) for k, v in out.items()}


def _dna_chain_ranges(candidates):
    """Infer per-chain residue ranges present in the duplex, for select_fixed_atoms / contig."""
    by_chain: dict[str, set[int]] = {}
    for c in candidates:
        by_chain.setdefault(c["chain"], set()).add(c["res_id"])
    ranges = {}
    for ch, ids in by_chain.items():
        ranges[ch] = (min(ids), max(ids))
    return ranges


def build_spec(design_name, duplex_cif, protein_len, ori_xyz, hbond_donor, hbond_acceptor, dna_ranges):
    # Fix all DNA atoms; contig = each DNA chain range, chain break, then designed protein length.
    fixed = {f"{ch}{lo}-{hi}": "ALL" for ch, (lo, hi) in dna_ranges.items()}
    dna_contig = ",/0,".join(f"{ch}{lo}-{hi}" for ch, (lo, hi) in dna_ranges.items())
    contig = f"{dna_contig},/0,{protein_len}"

    spec_body = {
        "input": duplex_cif,
        "contig": contig,
        "length": protein_len,
        "select_fixed_atoms": fixed,
        "ori_token": [round(float(x), 3) for x in ori_xyz],
        "is_non_loopy": True,
    }
    if hbond_acceptor:
        spec_body["select_hbond_acceptor"] = hbond_acceptor
    if hbond_donor:
        spec_body["select_hbond_donor"] = hbond_donor

    return {design_name: spec_body}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conditioning", required=True, help="conditioning bundle JSON from compute_conditioning.py")
    ap.add_argument("--duplex-cif", required=True, help="path (as rfd3na will see it) to the folded target duplex")
    ap.add_argument("--protein-len", default="120-150", help="designed protein length range (paper: 120-150)")
    ap.add_argument("--design-name", default="binder")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cond = json.load(open(args.conditioning))
    cands = cond["hbond_candidates"]
    donor = _group_hbond(cands, "donor")
    acceptor = _group_hbond(cands, "acceptor")
    dna_ranges = _dna_chain_ranges(cands)

    os.makedirs(args.out_dir, exist_ok=True)
    manifest = []
    for k, tok in enumerate(cond["ori_tokens"], start=1):
        name = f"{args.design_name}_ori{k}"
        spec = build_spec(
            name, args.duplex_cif, args.protein_len, tok["ori_xyz"],
            donor, acceptor, dna_ranges,
        )
        path = os.path.join(args.out_dir, f"{name}.json")
        with open(path, "w") as f:
            json.dump(spec, f, indent=2)
        manifest.append({
            "spec": os.path.basename(path),
            "design_name": name,
            "ori_bp_range": [tok["bp_start"], tok["bp_end"]],
            "ori_token": [round(float(x), 3) for x in tok["ori_xyz"]],
        })
        print(f"wrote {path}  (ori bp{tok['bp_start']}-{tok['bp_end']})")

    with open(os.path.join(args.out_dir, "manifest.json"), "w") as f:
        json.dump({
            "design_name": args.design_name,
            "protein_len": args.protein_len,
            "duplex_cif": args.duplex_cif,
            "n_specs": len(manifest),
            "specs": manifest,
        }, f, indent=2)
    print(f"wrote manifest with {len(manifest)} specs (one per ori placement)")


if __name__ == "__main__":
    main()
