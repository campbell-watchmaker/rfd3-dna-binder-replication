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
import random
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


def purine_chain(candidates):
    """Which strand is the purine-rich one, counted from the data (not assumed).

    For a poly-purine target like PRNP (TGAGGAGAGGAG) one strand carries nearly all
    the major-groove information: every G contributes N7+O6 and every A contributes
    N7+N6, while the complementary strand offers only C-N4 / T-O4. Conditioning on
    both strands would doubly constrain the same base pairs.
    """
    counts = {}
    for c in candidates:
        counts.setdefault(c["chain"], [0, 0])
        counts[c["chain"]][0 if c["res_name"] in ("DA", "DG") else 1] += 1
    return max(counts, key=lambda ch: counts[ch][0])


def subset_for_ori(candidates, bp_start, bp_end, n_central, strand):
    """Pick the handful of atoms this ori's spec should condition on.

    WHY SUBSET AT ALL. The generator can emit every candidate atom (36 for a 12-bp
    PRNP duplex).

    NOTE ON PROVENANCE, corrected 2026-08-05. This docstring used to assert that "the
    paper conditions on a selected subset, e.g. the N7/O6 of the central G/A run".
    That is NOT in the paper. A full-text read found exactly one Methods sentence on
    the subject -- "Hydrogen bond conditioning was applied during generation on
    candidate major groove donor and acceptor atoms (Fig. S1)" -- with no count, no
    atom names, no strand, and Fig. S1 unreachable. The over-constraining rationale is
    ours, not theirs.

    What the paper DOES say, in Results, is that the constraint set is VARIED:
    "we sample a variety of placements of the protein center of mass ... and a diverse
    set of hydrogen bond (Hbond) condition constraints". So a single fixed subset is
    the wrong model regardless of its size, and the rule below is a stand-in until
    sampling is implemented. Measured: 6 vs 8 atoms was indistinguishable over 100
    refolds (see analysis, commit 76350e1), which is consistent with the count not
    being the thing that matters.

    THE RULE, in three parts:
      1. purine strand only (see purine_chain) -- the information-bearing face;
      2. drop the duplex's terminal base pairs, where the predicted duplex frays and
         where a real binder's contacts are least reliable;
      3. of what remains inside this ori's own bp window, take the CENTRAL n bases,
         so each spec's constraints sit inside the span its ori token describes.

    Both atoms of a base are kept together, never one alone: G N7+O6 and A N7+N6 are
    the bidentate pairs Arg and Asn actually form against a purine, so splitting them
    would specify a weaker and less physical constraint.
    """
    all_ids = [c["res_id"] for c in candidates]
    duplex_lo, duplex_hi = min(all_ids), max(all_ids)
    chains = {purine_chain(candidates)} if strand == "purine" else \
        {c["chain"] for c in candidates}

    usable = sorted({c["res_id"] for c in candidates
                     if c["chain"] in chains
                     and bp_start <= c["res_id"] <= bp_end
                     and duplex_lo < c["res_id"] < duplex_hi})
    if n_central and len(usable) > n_central:
        off = (len(usable) - n_central) // 2
        usable = usable[off:off + n_central]
    keep = set(usable)
    return [c for c in candidates if c["chain"] in chains and c["res_id"] in keep]


def sample_for_ori(candidates, bp_start, bp_end, rng, n_bases_range=(2, 5),
                   strand="either"):
    """Draw ONE design's H-bond constraint set at random from this ori's window.

    WHY SAMPLING RATHER THAN A FIXED SUBSET
    ---------------------------------------
    The paper's own description of the design campaign is that the constraint set is
    varied per design, not fixed:

        "we sample a variety of placements of the protein center of mass relative to
         the DNA target using the RFD3 ori token feature, and a diverse set of
         hydrogen bond (Hbond) condition constraints"

    So `subset_for_ori()` -- one deterministic subset reused for every design off an ori
    -- is the wrong SHAPE regardless of how many atoms it picks, and our measurement is
    consistent with that: 6 vs 8 atoms was indistinguishable over 100 refolds, i.e. the
    count was not the operative variable. It is kept as the `fixed` branch so the
    existing arms stay reproducible.

    THE MECHANISM IS OURS. The paper says the set is diverse; it does not say how it is
    drawn. What is sampled here:

      * how many base positions (uniform over n_bases_range, inclusive) -- so a design
        may be loosely or tightly constrained, which is the diversity that matters;
      * which positions, uniformly without replacement from this ori's window;
      * which strand, per design, when strand="either".

    Both atoms of a base are always kept together -- G N7+O6 and A N7+N6 are the
    bidentate pairs Arg and Asn actually form against a purine, so splitting them
    specifies a weaker and less physical constraint. Terminal base pairs of the duplex
    are excluded, where the predicted duplex frays.

    Two upstream (foundry, not paper) signals bound the count from above: training
    subsamples to roughly a third of candidate atoms, and the shipped 14-bp
    na_binder_design.json example specifies 16 atoms over 4 base positions across BOTH
    strands. The default 2-5 bases spans 4-10 atoms against PRNP's 36 candidates, which
    straddles that third.

    Returns (kept_candidates, draw_record) -- the record is written into the manifest so
    any design can be traced back to the exact constraint set that produced it.
    """
    all_ids = [c["res_id"] for c in candidates]
    duplex_lo, duplex_hi = min(all_ids), max(all_ids)

    if strand == "either":
        chains = sorted({c["chain"] for c in candidates})
        chosen = [rng.choice(chains)]
    elif strand == "purine":
        chosen = [purine_chain(candidates)]
    else:
        chosen = sorted({c["chain"] for c in candidates})
    chosen_set = set(chosen)

    usable = sorted({c["res_id"] for c in candidates
                     if c["chain"] in chosen_set
                     and bp_start <= c["res_id"] <= bp_end
                     and duplex_lo < c["res_id"] < duplex_hi})
    lo, hi = n_bases_range
    n = rng.randint(lo, hi) if usable else 0
    n = min(n, len(usable))
    keep_ids = sorted(rng.sample(usable, n)) if n else []
    kept = [c for c in candidates
            if c["chain"] in chosen_set and c["res_id"] in set(keep_ids)]
    return kept, {
        "mode": "random",
        "strands": chosen,
        "n_bases_drawn": n,
        "n_bases_available": len(usable),
        "base_positions": keep_ids,
        "n_atoms": len(kept),
        "atoms": sorted(f"{c['chain']}{c['res_id']}:{c['atom']}" for c in kept),
    }


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
    ap.add_argument("--hbond-central-bases", type=int, default=3,
                    help="condition on the central N bases of each ori window "
                         "(0 = every candidate, which over-constrains diffusion)")
    ap.add_argument("--hbond-strand", choices=["purine", "both", "either"],
                    default="purine",
                    help="purine (default for --hbond-sampling fixed): the "
                         "information-bearing strand only. 'either' picks one strand per "
                         "design at random and is only meaningful with sampling.")
    ap.add_argument("--hbond-sampling", choices=["fixed", "random"], default="fixed",
                    help="fixed (default): one deterministic subset per ori, the "
                         "historical behaviour, kept so existing arms stay reproducible. "
                         "random: draw an independent constraint set PER DESIGN, which is "
                         "the shape the paper describes ('a diverse set of hydrogen bond "
                         "condition constraints'). The drawing mechanism is ours -- the "
                         "paper states the diversity, not the method.")
    ap.add_argument("--designs-per-ori", type=int, default=1,
                    help="with --hbond-sampling random, how many independently-sampled "
                         "specs to emit per ori (paper scale: ~5100)")
    ap.add_argument("--hbond-bases-min", type=int, default=2)
    ap.add_argument("--hbond-bases-max", type=int, default=5,
                    help="sampled base-position count range, inclusive. 2-5 bases spans "
                         "4-10 atoms against PRNP's 36 candidates, straddling foundry's "
                         "~one-third training subsample rate.")
    ap.add_argument("--seed", type=int, default=42,
                    help="seeds the constraint draw. Recorded per spec, so any design "
                         "traces back to the exact set that produced it.")
    args = ap.parse_args()

    cond = json.load(open(args.conditioning))
    cands = cond["hbond_candidates"]
    dna_ranges = _dna_chain_ranges(cands)
    sampling = args.hbond_sampling == "random"
    if not sampling and args.designs_per_ori != 1:
        ap.error("--designs-per-ori > 1 only makes sense with --hbond-sampling random; "
                 "the fixed subset is deterministic, so every spec would be identical")
    if args.hbond_bases_min > args.hbond_bases_max:
        ap.error("--hbond-bases-min must not exceed --hbond-bases-max")
    rng = random.Random(args.seed)

    os.makedirs(args.out_dir, exist_ok=True)
    manifest = []
    for k, tok in enumerate(cond["ori_tokens"], start=1):
        for d in range(1, args.designs_per_ori + 1):
            name = f"{args.design_name}_ori{k}" + (f"_h{d}" if sampling else "")
            if sampling:
                sub, draw = sample_for_ori(
                    cands, tok["bp_start"], tok["bp_end"], rng,
                    (args.hbond_bases_min, args.hbond_bases_max), args.hbond_strand)
            else:
                # Per-ori subset: each spec conditions only on atoms inside its own span.
                sub = subset_for_ori(cands, tok["bp_start"], tok["bp_end"],
                                     args.hbond_central_bases, args.hbond_strand)
                draw = {"mode": "fixed", "n_atoms": len(sub),
                        "central_bases": args.hbond_central_bases,
                        "strand": args.hbond_strand,
                        "atoms": sorted(f"{c['chain']}{c['res_id']}:{c['atom']}"
                                        for c in sub)}
            donor = _group_hbond(sub, "donor")
            acceptor = _group_hbond(sub, "acceptor")
            print(f"{name}: conditioning on {len(sub)} of {len(cands)} candidate atoms "
                  f"-- {', '.join(sorted({c['chain'] + str(c['res_id']) + ' ' + c['res_name'][-1] for c in sub}))}")
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
                "ori_index": k,
                "ori_bp_range": [tok["bp_start"], tok["bp_end"]],
                "ori_token": [round(float(x), 3) for x in tok["ori_xyz"]],
                "hbond_draw": draw,
            })

    with open(os.path.join(args.out_dir, "manifest.json"), "w") as f:
        json.dump({
            "design_name": args.design_name,
            "protein_len": args.protein_len,
            "duplex_cif": args.duplex_cif,
            "hbond_sampling": args.hbond_sampling,
            "hbond_seed": args.seed if sampling else None,
            "n_specs": len(manifest),
            "specs": manifest,
        }, f, indent=2)
    counts = [m["hbond_draw"]["n_atoms"] for m in manifest]
    print(f"\nwrote manifest with {len(manifest)} specs "
          f"({len(cond['ori_tokens'])} ori x {args.designs_per_ori} draw(s))")
    print(f"  hbond_sampling={args.hbond_sampling}  atoms per spec: "
          f"min {min(counts)}, max {max(counts)}, mean {sum(counts) / len(counts):.1f} "
          f"of {len(cands)} candidates")
    if sampling:
        print(f"  seed {args.seed}; every spec's exact draw is recorded in "
              "manifest.json under hbond_draw")
    else:
        print("  NOTE: a single fixed subset is not the shape the paper describes "
              "('a diverse set of hydrogen bond condition constraints'). "
              "Use --hbond-sampling random.")


if __name__ == "__main__":
    main()
