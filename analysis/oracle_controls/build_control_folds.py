#!/usr/bin/env python3
"""Emit all-by-all fold inputs for the oracle control panel (rf3 and protenix).

Builds one fold input per (control protein x DNA target) and a manifest that
compute_control_metrics.py consumes once the folds return.

SCHEMA NOTES -- both verified against pecli/foundry/protenix source, not assumed
-------------------------------------------------------------------------------
rf3 (RosettaFold3, via `rf3 fold inputs=<spec>`):
    Top level is a LIST of entries; each entry is {name, components}. A DNA
    strand is an ordinary sequence component distinguished by `chain_type`:
        {"seq": "...", "chain_type": "polydeoxyribonucleotide", "chain_id": "B"}
    A metal is a component with no `seq`: {"ccd_code": "ZN", "chain_id": "D"}.
    Copies of one chain are REPEATED components with distinct chain_id (there is
    no count field).

    TRAP: `chain_type` is optional and otherwise inferred from the alphabet, but
    a DNA strand containing no T is inferred as RNA (the all-RNA branch is
    tested before the all-DNA branch) and would silently fold as RNA. We always
    write chain_type explicitly.

protenix (AlphaFold3-class, via `protenix pred`):
    Top level is a LIST. Entities are {"proteinChain"|"dnaSequence"|"ion": {...}}
    and `count` is required on every entity. Copies are a count, not a list of
    ids. A duplex is TWO dnaSequence entities of count 1 each -- count 2 would
    duplicate one strand rather than add its complement.

    Chain letters are assigned in `sequences` order when "id" is omitted (we omit
    it, since "id" is a newer field and the container pins no protenix version).
    The manifest therefore records the expected chain order so the PAE can be
    sliced without guessing.

MSA
---
Both oracles are run MSA-FREE, and this is a deliberate scientific choice, not a
cost shortcut (though it is also much cheaper). The designs this metric will be
applied to are de novo sequences with no meaningful alignment. Giving natural
TFs deep MSAs while designs get none would make the controls strictly easier
than the real task and inflate the metric's apparent discriminative power.

  rf3       carries an inline single-sequence a3m via pecli's private
            `_pecli_rf3_msa_a3m` carrier. Supplying any MSA keeps the run a
            single bare fold; omitting it would auto-route to a paid `msa -> rf3`
            pipeline.
  protenix  carries `unpairedMsa` (see build_protenix; `--use-msa false` alone
            does NOT suppress the paid MSA routing).
  openfold3 carries `_pecli_main_msa_a3m` per protein chain (see build_openfold3).

Usage:
    python build_control_folds.py --oracle rf3       --out-dir folds/rf3
    python build_control_folds.py --oracle protenix  --out-dir folds/protenix
    python build_control_folds.py --oracle openfold3 --out-dir folds/openfold3
    python build_control_folds.py --oracle rf3 --only Engrailed:hd_taatta   # one probe
"""
from __future__ import annotations
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from control_panel import build_panel, verify_panel, load_controls  # noqa: E402

DNA_CHAIN_TYPE = "polydeoxyribonucleotide"


def _chain_letters(n: int):
    """A, B, ... Z, AA, AB, ... in the same order protenix letters them."""
    out = []
    for i in range(n):
        s, k = "", i
        while True:
            s = chr(ord("A") + k % 26) + s
            k = k // 26 - 1
            if k < 0:
                break
        out.append(s)
    return out


def build_rf3(fold_id, protein_seq, copies, ligands, sense, anti):
    """rf3 spec: list of one entry, components in a fixed order."""
    comps, chain_ids = [], []
    letters = iter(_chain_letters(64))
    a3m = f">query\n{protein_seq}\n"  # single-sequence a3m == MSA-free
    prot_chains = []
    for _ in range(copies):
        cid = next(letters)
        comps.append({"seq": protein_seq, "chain_type": "polypeptide(L)",
                      "chain_id": cid, "_pecli_rf3_msa_a3m": a3m})
        prot_chains.append(cid)
        chain_ids.append(cid)
    dna_chains = []
    for strand in (sense, anti):
        cid = next(letters)
        comps.append({"seq": strand, "chain_type": DNA_CHAIN_TYPE, "chain_id": cid})
        dna_chains.append(cid)
        chain_ids.append(cid)
    for ccd in ligands:
        cid = next(letters)
        comps.append({"ccd_code": ccd, "chain_id": cid})
        chain_ids.append(cid)
    spec = [{"name": fold_id, "components": comps}]
    return spec, prot_chains, dna_chains


def build_protenix(fold_id, protein_seq, copies, ligands, sense, anti):
    """protenix spec: list of one entry; chains lettered in sequences order.

    `unpairedMsa` carries a single-sequence a3m. This is load-bearing for cost as
    well as for the MSA-free design: pecli's _should_route_to_inhouse_msa()
    auto-routes any protenix fold WITHOUT a precomputed MSA to a paid
    `msa -> protenix` pipeline, and `--use-msa false` does NOT suppress that
    (verified: it still prepared an in-house MSA pipeline). Supplying unpairedMsa
    makes input_has_precomputed_msa() true, so the fold runs bare.
    """
    a3m = f">query\n{protein_seq}\n"
    seqs = [{"proteinChain": {"sequence": protein_seq, "count": copies,
                              "unpairedMsa": a3m}}]
    seqs.append({"dnaSequence": {"sequence": sense, "count": 1}})
    seqs.append({"dnaSequence": {"sequence": anti, "count": 1}})
    if ligands:
        # one entity per distinct ion so counts stay explicit
        seqs.append({"ion": {"ion": ligands[0], "count": len(ligands)}})
    # chain letters follow entity order, expanded by count
    letters = iter(_chain_letters(64))
    prot_chains = [next(letters) for _ in range(copies)]
    dna_chains = [next(letters), next(letters)]
    spec = [{"name": fold_id, "sequences": seqs}]
    return spec, prot_chains, dna_chains


def build_openfold3(fold_id, protein_seq, copies, ligands, sense, anti):
    """openfold3 query-set spec: {"seeds": [...], "queries": {name: {"chains": [...]}}}.

    Unlike rf3/protenix this top level is an OBJECT, not a list, and a chain is
    typed by `molecule_type` in {PROTEIN, DNA, LIGAND} with an explicit
    `chain_ids` LIST (a chain entry may cover several copies, but we write one
    entry per chain id so each protein copy carries its own MSA and the
    manifest's chain order is unambiguous). Ligands take `ccd_codes`, not
    `ccd_code`.

    `_pecli_main_msa_a3m` is the openfold3 MSA carrier (the analogue of rf3's
    `_pecli_rf3_msa_a3m` / protenix's `unpairedMsa`): without it the fold has no
    precomputed MSA and pecli auto-routes it to a paid `msa -> openfold3`
    pipeline via the ColabFold server (msa.use_msa_server defaults true). A
    single-sequence a3m keeps the run bare and MSA-free, matching the other two
    oracles. The seed lives in the input (`seeds`) as well as in the CLI flags.
    """
    a3m = f">query\n{protein_seq}\n"
    letters = iter(_chain_letters(64))
    chains, prot_chains, dna_chains = [], [], []
    for _ in range(copies):
        cid = next(letters)
        chains.append({"molecule_type": "PROTEIN", "chain_ids": [cid],
                       "sequence": protein_seq, "_pecli_main_msa_a3m": a3m})
        prot_chains.append(cid)
    for strand in (sense, anti):
        cid = next(letters)
        chains.append({"molecule_type": "DNA", "chain_ids": [cid], "sequence": strand})
        dna_chains.append(cid)
    for ccd in ligands:
        cid = next(letters)
        chains.append({"molecule_type": "LIGAND", "chain_ids": [cid], "ccd_codes": [ccd]})
    spec = {"seeds": [42], "queries": {fold_id: {"chains": chains}}}
    return spec, prot_chains, dna_chains


BUILDERS = {"rf3": build_rf3, "protenix": build_protenix, "openfold3": build_openfold3}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--oracle", required=True, choices=sorted(BUILDERS))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--controls", default=None,
                    help="curated_controls.json (default: alongside this script)")
    ap.add_argument("--only", default=None,
                    help="restrict to Protein:dna_id (for a single cheap probe fold)")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    controls = load_controls(args.controls or os.path.join(here, "curated_controls.json"))
    panel = build_panel()

    problems = verify_panel(panel)
    if problems:
        print("PANEL VERIFICATION FAILED -- refusing to emit folds:")
        for p in problems:
            print("  ! " + p)
        return 1

    os.makedirs(args.out_dir, exist_ok=True)
    manifest, n = [], 0
    for c in controls:
        for dna_id, d in panel.items():
            fold_id = f"{c['label']}__{dna_id}"
            if args.only and args.only != f"{c['label']}:{dna_id}":
                continue
            builder = BUILDERS[args.oracle]
            spec, prot_chains, dna_chains = builder(
                fold_id, c["protein_sequence"], c["copies"], c["ligands"],
                d["sense"], d["antisense"])
            path = os.path.join(args.out_dir, f"{fold_id}.json")
            with open(path, "w") as f:
                json.dump(spec, f, indent=2)
            manifest.append({
                "fold_id": fold_id,
                "protein": c["label"],
                "klass": c["klass"],
                "dna_id": dna_id,
                "is_on_target": (c.get("on_target") == dna_id),
                "on_target_for_protein": c.get("on_target"),
                "oracle": args.oracle,
                "fold_input": os.path.basename(path),
                "protein_chains": prot_chains,
                "dna_chains": dna_chains,
                "protein_copies": c["copies"],
                "ligands": c["ligands"],
                "protein_len": c["protein_length"],
                "dna_len": len(d["sense"]),
                "expected_tokens": c["protein_length"] * c["copies"] + 2 * len(d["sense"]),
                "pae_path": "FILL_AFTER_FOLD",
                "run_id": None,
            })
            n += 1

    mpath = os.path.join(args.out_dir, "folds_manifest.json")
    with open(mpath, "w") as f:
        json.dump(manifest, f, indent=2)

    tok = [m["expected_tokens"] for m in manifest]
    print(f"{args.oracle}: {n} fold inputs -> {args.out_dir}")
    print(f"  tokens per complex: min {min(tok)}, max {max(tok)} "
          f"(A10G is fine below ~700; above that pass gpu=large)")
    print(f"  manifest: {mpath}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
