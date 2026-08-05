#!/usr/bin/env python3
"""Build templated all-by-all fold inputs for the specificity block (rf3).

For each specificity-resampled design, emit one rf3 fold spec per DNA target in the
off-target set, so folding them all and comparing minPAE reveals which sites the
protein reads confidently. Also emits a folds_manifest.json for
scripts/compute_delta_minpae.py.

TEMPLATED BY DEFAULT
--------------------
The design's protein chain is supplied as a structural template, so every target is
scored on the same pose and only the DNA differs. This is the paper's protocol -- the
one place in the whole pipeline it uses templates:

    "Templates were not used throughout the design campaign with the exception of the
     all-by-all folding step in the specificity block ... the most recent AF3
     prediction before the all-by-all folding was used as the template for the
     protein chain."

and it is corroborated on our own natural-TF control panel
(analysis/oracle_controls/RESULTS.md): templating raised ΔminPAE for both TFs tested
(LambdaRep +1.93 -> +3.43, Engrailed +0.26 -> +0.72) while argmin held 3/3 and every
on-target interface stayed inside the motif window -- i.e. it buys ranking quality
without buying it by making everything look confident.

Only the PROTEIN is templated. The DNA is always free: templating the duplex would
hand the fold the protein-DNA docking geometry it is supposed to predict.

Templates come from scripts/../analysis/oracle_controls/make_predicted_templates.py --
a protein-only CIF extracted from that design's own prior prediction. Never a crystal
structure, and never the complex.

WHY rf3 AND NOT protenix
------------------------
Settled empirically over a 4-oracle control panel (see RESULTS.md): on argmin, rf3 4/5
vs protenix 2/5 vs openfold3 0/5, and protenix's binder/non-binder minPAE ranges
actually OVERLAP. rf3 is also ~2x cheaper than protenix and ~14x cheaper than
openfold3. An earlier version of this script emitted a `{"id", "chains":[...]}` shape
for protenix that no pecli oracle actually accepts.

DO NOT reuse this for the BINDER block's refold. That stage's gate is a DNA-aligned
Ca-RMSD self-consistency check -- "does the designed sequence INDEPENDENTLY fold back
into the intended backbone?" -- and templating it would hand the fold the answer,
collapsing RMSD toward zero and passing everything. The paper does not template there
either.

Usage:
    python build_allbyall_inputs.py \
        --design-fasta survivors.fasta \
        --offtargets   specs/specificity_block/offtargets.json \
        --template-dir templates_predicted \
        --out-dir      specs/specificity_block/allbyall_inputs
"""
from __future__ import annotations
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "analysis", "oracle_controls"))
from build_fold_inputs import read_fasta  # noqa: E402  (reuse the FASTA reader)

DNA_CHAIN_TYPE = "polydeoxyribonucleotide"


def build_templated(fold_id, seq, sense, anti, template_cif):
    """rf3 spec: templated protein (from CIF) + free DNA duplex.

    The protein rides as a `path` component whose chain ids come from the CIF, and
    `template_selection` names those chains so only the protein is templated.
    """
    comps = [{"path": template_cif}]
    comps.append({"seq": sense, "chain_type": DNA_CHAIN_TYPE, "chain_id": "B"})
    comps.append({"seq": anti, "chain_type": DNA_CHAIN_TYPE, "chain_id": "C"})
    spec = [{"name": fold_id, "components": comps, "template_selection": ["A"]}]
    return spec, ["A"], ["B", "C"]


def build_untemplated(fold_id, seq, sense, anti, _template_cif=None):
    """Fallback: sequence-only protein, MSA-free. NOT the paper's protocol here."""
    a3m = f">query\n{seq}\n"
    comps = [{"seq": seq, "chain_type": "polypeptide(L)", "chain_id": "A",
              "_pecli_rf3_msa_a3m": a3m},
             {"seq": sense, "chain_type": DNA_CHAIN_TYPE, "chain_id": "B"},
             {"seq": anti, "chain_type": DNA_CHAIN_TYPE, "chain_id": "C"}]
    return [{"name": fold_id, "components": comps}], ["A"], ["B", "C"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--design-fasta", required=True)
    ap.add_argument("--offtargets", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--template-dir", default=None,
                    help="protein-only CIFs, one per design, named <design_id>_template.cif "
                         "(from analysis/oracle_controls/make_predicted_templates.py). "
                         "REQUIRED unless --no-template.")
    ap.add_argument("--no-template", action="store_true",
                    help="emit sequence-only specs. This is NOT the paper's protocol for "
                         "the all-by-all and measurably lowers ΔminPAE; use only to "
                         "reproduce the untemplated comparison.")
    ap.add_argument("--skip-wt", action="store_true")
    args = ap.parse_args()

    if not args.no_template and not args.template_dir:
        ap.error("--template-dir is required unless --no-template is given "
                 "(the all-by-all is templated by default; see the module docstring)")

    off = json.load(open(args.offtargets))
    targets = off["offtargets"]
    if off.get("panel") == "sbs":
        print("WARNING: this off-target set is the single-base-variant sweep. ΔminPAE is a "
              "MINIMUM over off-targets, so ranking on near-identical variants makes it a "
              "near-worst-case statistic. Use --panel ranking to rank.", file=sys.stderr)
    designs = read_fasta(args.design_fasta)
    if args.skip_wt and designs:
        designs = designs[1:]

    os.makedirs(args.out_dir, exist_ok=True)
    builder = build_untemplated if args.no_template else build_templated
    manifest, n, missing = [], 0, []
    for dname, seq in designs:
        tcif = None
        if not args.no_template:
            tcif = os.path.join(os.path.abspath(args.template_dir), f"{dname}_template.cif")
            if not os.path.isfile(tcif):
                missing.append(dname)
                continue
        for t in targets:
            fold_id = f"{dname}__{t['id']}"
            # rf3 resolves the template path inside the container
            spec, prot, dna = builder(fold_id, seq, t["sense"], t["antisense"],
                                      f"/workspace/templates/{os.path.basename(tcif)}"
                                      if tcif else None)
            with open(os.path.join(args.out_dir, f"{fold_id}.json"), "w") as f:
                json.dump(spec, f, indent=2)
            manifest.append({
                "fold_id": fold_id, "design_id": dname, "dna_id": t["id"],
                "kind": t["kind"], "is_on_target": t["kind"] == "on_target",
                "oracle": "rf3", "templated": not args.no_template,
                "template_cif": tcif,
                "fold_input": f"{fold_id}.json",
                "protein_chain": prot[0], "protein_chains": prot, "dna_chains": dna,
                "protein_copies": 1, "ligands": [],
                "protein_len": len(seq), "dna_len": len(t["sense"]),
                "pae_path": "FILL_AFTER_FOLD", "run_id": None,
            })
            n += 1

    with open(os.path.join(args.out_dir, "folds_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    mode = "sequence-only (NOT the paper's protocol)" if args.no_template \
        else "protein templated, DNA free"
    print(f"{len(designs) - len(missing)} designs x {len(targets)} targets = {n} rf3 fold "
          f"inputs -> {args.out_dir}")
    print(f"  mode: {mode}")
    print(f"  manifest: folds_manifest.json ({len(manifest)} entries; fill pae_path after folds)")
    if missing:
        print(f"  SKIPPED {len(missing)} design(s) with no template CIF: {missing[:5]}")
        print("  -> run analysis/oracle_controls/make_predicted_templates.py first")
        return 1
    if not args.no_template:
        print("  NOTE: a templated rf3 input cannot be `pecli prepare`d directly -- see")
        print("        analysis/oracle_controls/submit_templated_folds.py for the")
        print("        prepare-untemplated-then-swap flow it requires.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
