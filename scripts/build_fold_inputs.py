#!/usr/bin/env python3
"""Build rf3 refold inputs from LigandMPNN output + the target DNA (binder block, Stage 6).

Each designed sequence is folded WITH the target duplex so the oracle predicts the
protein-DNA complex, giving both the refold structure (for the DNA-aligned Ca-RMSD
self-consistency gate) and the interface ipTM.

WHAT CHANGED, 2026-08-06
------------------------
This script used to emit `{"id": ..., "chains": [{"id","type","sequence"}, ...]}`.
**No pecli oracle accepts that shape.** The smoke test worked around it with an ad-hoc
emitter and never committed the fix, so `PIPELINE.md` pointed at a script that could
not have produced a single successful fold. It now emits the rf3 component shape that
the smoke test's 50 refolds actually ran on, and the manifest
`scripts/build_filter_manifest.py` consumes.

Two schema facts, both verified against pecli/foundry source rather than assumed
(the long form is in analysis/oracle_controls/build_control_folds.py):

  * `chain_type` is written EXPLICITLY on every DNA strand. It is optional and
    otherwise inferred from the alphabet, but the all-RNA branch is tested before the
    all-DNA branch, so a T-free strand silently folds as RNA.
  * `_pecli_rf3_msa_a3m` carries a single-sequence a3m. This is load-bearing: omitting
    any MSA makes pecli auto-route the fold to a paid `msa -> rf3` pipeline. It is also
    the scientifically correct setting here -- these are de novo sequences with no
    meaningful alignment.

NO TEMPLATE HERE, DELIBERATELY. The binder block's gate asks whether the designed
sequence INDEPENDENTLY folds back into the backbone it was designed for. Supplying that
backbone as a template hands the fold the answer: RMSD collapses toward zero and
everything passes. Templating belongs only in the specificity block's all-by-all
(see scripts/build_allbyall_inputs.py), which is also where the paper uses it.

THE WT RECORD. LigandMPNN writes the input sequence as record 0 of each .fa, and it is
NOT a design -- folding it wastes money and pollutes the ranking. It is identified
structurally, by the absence of `id=` in its header, rather than by position: the old
`--skip-wt` flag dropped record 0 unconditionally, which silently deletes a real design
from any FASTA that has already been filtered.

Usage:
    # straight from a LigandMPNN result tree (one .fa per backbone)
    python build_fold_inputs.py --ligandmpnn-dir raw/stage5 --dna TGAGGAGAGGAG \
        --out-dir folds/refold

    # or from a single FASTA / an explicit list
    python build_fold_inputs.py --fasta designs.fasta --dna TGAGGAGAGGAG \
        --out-dir folds/refold --backbone ori1_0_3
"""
from __future__ import annotations
import argparse
import glob
import json
import os
import sys

DNA_CHAIN_TYPE = "polydeoxyribonucleotide"
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]


def read_fasta(path):
    """[(header_name, seq)] -- header_name is the first whitespace-delimited token.

    Kept as the module's public FASTA reader because scripts/build_allbyall_inputs.py
    imports it.
    """
    recs = []
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


def parse_ligandmpnn_fasta(path, backbone=None):
    """Parse a LigandMPNN .fa into [{backbone, seq_id, seq, overall_confidence, ...}].

    Header form (verified against real output):
        >ori1_0_3, id=1, T=0.1, seed=1, overall_confidence=0.3918,
         ligand_confidence=0.3433, seq_rec=0.1382
    and for the input sequence, which carries no `id=`:
        >ori1_0_3, T=0.1, seed=1, num_res=123, ...

    Records with no `id=` are the WT input and are dropped (see module docstring).
    """
    out = []
    header, seq = None, []

    def flush():
        if header is None:
            return
        fields = {}
        parts = [p.strip() for p in header.split(",")]
        name = parts[0].split()[0] if parts else ""
        for p in parts[1:]:
            if "=" in p:
                k, _, v = p.partition("=")
                fields[k.strip()] = v.strip()
        if "id" not in fields:
            return                              # the WT input sequence, not a design
        def _f(k):
            try:
                return float(fields[k])
            except (KeyError, ValueError):
                return None
        out.append({
            "backbone": backbone or name,
            "seq_id": int(fields["id"]),
            "seq": "".join(seq),
            "overall_confidence": _f("overall_confidence"),
            "ligand_confidence": _f("ligand_confidence"),
            "seq_rec": _f("seq_rec"),
        })

    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            flush()
            header, seq = line[1:], []
        elif line:
            seq.append(line)
    flush()
    return out


def build_rf3_complex(fold_id, seq, sense, anti):
    """rf3 spec: one entry, protein chain A + DNA chains B/C. Returns (spec, prot, dna)."""
    a3m = f">query\n{seq}\n"                    # single-sequence a3m == MSA-free
    comps = [
        {"seq": seq, "chain_type": "polypeptide(L)", "chain_id": "A",
         "_pecli_rf3_msa_a3m": a3m},
        {"seq": sense, "chain_type": DNA_CHAIN_TYPE, "chain_id": "B"},
        {"seq": anti, "chain_type": DNA_CHAIN_TYPE, "chain_id": "C"},
    ]
    return [{"name": fold_id, "components": comps}], ["A"], ["B", "C"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--ligandmpnn-dir",
                     help="a LigandMPNN result tree; every */seqs/*.fa under it is read, "
                          "one file per designed backbone")
    src.add_argument("--fasta", help="a single FASTA (LigandMPNN-style headers preferred)")
    ap.add_argument("--dna", required=True, help="target DNA sense strand 5'->3'")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--backbone", default=None,
                    help="override the backbone name (--fasta only; defaults to the "
                         "header's first token)")
    ap.add_argument("--top-n-per-backbone", type=int, default=None,
                    help="keep only the N highest-overall_confidence sequences per "
                         "backbone. The triage lever for scale: at 1 seq/backbone a "
                         "refold costs $0.020/design instead of $0.10.")
    ap.add_argument("--skip-wt", action="store_true",
                    help="also drop record 0 of a plain FASTA. Not needed for LigandMPNN "
                         "output, where the WT record is identified by its missing id=.")
    args = ap.parse_args()

    dna_sense = args.dna.upper()
    if set(dna_sense) - set("ACGT"):
        ap.error(f"--dna must be ACGT only, got {sorted(set(dna_sense) - set('ACGT'))}")
    dna_anti = revcomp(dna_sense)

    designs = []
    if args.ligandmpnn_dir:
        fas = sorted(glob.glob(os.path.join(args.ligandmpnn_dir, "**", "seqs", "*.fa"),
                               recursive=True))
        if not fas:
            fas = sorted(glob.glob(os.path.join(args.ligandmpnn_dir, "**", "*.fa"),
                                   recursive=True))
        if not fas:
            print(f"ERROR: no .fa files under {args.ligandmpnn_dir}", file=sys.stderr)
            return 1
        for p in fas:
            designs += parse_ligandmpnn_fasta(p)
        print(f"read {len(designs)} designed sequences from {len(fas)} LigandMPNN file(s)")
    else:
        designs = parse_ligandmpnn_fasta(args.fasta, backbone=args.backbone)
        if not designs:
            # a plain FASTA with no LigandMPNN headers: fall back to positional ids
            recs = read_fasta(args.fasta)
            if args.skip_wt and recs:
                recs = recs[1:]
            designs = [{"backbone": args.backbone or n, "seq_id": i + 1, "seq": s,
                        "overall_confidence": None, "ligand_confidence": None,
                        "seq_rec": None}
                       for i, (n, s) in enumerate(recs)]
            print(f"{args.fasta}: no LigandMPNN headers found; read {len(designs)} "
                  "record(s) positionally")
        elif args.skip_wt and designs:
            designs = designs[1:]

    if args.top_n_per_backbone:
        by_bb = {}
        for d in designs:
            by_bb.setdefault(d["backbone"], []).append(d)
        kept = []
        for bb, ds in by_bb.items():
            # None sorts last, so backbones with no confidence keep their file order
            ds.sort(key=lambda d: (d["overall_confidence"] is None,
                                   -(d["overall_confidence"] or 0.0), d["seq_id"]))
            kept += ds[:args.top_n_per_backbone]
        print(f"triage: kept {len(kept)}/{len(designs)} sequences "
              f"({args.top_n_per_backbone}/backbone over {len(by_bb)} backbones)")
        designs = kept

    os.makedirs(args.out_dir, exist_ok=True)
    manifest = []
    for d in designs:
        fold_id = f"{d['backbone']}_s{d['seq_id']}"
        spec, prot, dna = build_rf3_complex(fold_id, d["seq"], dna_sense, dna_anti)
        with open(os.path.join(args.out_dir, f"{fold_id}.json"), "w") as f:
            json.dump(spec, f, indent=2)
        manifest.append({
            "fold_id": fold_id, "oracle": "rf3",
            "fold_input": f"{fold_id}.json",
            "backbone": d["backbone"], "seq_id": d["seq_id"],
            "protein_seq": d["seq"],
            "overall_confidence": d["overall_confidence"],
            "ligand_confidence": d["ligand_confidence"],
            "protein_chains": prot, "protein_chain": prot[0], "dna_chains": dna,
            "protein_copies": 1, "ligands": [],
            "protein_len": len(d["seq"]), "dna_len": len(dna_sense),
            "pae_path": "FILL_AFTER_FOLD", "run_id": None,
        })

    with open(os.path.join(args.out_dir, "folds_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    n_bb = len({m["backbone"] for m in manifest})
    print(f"wrote {len(manifest)} rf3 fold inputs ({n_bb} backbones) -> {args.out_dir}")
    print(f"  manifest: folds_manifest.json (fill run_id/pae_path after the folds return)")
    print(f"  DNA sense 5'->3': {dna_sense}  |  antisense: {dna_anti}")
    print(f"  est. rf3 cost: ${0.020 * len(manifest):.2f} at the measured $0.020/fold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
