# Specificity block — pipeline runbook (PRNP-site)

The specificity block does explicit **negative design**: it takes binder-block
passers and re-optimizes them to bind the on-target site while *rejecting*
off-target sites, ranking by **ΔminPAE**. This is the step that took the paper
from ~0.5% (binder block alone) to ~3% specific designs (~6× improvement).

Prerequisite: a completed binder block (`specs/binder_block/`) with passing
designs and their on-target minPAE recorded.

## Stage 0 — select entrants (CPU, here)

From the binder-block passers, keep those with **on-target minPAE < 1.25**
(paper gate). These backbones enter the specificity block.

## Stage 1 — build the off-target set (CPU, here)

```bash
python scripts/make_offtarget_set.py \
    --on-target TGAGGAGAGGAG \
    --out specs/specificity_block/offtargets.json
```

Produces 46 off-targets for PRNP: the on-target itself (ΔminPAE reference) +
36 single-base-substitution variants (the paper's fine-grained specificity test,
DBS5 was "specific over 35/40 single-base variants") + 9 unrelated Table 1
decoys.

## Stage 2 — specificity resample (GPU, pecli, per backbone)

Deeper LigandMPNN sampling (100 seq/backbone vs 5 in the binder block) to find
the specificity-optimal sequence:

```bash
for pdb in <entrant backbones>/*_relaxed.pdb; do
    pecli prepare ligandmpnn --input "$pdb" \
        --config specs/specificity_block/ligandmpnn_resample_config.json
    pecli submit <run>
done
```

## Stage 3 — pre-filter fold (GPU, pecli)

Fold resampled sequences against the **on-target** and keep the good ones
(paper: DNA-aligned RMSD < 1.5 Å, ipTM > 0.9) before the expensive all-by-all:

```bash
python scripts/build_fold_inputs.py --fasta <resample>.fasta \
    --dna TGAGGAGAGGAG --out-dir specs/specificity_block/on_fold_inputs
# submit protenix folds; filter with scripts/filter_binder_block.py at the
# tighter specificity gates (--rmsd-gate 1.5 --iptm-gate 0.9)
```

## Stage 4 — templated all-by-all fold (GPU, pecli)

For each surviving design, fold against the on-target **and every off-target**,
using the on-target complex as a structural template so off-targets are scored on
the same pose:

```bash
python scripts/build_allbyall_inputs.py \
    --design-fasta <survivor>.fasta \
    --offtargets specs/specificity_block/offtargets.json \
    --out-dir specs/specificity_block/allbyall_inputs
# submit protenix (or openfold3) folds -- PAE-emitting oracle REQUIRED
# (esmfold2 has no PAE, cannot be used here)
```

## Stage 5 — ΔminPAE ranking (CPU, here)

```bash
python scripts/compute_delta_minpae.py \
    --manifest specs/specificity_block/folds_manifest.json \
    --out results/specificity_block/delta_minpae.csv \
    --per-complex-out results/specificity_block/min_pae_all.csv
```

`folds_manifest.json` lists, per (design, dna_target) fold: the PAE output path,
the target `kind` (on_target / sbs / decoy), and the oracle. The script computes
minPAE per complex and ΔminPAE per design, ranked descending. **Take the top 96
per target** (paper).

## Oracle note

The binder block compares three folders (protenix / openfold3 / esmfold2). The
specificity block **cannot use esmfold2** — ΔminPAE needs a PAE matrix, which
only the AF3-class folders (protenix / openfold3) emit. Use protenix as primary,
openfold3 as cross-check.

## Reference for PRNP

Paper reported the **highest specificity-block hit rate for PRNP: 13/96 specific
designs**. DBS5 (specificity-block design for this site) is specific over 35/40
single-base variants — our SBS panel is built to reproduce exactly that test.
