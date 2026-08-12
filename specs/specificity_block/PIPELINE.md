# Specificity block — pipeline runbook (PRNP-site)

The specificity block does explicit **negative design**: it takes binder-block
passers and re-optimizes them to bind the on-target site while *rejecting*
off-target sites, ranking by **ΔminPAE**. This is the step that took the paper
from ~0.5% (binder block alone) to ~3% specific designs (~6× improvement).

Prerequisite: a completed binder block (`specs/binder_block/`) with passing
designs and their on-target minPAE recorded.

## Stage 0 — select entrants (CPU, here)

From the binder-block passers, keep those with on-target **minPAE < 6.6 Å**.

> **This is not the paper's number, and the substitution is deliberate.** The paper
> gates at `minPAE < 1.25`, measured on AF3. Applied literally to rf3 output it admits
> almost nobody: of five real, crystallographically-characterised TFs folded against
> their own cognate sites, **1/5** clears 1.25. The offset is rf3 calibration rather
> than a modelling error — it survived adding MSAs (−0.11 Å), 5× sampling (≈0) and
> protein templating (−0.93 to +0.11 Å). Measured rf3 on-target minPAE: specific TFs
> 1.07–3.88 Å, non-specific binder (Sac7d) 6.60 Å, non-binders (Ubiquitin/GFP)
> 15.7–16.2 Å, so 6.6 Å is where binders separate from non-binders on measured data.
> Full derivation in `docs/replication_log.md`. **Provisional — it rests on 5
> proteins.** The *success* criterion downstream is the paper's own calibration-free
> `ΔminPAE > 0`, which does transfer (4/5 TFs).

minPAE for a binder-block refold comes from `scripts/filter_binder_block.py`
(`--min-pae-out`), which reads the PAE already on disk in each rf3 refold's
`*_confidences.json`.

## Stage 1 — build the off-target set (CPU, here)

```bash
python scripts/make_offtarget_set.py \
    --on-target TGAGGAGAGGAG \
    --panel ranking \
    --out specs/specificity_block/offtargets.json
```

Produces **10** folds per design for PRNP: the on-target (the ΔminPAE reference)
plus the 9 other Table 1 targets, **every one padded to a common 12 bp**.

> **Padding, added 2026-08-06 — a deviation from the paper.** 4 of the 9 Table 1
> decoys are 10 bp against a 12-bp on-target. minPAE is a *minimum* over protein×DNA
> token pairs, so a shorter duplex simply offers fewer pairs to minimise over and is
> systematically disadvantaged as an off-target — biasing ΔminPAE upward for reasons
> that have nothing to do with specificity. Every target is therefore centred in a
> verified-neutral flank, reusing `build_duplex()`/`verify_panel()` from
> `analysis/oracle_controls/control_panel.py`; `verify_panel()` also refuses the panel
> if a padded decoy picks up another target's motif from its flank. This is the same
> confound the oracle control panel was padded to 24 bp to remove, and padding there
> did not degrade discrimination (rf3 argmin 4/5 on padded duplexes). The paper appears
> to fold Table 1 sites at native length. `--no-pad` reproduces that.

> **Corrected 2026-07-31.** This stage used to emit 46 targets, folding the 36
> single-base-substitution variants alongside the decoys. That was wrong on two
> counts. It inflated the all-by-all ~4.6×, and — worse — it changed what the
> metric means: ΔminPAE is a *minimum* over off-targets, so including sequences one
> base from the on-target turns it into a near-worst-case statistic instead of the
> discrimination-against-unrelated-sites statistic the paper reports. The paper's
> all-by-all folds against the on-target plus the other Table 1 targets; its
> "specific over 35/40 single-base variants" claim is wet-lab characterisation of
> one already-selected binder (DBS5), not the ranking panel.
>
> The SBS sweep is still available as `--panel sbs`, for characterising a design
> *after* it has been ranked. Do not rank on it.

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
python scripts/build_fold_inputs.py \
    --ligandmpnn-dir <resample raw dir> \
    --dna TGAGGAGAGGAG --out-dir specs/specificity_block/on_fold_inputs
# submit rf3 folds, then:
python scripts/build_filter_manifest.py \
    --refold-manifest specs/specificity_block/on_fold_inputs/folds_manifest.json \
    --relaxed-dir <entrant backbones> --raw-dir <downloaded results> \
    --out filter_manifest.json
python scripts/filter_binder_block.py --manifest filter_manifest.json \
    --target-dna TGAGGAGAGGAG --stage post_resample \
    --rmsd-gate 1.5 --iptm-gate 0.9 --out results/specificity_block/prefilter.csv
```

`--stage post_resample` is **required** for the ipTM gate to apply at all
(`filter_binder_block.py`), and `--ligandmpnn-dir` reads the resample's `.fa` files
directly, dropping LigandMPNN's WT input record by its missing `id=` rather than by
position.

## Stage 4 — templated all-by-all fold (GPU, pecli)

For each surviving design, fold against the on-target **and every off-target**,
with the design's **protein chain templated** so every target is scored on the same
pose and only the DNA differs.

```bash
python scripts/build_allbyall_inputs.py \
    --design-fasta <survivor>.fasta \
    --offtargets specs/specificity_block/offtargets.json \
    --out-dir specs/specificity_block/allbyall_inputs
# then submit rf3 folds (see the oracle note below)
```

**Templating — what to template, and with what.** The paper templates the protein
chain here and nowhere else in the pipeline: *"the most recent AF3 prediction before
the all-by-all folding was used as the template for the protein chain"*. So the
template is the design's own prior prediction (its Stage 3 on-target fold), protein
chain only — never the DNA, or the fold would be handed its own answer.

Measured on the natural-TF control panel (`analysis/oracle_controls/RESULTS.md`),
templating **raises ΔminPAE without costing discrimination**:

| TF | ΔminPAE untemplated | templated | argmin | interface still in-motif |
|---|---|---|---|---|
| LambdaRep | +1.93 | **+3.43** | held | yes |
| Engrailed | +0.26 | **+0.72** | held | yes |

It does **not** fix the absolute-minPAE offset (see the oracle note), so it is worth
doing for ranking quality, not for threshold calibration.

Mechanics for rf3 (`analysis/oracle_controls/` has working implementations):
`make_predicted_templates.py` extracts a protein-only CIF from a prior prediction;
the spec then carries a `{"path": "...cif"}` component plus a top-level
`"template_selection": ["A"]` alongside the DNA `seq` components.
`submit_templated_folds.py` handles the delivery, which is **not** obvious: a
templated rf3 input cannot be prepared directly, because pecli decides bare-fold
vs paid-MSA-pipeline by looking for an MSA carrier on a *sequence* component and a
templated input has no protein `seq` component to hang one on. The workaround is
pecli's own prepare → edit → submit flow (prepare the untemplated input to get a
bare run, swap in the templated spec plus the CIF as an `aux_files` companion, then
submit).

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

## Oracle note — use rf3

**Settled empirically.** All four open oracles were run over a natural-TF control
panel (5 sequence-specific TFs, 1 non-specific duplex binder, 2 non-binders × 8 DNA
targets). The metric of record is *argmin*: does the panel-wide lowest minPAE land on
each TF's own cognate site? See `analysis/oracle_controls/RESULTS.md`.

| oracle | argmin | binder vs non-binder gap | $/fold |
|---|---|---|---|
| **rf3** | **4/5** | **+9.12 Å** (excl. TBP) | **$0.022** |
| esmfold2 | 4/5 (double-centred) | not measured | $0.22 |
| protenix | 2/5 | −2.97 Å (**overlap**) | $0.042 |
| openfold3 | 0/5 | not measured | $0.31 |

**rf3 is primary** — jointly the most discriminative and by far the cheapest, and the
only oracle with a measured binder/non-binder separation. **esmfold2 is the
cross-check** (it ties rf3 on discrimination and is the only oracle that models TBP
at all, but costs 10×, so use it for spot checks on top-ranked designs rather than
panel-wide).

Three earlier claims in this file were wrong and are corrected:

- ~~"protenix as primary"~~ — protenix writes **no per-token PAE** unless the run
  passes `--need-atom-confidence true`, and on this panel its ranges *overlap*
  between real binders and non-binders (the non-specific binder Sac7d scores the
  lowest minPAE of any protein, 0.61 Å).
- ~~"esmfold2 cannot be used"~~ — true when written, fixed upstream by pecli #175.
  `--emit-pae true` now writes a per-token PAE.
- ~~"openfold3 as cross-check"~~ — 0/5, worse than chance.

**Absolute thresholds do not transfer.** Only 1 of 5 real TFs clears the paper's
`minPAE < 1.25` gate on rf3, and that offset survived a deep MSA (−0.11 Å), 5×
sampling (≈0) and protein templating (−0.93 to +0.11 Å) — it is oracle calibration,
not a config error. Recalibrate the gate empirically on rf3 output. The paper's
calibration-free **`ΔminPAE > 0`** criterion — which it reports as the one that
"enriched for successful designs experimentally" — does transfer, and is met by 4/5
TFs on rf3.

**Scope limit:** the metric is validated for **major-groove readers**, which is the
mode rfd3na's ori-token + H-bond conditioning targets. TBP, a minor-groove reader
that kinks DNA ~80°, fails on rf3, protenix and openfold3 alike — flagged as a stress
test in `curated_controls.json` before any fold was run.

## Reference for PRNP

Paper reported the **highest specificity-block hit rate for PRNP: 13/96 specific
designs**. DBS5 (specificity-block design for this site) is specific over 35/40
single-base variants — our SBS panel is built to reproduce exactly that test.
