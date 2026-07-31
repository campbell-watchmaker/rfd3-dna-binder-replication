# Replication log

## Scope

Replicate the Sehgal et al. 2026 DNA-binder pipeline end-to-end against a single target, the **PRNP-site**
(`TGAGGAGAGGAG`, target T1 in the paper's Table 1). In-silico only.

## Target rationale

The PRNP-site is the paper's best-characterized target: it reported the highest specificity-block hit rate
(13/96 specific designs for this site) and its strongest-affinity characterized binders bind here (DBB5 at
3 nM, DBB3 at 10 nM), recognizing the poly-purine tract via Asn/Arg major-groove contacts. That gives us
concrete reference designs to benchmark our returned designs against.

## Key parameters (from the papers, to hold fixed)

**RFdiffusion3 / rfd3na sampler** (paper Fig. S4f; pecli `rfd3na` defaults match):
- protein length 120–150
- step_scale (η) = 1.5, num_timesteps = 200, gamma_0 (γ₀) = 0.6
- classifier-free guidance available (cfg_scale); DNA held fixed during diffusion
- `is_non_loopy = True`
- ori (center-of-mass) tokens: one per 6 consecutive bp, placed 3 Å toward the major groove from the
  stretch centroid, perpendicular to the helical axis
- H-bond conditioning on candidate major-groove donor/acceptor atoms

**LigandMPNN:** temperature 0.1; 5 seq/backbone (binder block), 100 seq/backbone (specificity resample);
the paper relaxes the rfd3na output with Rosetta FastRelax before sampling.

**Relaxation substitution (open-source requirement):** Rosetta is free for academic use but
is not permissively licensed, so this replication uses **OpenMM** (MIT/LGPL) instead. The
diffused protein–DNA complex is energy-minimized with an Amber ff14SB (protein) +
OL15/bsc1 (DNA) force field combination, with **DNA atoms under a positional restraint** and
the protein free to relax — consistent with rfd3na treating the DNA as fixed throughout
diffusion. This is the same class of step AlphaFold2's Amber-relax post-processing performs
(clash/stereochemistry cleanup after generation), applied here to the rfd3na output instead
of Rosetta FastRelax. Runs CPU-side (`scripts/relax_openmm.py`); no GPU hop needed for a
single structure. See README "Substitutions vs. the original" for the full list of
open-source swaps.

Note: pecli's own `gromacs` tool was considered and rejected for this step — it is scoped to
protein-only PDBs and rejects nucleic acids, ligands, and metals at prepare time (see
pecli ADR 0050), so it cannot see or restrain the DNA half of the complex being relaxed here.

**Binder-block filters:** DNA-aligned protein Cα-RMSD < 8 Å → resample → < 3 Å, ipTM > 0.7, high H-bond counts.
**Specificity-block filters:** binder-block passers with minPAE < 1.25 → resample (100) → < 1.5 Å RMSD,
ipTM > 0.9 → templated all-by-all fold → rank by ΔminPAE, take top 96.

**ΔminPAE** = min over off-targets of (minPAE_offtarget) − minPAE_ontarget, where
minPAE = min over protein–DNA residue pairs of PAE(i, j).

> **Corrections from a close read of the paper's Methods (2026-07-31).** Three
> things our specs get wrong or omit:
>
> 1. **Templating.** Templates are used nowhere in the paper's pipeline *except*
>    the specificity block's all-by-all fold, where "the most recent AF3
>    prediction before the all-by-all folding was used as the template for the
>    protein chain" (a self-template). Their native-TF minPAE benchmark is
>    likewise "with the protein templated and run in single-sequence mode". Our
>    specificity spec does not template. This matters: templating removes protein
>    fold uncertainty so minPAE reflects the interface. See
>    `analysis/oracle_controls/RESULTS.md`.
> 2. **Off-target panel is over-built.** The paper's ΔminPAE all-by-all runs
>    against the on-target plus the *other Table 1 targets* (6 core, +10
>    additional where applicable) — **not** single-base variants. The "specific
>    over 35/40 single-base variants" claim is separate wet-lab characterisation
>    of one binder (DBS5), not the ranking panel. `scripts/make_offtarget_set.py`
>    builds 46 targets including 36 single-base substitutions, which inflates the
>    all-by-all ~3× and changes what ΔminPAE means (a minimum taken over
>    near-identical variants is a much harsher denominator than one taken over
>    unrelated sites).
> 3. **Binder block has an earlier gate we omit.** The sequence is: fold →
>    **DNA-aligned RMSD < 8 Å** → LigandMPNN resample → fold → RMSD < 3 Å,
>    ipTM > 0.7, high H-bond counts. Our spec starts at the 3 Å gate. (The paper
>    also states `ΔminPAE > 0` "enriched for successful designs experimentally",
>    which is a usable criterion that needs no absolute calibration.)
>
> Not acted on: the extracted text renders the metric as "Cε-RMSD" throughout.
> Cε is not a backbone atom, so this is almost certainly a text-extraction
> artifact of "Cα"; our Cα implementation stands. AF3 seeds/samples/recycles for
> the *design* folds are not reported anywhere in either paper — only the DNA-only
> starting duplex is specified (seed 42, single diffusion sample).

**Interaction counting** (paper used DSSR v1.7.8): total protein–DNA H-bonds, major-groove H-bonds,
and "supporting" (buttressing) intra-protein H-bonds to DNA-contacting residues. Native reference =
357 JASPAR TF–DNA PDB structures with info content > 1.5.

## Released assets (from the paper)

- RFD3 DNA checkpoint: `https://files.ipd.uw.edu/pub/dna_binder_rfd3/rfd3-1030-foundry.ckpt`
- Design summary metrics: `https://files.ipd.uw.edu/pub/dna_binder_rfd3/summary_data.csv`

## Division of labour

- **Claude Science (CPU):** target prep, pipeline specs, all downstream analysis, figures, this repo.
- **pecli + Claude Code (GPU/AWS):** rfd3na generation, ligandmpnn, folding. Prepare→approve→submit gate.

## Binder-block architecture (specs/binder_block/)

Authored the binder-block pipeline: rfd3na → OpenMM relax → ligandmpnn → three-oracle
refold → filter. See `specs/binder_block/PIPELINE.md` for the stage-by-stage runbook.

**rfd3na input schema — verified against the upstream foundry reference**
(rosettacommons.github.io/foundry/models/rfd3/input.html + NA binder tutorial), not
assumed. Findings that shaped the spec generator (`scripts/make_rfd3na_specs.py`):

- `ori_token` is a **single `[x,y,z]`** per spec (COM-placement override), not a list.
  The paper's "~5100 scaffolds per ori" therefore means **one diffusion run per ori
  placement**, swept over positions. The generator emits one spec per ori (2 for the
  12-bp PRNP target) + a manifest.
- H-bond conditioning uses two `InputSelection` dicts — `select_hbond_donor` /
  `select_hbond_acceptor` — keyed by DNA residue id (`"A6"`, `"B13-24"`) with
  comma-joined atom-name strings (`"N7,O6"`). Requires **HBPLUS** installed GPU-side.
- DNA is fixed via `select_fixed_atoms: {"<dna range>": "ALL"}`; `contig` lists the
  fixed DNA chains + the designed protein length via the InputSelection mini-language.
- CFG: `use_classifier_free_guidance` + `cfg_features` (subset of `active_donor`,
  `active_acceptor`, `ref_atomwise_rasa`) + `cfg_scale` (default 1.5).
- **Caveat to apply before submit:** the generator emits *all* candidate major-groove
  atoms; conditioning on all of them over-constrains diffusion. Subset to the handful
  of major-groove acceptors/donors on the poly-purine core actually being read (the
  paper conditions on a selected subset). Documented in PIPELINE.md.

**Sampler config** (`sampler_config.json`): `_smoke_test` (~10 designs, first pass per
user decision) and `_full_run` (~1000 backbones/ori, paper scale) profiles. Params:
num_timesteps 200, step_scale 1.5, gamma_0 0.6 (paper Fig. S4f; pecli rfd3na defaults).

**Refold oracle: three-way comparison** (user decision) — protenix + openfold3 +
esmfold2 on the same designs, comparing fold quality (DNA-aligned RMSD, ipTM) AND
runtime/cost. esmfold2 needs both DNA strands listed explicitly (no auto-complement);
`scripts/build_fold_inputs.py` writes both strands so one input serves all three.

**Filtering** (`scripts/filter_binder_block.py`, CPU, here): DNA-aligned protein
Cα-RMSD (superpose refold onto design by DNA atoms, measure protein Cα displacement —
the paper's self-consistency metric), ipTM (from oracle output), and protein–DNA
H-bond counts (open geometric reimplementation replacing DSSR). Gates: RMSD < 3 Å,
ipTM > 0.7. Validated on a real complex (λ repressor–operator, PDB 1LMB): identity
pair → 0.0 Å RMSD; a 3°-rotated protein → 1.78 Å; 15–16 interface H-bonds (4
major-groove), consistent with a HTH major-groove reader. Unit-tested in
`tests/test_binder_block.py`.

## Specificity-block architecture (specs/specificity_block/)

Authored the specificity block — the ΔminPAE negative-design half that took the
paper from ~0.5% (binder block) to ~3% specific designs. See
`specs/specificity_block/PIPELINE.md`.

Flow: binder-block passers with on-target **minPAE < 1.25** → LigandMPNN resample
(100 seq/backbone, temp 0.1) → on-target pre-filter fold (RMSD < 1.5 Å, ipTM >
0.9) → **templated all-by-all** fold vs on-target + off-targets → rank by ΔminPAE,
top 96/target.

**Off-target panel** (`scripts/make_offtarget_set.py`): for PRNP, 46 targets —
on-target (ΔminPAE reference) + 36 single-base-substitution variants (3 × 12 bp;
reproduces the paper's "specific over 35/40 single-base variants" test for DBS5)
+ 9 unrelated Table 1 decoys. Decoy sequences transcribed from the paper's Table 1
(flagged in-code to verify before a production run).

**ΔminPAE** (`scripts/compute_delta_minpae.py`, CPU, here): minPAE = min over
protein-residue × DNA-residue pairs of PAE(i,j), checked in **both** PAE
orientations; ΔminPAE = min over off-targets of minPAE(off) − minPAE(on). Ranked
descending. Validated on synthetic PAE matrices: a specific design (on-target low,
off-targets high) ranks above a promiscuous one (off-target also low), and minPAE
correctly takes the global protein–DNA block minimum.

**Oracle constraint (important):** the specificity block **cannot use esmfold2** —
ΔminPAE needs a PAE matrix, which only the AF3-class folders emit. This differs
from the binder block's three-oracle comparison (which only needs RMSD/ipTM).

> **Corrected 2026-07-30 (was: "protenix is primary, openfold3 the cross-check").**
> That plan assumed protenix returns a per-token PAE matrix. It does not, as
> wrapped by pecli: a completed run retains only
> `*_summary_confidence_sample_0.json` — scalars plus 2×2 `chain_pair_*`
> aggregates — and nothing else is even written to S3. The full matrix requires
> `--need-atom-confidence true`, which additionally emits
> `*_full_data_sample_<rank>.json` with `token_pair_pae`; the array is always
> computed in memory but discarded otherwise.
>
> Meanwhile **rf3** (RosettaFold3) emits a full PAE *natively* —
> `*_confidences.json` with `pae [N,N]`, `token_chain_ids`, `token_res_ids`,
> already in the shape `scripts/compute_delta_minpae.py` parses — at roughly half
> protenix's realised cost (~$0.06 vs ~$0.12 per fold, pecli's own figures over
> ~90 runs each).
>
> **So rf3 becomes the primary specificity oracle, with protenix (PAE flag on) as
> the cross-check.** Note protenix's PAE carries no chain labels, only integer
> `token_asym_id`, so protein-vs-DNA tokens must be resolved positionally from
> input entity order and the resulting counts asserted against the submitted
> sequence lengths. See `analysis/oracle_controls/`.

`scripts/build_allbyall_inputs.py` builds one complex-JSON per (design × DNA
target) plus a `folds_manifest.json` skeleton for `compute_delta_minpae.py`.
Unit-tested in `tests/test_specificity_block.py` (3 tests).

## First-pass scale & sequencing decisions

- **Smoke test first** (~10 designs) to validate the spec end-to-end before GPU budget.
- **Analysis sequencing (revised).** Of the three planned analyses, only the
  DNA-similarity premise is a genuine *pre-generation* baseline (it depends on B-DNA
  geometry alone, not on any design) — done, PR #6. The other two are really
  *post-generation* design analyses and are deferred until designs have been returned
  from the binder + specificity blocks:
  - **TF sequence-space embedding map** — its scientific payload is whether *our*
    designs land in novel regions of DNA-binder sequence space, which requires the
    returned designs. The natural-set backdrop (Evo-1 on JASPAR TFs, ESM-2 on PDB
    complex chains) is design-independent and will be batched into the generation GPU
    session so that, once designs return, only the ~73 designs need embedding.
  - **ΔminPAE re-derivation from released data** — an independent check of
    `scripts/compute_delta_minpae.py` against the paper's released `summary_data.csv`.
    Runs on CPU with no designs, but grouped with the post-generation analysis phase
    so the specificity metric is validated right before it is applied to our designs.

## Progress

- [x] Repo scaffolded.
- [x] PRNP-site target prepared.
- [x] Binder-block spec authored.
- [x] Specificity-block spec authored.
- [x] Pre-generation analysis: DNA-similarity premise (analysis/dna_similarity/, PR #6).
- [x] Pre-generation analysis: **ΔminPAE oracle controls** (analysis/oracle_controls/).
      128 folds of 8 natural controls (5 specific TFs / 1 non-specific duplex binder /
      2 non-binders) × 8 DNA targets × {rf3, protenix}, MSA-free. $4.06, 0 failures.
      On rf3 the classes separate as designed (argmin on the correct cognate site for
      4/5 TFs; +9.1 Å binder/non-binder gap excluding TBP); on protenix they do not
      (2/5; ranges overlap). Established rf3 as the primary specificity oracle and
      measured the real per-fold cost. See analysis/oracle_controls/RESULTS.md.
- [x] Off-target decoy panel verified against Sehgal et al. Table 1 (decoy-controlled;
      note GGGCTTGCGA is labelled both Oct4-gRNA2 and Dux4-gRNA2 in the paper).
- [ ] Generation run via pecli (binder block → specificity block).
- [ ] Post-generation analysis: returned designs (DNA-aligned RMSD, ipTM, interactions).
- [ ] Post-generation analysis: ΔminPAE re-derivation from released data (analysis/delta_minpae/) — validates the metric before applying it to our designs.
- [ ] Post-generation analysis: TF sequence-space embedding map (analysis/tf_embedding/) — natural backdrop batched into the generation GPU session; designs overlaid after they return.
- [ ] Figures + public writeup.
- [ ] Reusable campaign-analysis skill.
