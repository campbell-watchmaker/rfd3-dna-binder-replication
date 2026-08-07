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

**RFdiffusion3 / rfd3na sampler** (see the provenance correction below; pecli `rfd3na`
defaults match):
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
  atoms; conditioning on all of them may over-constrain diffusion. Subset to the
  handful of major-groove acceptors/donors on the poly-purine core actually being read.
  Documented in PIPELINE.md. NB the parenthetical "(the paper conditions on a selected
  subset)" that used to close this line was NOT verified and is now known to be
  unsupported -- see the provenance correction below.

**Sampler config** (`sampler_config.json`): `_smoke_test` (~10 designs, first pass per
user decision) and `_full_run` (~1000 backbones/ori, paper scale) profiles. Params:
num_timesteps 200, step_scale 1.5, gamma_0 0.6 (pecli/foundry defaults; NOT paper-stated
-- see the provenance correction below).

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


## Provenance correction — H-bond conditioning and CFG (2026-08-05)

Two claims in this log and in `specs/binder_block/sampler_config.json` asserted paper
backing they do not have. Both predate the smoke test, and 20 designs were generated
under them. Recording the correction rather than quietly editing it away.

**What the paper actually says about H-bond conditioning.** One sentence in Methods:

> "Hydrogen bond conditioning [18] was applied during generation on candidate major
> groove donor and acceptor atoms (Fig. S1)."

No count. No atom names — `N7`/`O6`/`O4`/`N6`/`N4` appear nowhere in the paper. No
statement of strand, and no statement tying the selection to an ori token's 6-bp span.
Fig. S1 is in the supplement, which is unreachable (bioRxiv 403 direct, 429 through a
text proxy across repeated attempts). The `[65]` that appears mid-sentence in the
rendered text is a bibliography marker, not a count.

So **"the paper conditions on a selected subset, e.g. the N7/O6 of the central G/A
run" was never verified and is not supported.** Our 6-atom purine-strand rule is our
construction.

**The structural point, which matters more than the count.** From Results:

> "we sample a variety of placements of the protein center of mass relative to the DNA
> target using the RFD3 ori token feature, and **a diverse set of hydrogen bond (Hbond)
> condition constraints**"

The constraint set is *sampled to be diverse across designs*. A single fixed subset is
therefore the wrong model whatever its size. The sampling mechanism is not stated, so
it has to be chosen by us either way. This is consistent with the measured result that
6 vs 8 atoms was indistinguishable over 100 refolds.

**CFG has no paper basis at all.** The strings `cfg`, `cfg_scale`, `cfg_features`,
`guidance` and `classifier-free` do not occur anywhere in Sehgal et al. The
`(paper Fig. S4f)` attribution was wrong. Note also that foundry ships
`use_classifier_free_guidance: False` as the default, so having it on is a deviation,
not a match. It is kept on, re-attributed to Butcher et al. 2025's ablation — but that
ablation's DNA effect is marginal (11% → 11.3% → 12.5% H-bond satisfaction), so this is
a weakly-supported choice and a live candidate if generation underperforms.

**What the paper DOES state for the sampler**, and which we do follow: one ori token
per six consecutive base pairs, placed 3 Å toward the major groove from the 6-bp
centroid perpendicular to the helical axis; 5100 scaffolds per ori; protein length
120–150; `is_non_loopy` True; target DNA held fixed during diffusion.

**Implementation signals on magnitude** (foundry, not the paper, and flagged as such):
training subsamples H-bond atoms with the kept fraction interpolating 0.9 → 0.1 as the
true H-bond count rises to 50, i.e. roughly a third of candidates for a 36-atom target;
and the shipped 14-bp `na_binder_design.json` example specifies 16 atoms across 4 base
positions on **both** strands, mixing base *and* phosphate/sugar atoms. Both suggest our
6-atom, purine-strand-only, base-edge-only rule is narrower than upstream practice.


## Phase 1 — specificity block made runnable (2026-08-06)

The specificity block had never run. Fixing it turned up two deviations that need
recording, one recalibrated threshold, and a script that could not have worked.

### `build_fold_inputs.py` emitted a shape no oracle accepts

`specs/binder_block/PIPELINE.md` Stage 6 and `specs/specificity_block/PIPELINE.md`
Stage 3 both pointed at `scripts/build_fold_inputs.py`, which emitted
`{"id": ..., "chains": [{"id","type","sequence"}, ...]}`. **No pecli oracle accepts
that.** The smoke test worked around it with an ad-hoc emitter and never committed the
fix, so the committed pipeline could not have folded a single design.

It now emits the rf3 component shape, reads LigandMPNN `.fa` files directly, and writes
the `folds_manifest.json` that `build_filter_manifest.py` consumes. Verified by
re-emitting the smoke test's 50 Stage-6 inputs from the original LigandMPNN output: **all
50 specs and every science-bearing manifest field are byte-identical** to the ones that
actually ran.

One behaviour changed deliberately. LigandMPNN writes its input sequence as record 0 of
each `.fa`; the old `--skip-wt` dropped record 0 *positionally*, which silently deletes a
real design from any FASTA that has already been filtered. The WT record is now
identified by the **absence of `id=`** in its header, which is what actually
distinguishes it. `build_allbyall_inputs.py --skip-wt` keeps the positional behaviour and
now documents the hazard.

### Deviation: off-target panel padded to a common length

4 of the 9 Table 1 decoys are 10 bp against a 12-bp on-target. minPAE is a **minimum**
over protein×DNA token pairs, so a shorter duplex simply offers fewer pairs to minimise
over and is systematically disadvantaged as an off-target — biasing ΔminPAE upward for
reasons unrelated to specificity. Every target is now centred in the verified-neutral
flank at one length (12 bp for PRNP, so padding is 1 bp per side at most), reusing
`build_duplex()`/`verify_panel()` from `analysis/oracle_controls/control_panel.py`.
`verify_panel()` additionally refuses the panel if a padded decoy picks up another
target's motif from its flank.

Same confound the oracle control panel was padded to 24 bp to remove, where padding did
not degrade discrimination (rf3 argmin 4/5 on padded duplexes). **The paper appears to
fold Table 1 sites at native length**; `--no-pad` reproduces that. A previously emitted
`same_length_as_on` field recorded the problem and was read by nothing.

### Recalibrated: the specificity-block entry gate, minPAE < 1.25 → < 6.6

**This replaces a paper number and the derivation is here, not buried in a config.**

The paper gates entry at `minPAE < 1.25`, measured on AF3. Applied literally to rf3
output it admits almost nobody: of five real, crystallographically-characterised TFs
folded against their own cognate sites, **1/5** clears 1.25.

The offset is rf3 calibration rather than a modelling error. It survived every
intervention tried: adding deep MSAs (−0.11 Å), 5× sampling (≈0), and protein templating
(−0.93 to +0.11 Å). Measured rf3 on-target minPAE across the control panel:

| class | protein(s) | on-target minPAE |
|---|---|---|
| sequence-specific | Zif268, LambdaRep, MAX_bHLH, Engrailed, TBP | 1.07 – 3.88 Å |
| non-specific duplex binder | Sac7d | 6.60 Å |
| non-binders | Ubiquitin, GFP | 15.7 – 16.2 Å |

A cut at **6.6 Å** is where binders separate from non-binders on measured data.
**Provisional: it rests on 5 proteins.** The *success* criterion downstream stays the
paper's own calibration-free **`ΔminPAE > 0`** — the criterion it reports as the one that
"enriched for successful designs experimentally", and the one control-panel result that
did transfer to rf3 (4/5 TFs). A sign test needs no recalibration; an absolute cut does.

### minPAE is now emitted by the binder block

The entry gate above was unevaluable: `filter_binder_block.py` emitted RMSD/ipTM/H-bonds
only, though the PAE was already on disk in each rf3 refold's `*_confidences.json`. It
now reports `min_pae` per design and gates on it with `--min-pae-gate`, reusing
`compute_delta_minpae.py`'s loader and masks so the number is defined identically in both
places. A requested gate **fails** rows with no readable PAE rather than passing them.

Measured on the smoke test's 50 real refolds: minPAE recovered 50/50, min 2.84, median
6.23, max 13.01. Note the glob must exclude `*_summary_confidences.json`, which holds
scalars only and no matrix.

### Two silent failures in the ΔminPAE ranking

`compute_delta_minpae.py` had a **guaranteed crash**: `tuple(j["protein_len"])` against
the int that `build_allbyall_inputs.py` emits. Python evaluates call arguments eagerly,
so it fired on record 1 of every real manifest even when chain labels were present and
the range was never consulted. Fixed in the consumer, not the emitters, because the
scalar form is the repo's majority convention. It went unseen because every existing test
hand-built a manifest that omitted the field; there is now an integration test that feeds
the real emitter's output straight into the consumer.

Worse, it **silently dropped** any design whose on-target fold was missing or
uncollected — no warning, no count, no row — so a partially drained batch produced a
shorter, entirely credible, wrong CSV with nothing downstream able to tell a design that
ranked badly from one that was never scored. It now reports them, and warns rather than
picks arbitrarily when a design has duplicate on-target rows.

### Chain layout, measured not assumed

`chains_to_design` was the literal placeholder `"SET_FROM_RFD3NA_OUTPUT"` in both
LigandMPNN configs. Read off a real relaxed backbone it is **`C`**: rfd3na emits the
fixed target duplex first, so the relaxed PDB is DNA A(12 nt) + DNA B(12 nt) + protein
C(123 aa). This is the *opposite* of the refold layout (protein A + DNA B,C), and that
mismatch is what caused the Stage-7 bug where protein Cα found zero overlap and the DNA
aligned on the wrong strand.

### Also corrected

- `specs/specificity_block/fold_config.json` still specified protenix and the SBS panel;
  `specs/binder_block/fold_config.json` still specified a three-oracle comparison and the
  dead complex-JSON shape. Both rewritten to the measured rf3 configuration, with the
  superseded content noted rather than deleted.
- `PIPELINE.md` showed `--iptm-gate 0.9` without `--stage post_resample`, which silently
  disables the ipTM gate.
- Binder-block Stage 0 said `pecli prepare protenix` for the DNA-only duplex; a
  protein-free protenix fold is impossible in pecli (issue #189). Now rf3.
- `submit_templated_folds.py` indexed `rec["protein"]`, a KeyError on any specificity
  manifest (which carries `design_id`); it now accepts either, prefers the manifest's
  recorded `template_cif` over reconstructing the filename, and takes `--group` /
  `--description` so it cannot mix two campaigns' runs into the control panel's group.
- `make_predicted_templates.py` gained a `--design-fasta` mode (`expect_chains=1`,
  residue count from the sequence) sharing `extract()` verbatim, so the design path keeps
  the zero-nucleotide and residue-count gates that the crystal-template attempt paid for.

Tests: 53 → 79.


## Phase 2 (partial) — H-bond conditioning is now sampled (2026-08-06)

`make_rfd3na_specs.py --hbond-sampling random` draws an independent constraint set per
design instead of reusing one fixed subset per ori. This closes the shape mismatch the
2026-08-05 provenance correction identified: the paper varies the constraint set across
designs, so a single fixed subset was wrong regardless of its size — which is also why
the 6-vs-8-atom comparison came out indistinguishable over 100 refolds. No fixed count
is the fix.

**What is sampled, and what is ours.** The paper states the diversity and not the
method, so the mechanism is our construction and is labelled as such in the code:

- base-position count, uniform over `--hbond-bases-min/max` (default 2–5, spanning 4–10
  atoms against PRNP's 36 candidates — which straddles foundry's ~one-third training
  subsample rate, an implementation signal, not paper evidence);
- which positions, uniform without replacement inside that ori's own bp window;
- which strand, drawn per design under `--hbond-strand either` (previously purine-only;
  foundry's shipped 14-bp example uses both strands).

Two invariants are preserved from the fixed rule and tested: both atoms of a purine are
always kept together (G N7+O6, A N7+N6 are the bidentate pairs Arg and Asn form against
a purine, so splitting them specifies a weaker and less physical constraint), and the
duplex's terminal base pairs are never drawn.

Every draw is recorded per spec in `manifest.json` under `hbond_draw`, and the manifest
carries the seed, so a design traces back to the exact constraint set that produced it.
`--designs-per-ori > 1` without sampling is refused rather than silently emitting N
identical specs. `fixed` remains the default so the existing arms stay reproducible.

Also from Phase 2: the specificity entry gate recalibration (minPAE < 6.6) and minPAE
emission from the binder block landed with Phase 1 above, because Stage 0 was
unexecutable without them.

Tests: 79 → 94.


## Stop point — 2026-08-06

Working tree clean, 94 tests passing, two commits on
`analysis/delta-pae-oracle-controls`. Everything below is CPU-only work; **no GPU spend
this session**.

### Done

Phase 1 complete except its validation run (below). Phase 2: sampling, the recalibrated
entry gate, and minPAE emission are in; the two A/B arms (~$1.3) are not.

### Blocked, and the exact blocker

**The control-panel validation cannot run until the rf3 PAE files are back on disk.**
The scratchpad holding them was cleaned. Re-downloading is free — the 64 runs are still
`SUCCEEDED` — but the route is not the obvious one:

1. **The `run_id`s committed in `analysis/oracle_controls/folds/rf3_manifest.json` do not
   resolve.** `pecli status 5347ff` → `unknown run: 5347ff`. These are the IDs recovered
   by regex from the submitter log after the manifest race, and at least some are wrong.
   Do not trust that file's `run_id` field.
2. `pecli runs --group oracle-controls-rf3` lists **63** SUCCEEDED rf3 runs (not 64 —
   worth resolving), but the table **truncates the description**, so the fold_id cannot
   be read from it. Real IDs: `analysis/oracle_controls/folds/` has none usable; a
   scraped list is not committed.
3. The mapping route that does work: `pecli results <id> --out <dir>` for each of the 63,
   then map each directory to its fold by the emitted `<fold_id>_confidences.json`
   filename. ~14 s per `pecli status` call, so avoid status entirely and download
   directly.
4. `collect_control_results.py --skip-status` is a trap here: it **refuses to download
   anything not already present** (it only assumes-OK what is on disk). Run it without
   the flag, or drive `pecli results` directly.

Then: `analysis/oracle_controls/validate_specificity_block.py --control-manifest <that>
--out-dir <scratch>` reshapes the panel into the specificity block's manifest schema,
runs `scripts/compute_delta_minpae.py` over it, and asserts the recorded answer (argmin
4/5, ΔminPAE > 0 for 4/5). It exits non-zero on disagreement. Its reshape logic is
already unit-tested; only the end-to-end run is outstanding.

### Next, in order

1. Restore the 63–64 rf3 control PAEs (free) and run `validate_specificity_block.py`.
   **Nothing should run on designs until this reproduces** — it is the only place where
   the block's answer is known independently.
2. Phase 2.5: the two A/B arms, 20 backbones each, 1 seq/backbone, rf3 — CFG on/off and
   sampled vs fixed conditioning, ~$0.64 each. Compare with
   `scripts/compare_conditioning_arms.py`. n=20 resolves only large effects.
3. Phase 3 paper scale (~$367) — gated on 1 and 2 being green.

### Carried forward, unresolved

- **Placement, not foldability, is the failure mode** (9/50 folded to <3 Å protein-only
  but sat >8 Å off after DNA alignment; median DNA-aligned RMSD ~30 Å). If Phase 2's arms
  leave that unchanged, the next suspects are ori-token placement and the relax's DNA
  restraint — not the H-bond conditioning.
- The recalibrated 6.6 Å gate rests on 5 proteins.
- CFG is on with no paper basis (foundry ships `False`).
- The H-bond sampling *mechanism* is ours; the paper states only that the set is diverse.


## Specificity block VALIDATED on the control panel (2026-08-06)

The block's own ranking code now reproduces the control panel's recorded answer exactly.
This was the gate on running anything on designs, and it is green.

`analysis/oracle_controls/validate_specificity_block.py` reshapes the 64-fold natural-TF
panel into the specificity block's manifest schema and runs
**`scripts/compute_delta_minpae.py`** over it — the block's own script, not the
control-panel one that produced the original numbers. Two implementations of a minimum
over the same matrix is exactly where a silent inconsistency would live.

| protein | ΔminPAE | on-target | best off-target | argmin |
|---|---|---|---|---|
| LambdaRep | **+1.93** | 3.87 | 5.80 | lambda_OL1 ✓ |
| MAX_bHLH | **+1.83** | 2.15 | 3.98 | ebox ✓ |
| Zif268 | **+1.60** | 1.07 | 2.67 | zif268_site ✓ |
| Engrailed | **+0.26** | 3.88 | 4.14 | hd_taatta ✓ |
| TBP | −2.16 | 15.73 | 13.57 | scramble ✗ |

argmin on own site **4/5**, ΔminPAE > 0 **4/5** — matching the recorded result value for
value, TBP included (the minor-groove reader flagged as a stress test before any fold was
run). Outputs committed to `results/specificity_block/validation_controls_*.csv`.

### Recovering the PAE files: a pecli resolver bug, not our data

The scratchpad holding the PAEs had been cleaned, and re-downloading turned out to be the
hard part. `pecli runs --group oracle-controls-rf3` lists all 63 runs as SUCCEEDED, but
`pecli results <id>` answers **`unknown run`** for every one of them.

The cause is in pecli (`pecli/runner.py::_resolve`), not in our ids:

    find_by_short_id -> runs.list_for_user(user_id, limit=100)   # 100 NEWEST, per user
    fallback         -> runs.list_all(limit=200)                 # 200 NEWEST, team-wide

Both are **recency windows**; `pecli runs --group` queries the sparse group index, which
has none. So any run older than roughly the 100 most recent is listed but unresolvable,
and the two commands disagree about whether it exists. These folds are from 2026-07-30.
Worth filing: `pecli results` should fall back to the same index `pecli runs` uses rather
than calling a listed run unknown.

**This retires the earlier diagnosis.** The stop-point note said the run_ids committed in
`folds/rf3_manifest.json` were corrupted by the manifest race. They are not — `5347ff`,
the id cited there as proof, is correct and resolves to the right fold. Neither those
ids, nor ids scraped from `pecli runs`, nor the full `YYYYMMDD-HHMMSS-xxxxxx` ids from
`~/.pecli/jobs.json` resolve, because the id was never the problem. (`pecli runs` also
displays only the 6-char suffix of a run id, which sent the first attempt down the wrong
path.) A `remap_run_ids.py` written against that wrong diagnosis was deleted rather than
kept.

`analysis/oracle_controls/fetch_group_results.py` is the working route: query the group
index, read `s3_output` off each Run record, and call the same S3 download `pecli results`
would have, skipping resolution entirely. Free — S3 GETs only. It recovered 63/63.

The 64th fold, **Zif268's on-target**, is not in the rf3 group at all: it was submitted
earlier as the schema probe, under `oracle-controls-probe`. Nothing was wrong with it —
but note the panel's most load-bearing single fold sits outside the group its 63 siblings
share, which is worth knowing before anyone treats the group as the panel.


## Phase 2.5 — the A/B arms (2026-08-06)

**Design: a 2×2 factorial, not two separate A/B pairs.** CFG (on/off) × H-bond
conditioning (fixed/sampled), 20 designs per arm, 80 total. Same cost as the planned two
pairs (~$0.96 diffusion + ~$1.60 refold), but each main effect is measured at n=40 per
level instead of n=20. Both factors are crossed, so an interaction is at least visible
rather than aliased.

Arm structure differs by necessity: the fixed arms are 2 specs (one per ori) × batch 5 ×
2 batches; the sampled arms are 10 specs (5 independent draws × 2 ori) × batch 2, because
a sampled arm needs a *different constraint set per spec* — that is the thing being
tested. Seed 2025, and every draw is recorded in each arm's `manifest.json`.

All 24 diffusion runs SUCCEEDED, ~$0.96.

**Caveat carried from `compare_conditioning_arms.py`, and it still applies:** rfd3na
exposes no seed control, so the arms are independent draws rather than paired. n=20 per
arm resolves only large effects. A shifted median with overlapping ranges is noise.

### Two documentation bugs found by trying to run the documented command

`specs/binder_block/PIPELINE.md` Stage 3 showed

    pecli prepare rfd3na --design-inputs "$spec" --config sampler_config.json:_smoke_test

**Neither flag exists.** pecli answers `unknown option(s) for rfd3na: config,
design_inputs`. The spec goes in via `--input`, and every sampler knob is its own flag.
The command as documented could never have run. It also omitted the target-CIF staging
step, without which a run starts and dies immediately on a CIF that was never uploaded.
`scripts/submit_arm_diffusion.py` now drives prepare → stage → submit as one sequence.

### `cfg_features` has never been applied

`sampler_config.json` carries `cfg_features: [active_donor, active_acceptor]` with a
paragraph of rationale. **pecli's rfd3na tool exposes no such field**
(`pecli/tools/rfd3na.py`), and `pecli prepare` rejects unknown options — so the setting
has never reached the sampler on any run this project has made, including the 20 designs
of the original smoke test. `use_classifier_free_guidance` and `cfg_scale` are real and
do apply, so the CFG arms are still a valid comparison; what was never true is the claim
that we guide on those two features specifically. Annotated in the config rather than
deleted, and the submitter reports it as skipped. A test now translates the committed
config so any future inert knob fails CI.

### A second pecli group-label trap

pecli slugifies a `--group` label at submit time, so `p25-fixed_cfgon` is stored as
`p25-fixed-cfgon`. The CLI slugifies the query too, so `pecli runs --group` matches
either spelling — but querying the group index directly does not, and returns a
confident, wrong `0 run(s)`. `fetch_group_results.py` now slugifies first. Related to,
but distinct from, the resolver bug fixed in pecli PR #191.

### Phase 2.5 results (2026-08-06)

80 designs, 4 cells of 20, all stages clean: diffusion 24/24, relax 80/80, LigandMPNN
80/80, rf3 refold 80/80. Spend ~$2.56 (diffusion $0.96, refold $1.60; relax and
LigandMPNN free). Per-arm CSVs in `results/phase25/`, and the sampled arms' exact
constraint draws in `results/phase25/sampled_hbond_draws.json`.

**Main effect: CFG** (n=40 per level)

| metric | CFG on (min/med/max) | CFG off (min/med/max) |
|---|---|---|
| protein-only Cα-RMSD | 1.89 / **6.99** / 19.92 | 2.29 / **6.10** / 20.34 |
| DNA-aligned Cα-RMSD | 16.31 / **32.04** / 73.95 | 10.45 / **30.91** / 41.02 |
| major-groove H-bonds | 0 / **2** / 11 | 0 / **2** / 9 |
| ipTM | 0.26 / **0.54** / 0.72 | 0.28 / **0.52** / 0.76 |
| passes 8 Å gate | 0/40 | 0/40 |
| folds correctly (<3 Å protein-only) | 5/40 | 6/40 |
| reads bases (≥3 major-groove H-b) | 15/40 | 16/40 |

**Main effect: H-bond conditioning** (n=40 per level)

| metric | fixed (min/med/max) | sampled (min/med/max) |
|---|---|---|
| protein-only Cα-RMSD | 2.44 / **6.91** / 19.92 | 1.89 / **6.38** / 20.34 |
| DNA-aligned Cα-RMSD | 12.26 / **32.14** / 73.95 | 10.45 / **31.75** / 52.81 |
| major-groove H-bonds | 0 / **2** / 9 | 0 / **2** / 11 |
| ipTM | 0.28 / **0.51** / 0.76 | 0.26 / **0.54** / 0.72 |
| passes 8 Å gate | 0/40 | 0/40 |
| folds correctly | 6/40 | 5/40 |
| reads bases | 15/40 | 16/40 |

**Four cells** (n=20): passes-8 Å 0/0/0/0; folds-correctly 2/4/3/2; reads-bases 8/7/7/9;
median DNA-aligned RMSD 32.47 / 31.80 / 32.04 / 30.17.

Every median difference sits well inside the overlapping ranges, and n=20 per cell with
unpaired arms (rfd3na exposes no seed control) resolves only large effects. The
comparison script prints that caveat with every table.

**The fixed/CFG-on cell is the smoke test's own configuration**, and it reproduces the
smoke test closely — median DNA-aligned RMSD 32.47 Å vs ~30 Å, median major-groove
H-bonds 2 vs 2. That is a consistency check on the rebuilt pipeline: the Stage-6 emitter
was rewritten, the manifest schema changed, and minPAE was added since those numbers
were taken.

**0/80 designs clear the paper's 8 Å pre-resample gate**, in every cell. Median
DNA-aligned RMSD ~31 Å against a median protein-only RMSD of ~6.4 Å. As the smoke test
found, the designs fold roughly as intended but are not placed on the duplex; neither
factor tested here moved that. Per the plan's own risk note, the remaining suspects are
ori-token placement and the relax's DNA restraint rather than the H-bond conditioning.

### Correction — `cfg_features` was never absent, only never ours (2026-08-07)

The Phase-2.5 note above says `cfg_features` "has never reached the sampler on any run".
That is right about *our value* and **wrong about the effect**, and the difference
matters. Read from foundry source (`models/rfd3na/`):

- pecli's rfd3na tool exposes no `cfg_features` field, and `containers/rfd3na/pipeline.py`
  overrides only `use_classifier_free_guidance` and `cfg_scale`. So our list is inert.
- But the field then falls back to foundry's own default in
  `configs/inference_engine/rfdiffusion3.yaml`:
  **`[active_donor, active_acceptor, ref_atomwise_rasa, bp_partners]`**
  (the checkpoint sampler config `configs/model/samplers/edm.yaml` carries the first
  three).

So **both H-bond features we wanted guided ARE guided**. What we additionally get,
unasked, is guidance on `ref_atomwise_rasa` (per-atom relative solvent accessibility) and
`bp_partners` (base-pair partners). Setting our two-item list is a no-op; *narrowing* to
only the H-bond features would need a pecli change.

**What CFG actually does here** (`model/inference_sampler.py:273-305`, `model/cfg_utils.py`):
at each diffusion step the model is run twice — once with the conditioning features, once
with the `cfg_features` entries **zeroed** (`strip_f`) — and the two update directions are
combined as

    delta = delta_cond + (cfg_scale - 1) * (delta_cond - delta_uncond)

i.e. it extrapolates *away* from what the model would have done knowing nothing about the
H-bond constraints. `cfg_scale = 1.0` is a no-op; our 1.5 pushes 50% past the conditional
prediction. `cfg_t_max` is null, so it applies at every timestep. **Cost: two forward
passes per step instead of one** whenever CFG is on.
