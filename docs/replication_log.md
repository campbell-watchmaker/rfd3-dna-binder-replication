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

## First-pass scale & sequencing decisions

- **Smoke test first** (~10 designs) to validate the spec end-to-end before GPU budget.
- **CPU baseline analyses deferred** until the binder + specificity block architecture
  is built (user decision).

## Progress

- [x] Repo scaffolded.
- [x] PRNP-site target prepared.
- [x] Binder-block spec authored.
- [ ] Specificity-block spec authored.
- [ ] Baseline CPU analyses (DNA-similarity, TF embedding, ΔminPAE re-derivation).
- [ ] Generation run via pecli.
- [ ] Returned designs analyzed.
- [ ] Figures + public writeup.
- [ ] Reusable campaign-analysis skill.
