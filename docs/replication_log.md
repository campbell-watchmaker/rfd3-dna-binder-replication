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

## Progress

- [x] Repo scaffolded.
- [ ] PRNP-site target prepared.
- [ ] Binder-block + specificity-block specs authored.
- [ ] Baseline CPU analyses (DNA-similarity, TF embedding, ΔminPAE re-derivation).
- [ ] Generation run via pecli.
- [ ] Returned designs analyzed.
- [ ] Figures + public writeup.
- [ ] Reusable campaign-analysis skill.
