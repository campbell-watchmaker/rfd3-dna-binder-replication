# De novo sequence-specific DNA-binder design — a replication

An open, end-to-end replication of the DNA-binding-protein design pipeline from
**Sehgal et al. 2026, *Generative design of sequence specific DNA binding proteins*** (bioRxiv `2026.04.27.720408`),
which builds on **Butcher et al. 2025, *De novo Design of All-atom Biomolecular Interactions with RFdiffusion3*** (bioRxiv `2025.09.18.676967`).

This repository reproduces the pipeline against one target — the **PRNP-site** (`TGAGGAGAGGAG`), the poly-purine
core target for which the original paper reported its best-characterized binders (e.g. DBB5, KD ≈ 3 nM).

> **Status:** work in progress. This is an *in-silico* replication: it reproduces the computational
> design and analysis, not the wet-lab characterization (yeast display / BLI / HT-SELEX) of the original.

## The method in one paragraph

DNA is a hard de-novo target because B-form DNA has nearly the same global shape for every sequence, and
base-to-base chemical differences are subtle — so designing for affinity alone yields promiscuous binders.
The pipeline addresses this in two stages. A **binder block** uses RFdiffusion3 (the all-atom successor to
RFdiffusion) to generate a protein wrapped around a fixed DNA target, conditioned on center-of-mass ("ori")
tokens and major-groove hydrogen-bond constraints, then designs its sequence with LigandMPNN and filters on
an AlphaFold3-class refold. A **specificity block** then does explicit negative design: each candidate is
folded against every off-target DNA and ranked by **ΔminPAE** = min(off-target minPAE) − on-target minPAE.
Testing only 96 designs per target, the original reported specific binders for 7 targets — roughly a 100×
improvement in specific-design success rate over prior scaffold-docking approaches.

## Pipeline

```
target DNA ──▶ [ AF3-class fold → B-DNA duplex ]
                        │  + ori tokens (1 / 6 bp, 3Å into major groove)
                        │  + major-groove H-bond donor/acceptor atoms
                        ▼
   BINDER BLOCK    rfd3na  ──▶  ligandmpnn (T=0.1)  ──▶  fold (protenix/openfold3)
                        │        └─ filter: DNA-aligned RMSD, ipTM, H-bond counts
                        ▼
   SPECIFICITY    ligandmpnn resample (100/backbone) ──▶ all-by-all fold (on + off targets)
   BLOCK                └─ ΔminPAE ranking ──▶ top 96 designs
```

Generation runs on GPUs via **[pecli](https://github.com/watchmaker-genomics/pecli)** (an MCP/CLI layer over
on-demand AWS GPUs); `rfd3na` is the nucleic-acid-capable RFdiffusion3. All CPU-side science —
target preparation, interaction analysis, the specificity metric, embeddings, and figures — is reproducible
from this repo.

## Layout

| Path | Contents |
|---|---|
| `targets/prnp/` | PRNP-site duplex structure, ori-token placements, H-bond candidate atoms |
| `specs/` | pecli pipeline specs for the binder block and specificity block |
| `analysis/dna_similarity/` | "why specificity is hard" — B-DNA cross-sequence similarity |
| `analysis/tf_embedding/` | TF binding-sequence-space map (Evo / ESM-2 + PCA) |
| `analysis/delta_minpae/` | ΔminPAE re-derivation from the released summary data |
| `results/` | returned designs + per-design metrics |
| `figures/` | publication-grade figures |
| `docs/` | method notes, replication log, limitations |

## Substitutions vs. the original (stated honestly)

| Original | Here | Why |
|---|---|---|
| AlphaFold3 (filtering + specificity oracle) | protenix / openfold3 (open AF3-class) | AF3 weights are license-restricted |
| DSSR v1.7.8 (interaction counting) | open reimplementation (Biotite-based H-bond geometry) | DSSR is not freely redistributable |
| Rosetta FastRelax (pre-LigandMPNN relax) | OpenMM energy minimization, DNA restrained | Rosetta is free-for-academic but not permissively licensed |
| Wet-lab validation (yeast display, BLI, SELEX) | *omitted* — in-silico self-consistency only | no wet lab in this replication |

Because the folding oracle differs, **absolute** ΔminPAE / ipTM values are not expected to match the paper;
the goal is to reproduce the *trends and mechanism*.

## Reproducing

```bash
conda env create -f environment.yml
conda activate rfd3-dna-replication
# CPU-side analyses run from analysis/ ; GPU generation runs via pecli (see specs/).
```

## Attribution & licensing

Method credit: Sehgal et al. 2026 and Butcher et al. 2025 (Baker lab, Institute for Protein Design, UW).
Both preprints are CC-BY 4.0. This replication is independent and not affiliated with the original authors.
Code in this repo is MIT-licensed (see `LICENSE`). The RFdiffusion3 checkpoints and pecli are governed by
their own upstream licenses.
