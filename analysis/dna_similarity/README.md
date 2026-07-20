# DNA-similarity premise

The structural motivation for the whole campaign, established quantitatively.

**Claim:** sequence-specific DNA binding cannot work by reading backbone shape —
the B-DNA sugar-phosphate backbone is essentially sequence-independent. The
discriminating signal lives on the base edges exposed in the **major groove**,
which is exactly what rfd3na's ori-token + H-bond conditioning targets.

## Contents

- `dna_similarity_premise.ipynb` — the executed notebook (outputs embedded).
- `dna_similarity_premise.py` — the notebook source in `# %%` percent format
  (the `.ipynb` is generated from this; edit here, regenerate, re-execute).
- `dna_similarity_premise.png` — the figure.

## What it computes

1. **Backbone degeneracy.** All-vs-all Kabsch RMSD of the 11-atom
   sugar-phosphate backbone unit across the interior nucleotides of the
   Drew–Dickerson dodecamer (PDB 1BNA). Same-base vs different-base pairs:
   **0.328 vs 0.320 Å** — indistinguishable. Max deviation < 0.7 Å.
2. **Major-groove readout.** The ordered donor/acceptor/methyl signature each
   Watson–Crick base pair presents on its major-groove edge. All four base pairs
   are distinct; A:T is the mirror of T:A (same features, opposite 5′→3′ order),
   which is why recognition depends on atom *coordinates*, not an atom list.

## Reproduce

```bash
# from the repo root, with the environment.yml env active
jupyter nbconvert --to notebook --execute \
    analysis/dna_similarity/dna_similarity_premise.ipynb \
    --output dna_similarity_premise.ipynb
```

The notebook fetches `1BNA.cif` from RCSB on run (not committed).
