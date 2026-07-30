# Oracle controls — does ΔminPAE work on an open oracle at all?

The specificity block's whole ranking rests on

    ΔminPAE = min over off-targets of minPAE(off) − minPAE(on)

computed from a fold oracle's PAE matrix. The paper used AlphaFold3; this
replication substitutes open AF3-class models. Before that metric is trusted on
our own designs, it needs to be shown to carry specificity information *on the
oracle we actually have*.

The paper's own binder sequences are unpublished, so we calibrate against
nature instead: proteins whose DNA-binding behaviour is already known.

## Design

Eight proteins in three classes, folded against a shared eight-member DNA panel
(all-by-all, 64 complexes per oracle).

| protein | class | PDB | cognate site | modelling notes |
|---|---|---|---|---|
| Zif268 | specific | 1AAY | `GCGTGGGCGT` | **+3 Zn²⁺** — the ββα fold does not exist without them |
| λ repressor | specific | 1LMB | `TATCACCGCCAGTGGTA` | 2 chains — binds the 17-bp pseudo-palindrome as a dimer |
| MAX bHLH-Zip | specific | 1AN2 | `CACGTG` | 2 chains — leucine zipper is an obligate dimer |
| Engrailed HD | specific | 3HDD | `TAATTA` | 1 chain (the 2nd copy in 3HDD is a lattice artefact) |
| TBP | specific | 1CDW | `TATAAAA` | 1 chain; **minor-groove reader, kinks DNA ~80°** |
| Sac7d | non-specific binder | 1AZP | — | binds any duplex; minor groove, kinks ~61° |
| Ubiquitin | non-binder | 1UBQ | — | pI ≈ 6.8, so no basic-patch artefact |
| GFP | non-binder | 1EMA | — | **`X`→`TYG`**: PDB pos. 65 is the CRO chromophore |

DNA panel: the five cognate sites above, plus `TGAGGAGAGGAG` (**this project's
real PRNP target**, so the oracle's behaviour on it is measured on the same
footing as the natural controls), a neutral scramble, and a poly-GC extreme.

Because the five TF sites are all in one panel, each specific TF gets 1
on-target and 4 *biologically real* off-targets — a much stronger negative set
than scrambled sequence alone.

Every duplex is padded to **24 bp** with the motif centred in a verified-neutral
flank. This is not cosmetic: minPAE is a minimum over protein×DNA token pairs,
so unequal duplex lengths would confound specificity with token count.
`control_panel.verify_panel()` re-scans every built duplex (both strands) for
every panel motif and refuses to emit folds if an intended off-target
accidentally contains an on-target site.

### Expected pattern

| class | ΔminPAE | absolute minPAE | argmin on own site |
|---|---|---|---|
| specific | large positive | low on-target | yes |
| non-specific binder | ~0 | low everywhere | n/a |
| non-binder | ~0 | high everywhere | n/a |

Two independent axes fall out of one all-by-all: **ΔminPAE/spread** separates
specific from the rest, **absolute minPAE** separates non-binders from actual
binders. The sharpest single test is **argmin correctness** — for each specific
TF, does the panel-wide minimum minPAE land on its own cognate site? Being a
rank statistic it needs no cross-oracle calibration, so it compares oracles
fairly even though their PAE scales differ.

### MSA: deliberately absent

Both oracles are run **MSA-free** (a single-sequence a3m). The de novo designs
this metric will ultimately rank have no meaningful alignment either. Giving
natural TFs deep MSAs while designs get none would make the controls strictly
easier than the real task and inflate the metric's apparent power.

## Findings that changed the pipeline

These came out of inspecting real pecli runs and reading pecli/foundry/protenix
source, and several contradict what `docs/replication_log.md` assumed.

1. **protenix does not emit a per-token PAE matrix by default.** A completed
   protenix run retains only `*_summary_confidence_sample_0.json` — scalars plus
   2×2 `chain_pair_*` aggregates. No PAE. The switch is
   `--need-atom-confidence true`, which additionally writes
   `*_full_data_sample_<rank>.json` containing `token_pair_pae`. The array is
   always computed in memory; it is simply never written to disk otherwise. So
   the replication log's plan ("protenix is primary") was viable only with a
   flag nobody had set.

2. **rf3 emits a full PAE natively, and is the cheapest oracle.**
   `*_confidences.json` carries `pae [N,N]`, `token_chain_ids`, `token_res_ids`
   — directly parseable by the existing `scripts/compute_delta_minpae.py`. At
   ~$0.06/fold vs protenix's ~$0.12 (pecli's own realised figures, n≈90 each),
   **rf3 is the better primary specificity oracle.**

3. **protenix's PAE has no chain labels** — only integer `token_asym_id`. Protein
   vs DNA must be resolved positionally from input entity order, so
   `compute_control_metrics.py` asserts the resulting token counts against the
   submitted sequence lengths rather than trusting the mapping.

4. **`targets/prnp/prnp_fold_input.json` was invalid and had never been folded.**
   It used the AF3-server spelling `{"dna": {...}}` inside a top-level object;
   protenix requires a top-level **list**, the key **`dnaSequence`**, and
   **`count`** on every entity. pecli's own prepare-time validator rejects the
   old form, so smoke-test Stage 0 would have failed immediately. Fixed.

5. **`scripts/build_allbyall_inputs.py` emits the wrong schema** — `{"id",
   "chains":[{"type": "dna"}]}`, which no pecli oracle accepts. It needs porting
   to the per-oracle emitters used here before the specificity block can run.

6. **rf3 silently folds T-less DNA as RNA.** `chain_type` is inferred from the
   alphabet and the all-RNA branch is tested before all-DNA, so a strand like
   `GCGCGCGC` becomes RNA with no error. Always write
   `chain_type: "polydeoxyribonucleotide"` explicitly.

7. **`--use-msa false` does not stop pecli auto-routing to a paid MSA step.**
   Only a precomputed MSA in the input does (`unpairedMsa` for protenix, the
   `_pecli_rf3_msa_a3m` carrier for rf3). Verified: `--use-msa false` still
   prepared an `msa → protenix` pipeline. At ~$2.46/MSA this matters.

8. **Sehgal et al. Table 1 decoys are verified.** All 9 transcribed sequences
   plus PRNP match the paper's full text, decoy-controlled (4 fabricated strings
   were correctly reported absent). One real ambiguity: `GGGCTTGCGA` appears
   labelled both *Oct4-gRNA2* and *Dux4-gRNA2*. See
   `curated_controls.json → table1_verification`.

## Files

| file | what |
|---|---|
| `curated_controls.json` | the 8 controls with verbatim RCSB sequences, motif provenance, and the full per-entry confound analysis; plus the Table 1 verification |
| `control_panel.py` | DNA panel construction, neutral-flank verification, per-protein modelling recipe (copies / ligands / sequence repairs) |
| `build_control_folds.py` | emits 64 fold inputs + manifest per oracle (rf3 and protenix schemas) |
| `compute_control_metrics.py` | minPAE, ΔminPAE, spread, argmin correctness; per-fold CSV + per-protein summary |

## Reproduce

```bash
python analysis/oracle_controls/control_panel.py            # build + verify the panel
python analysis/oracle_controls/build_control_folds.py --oracle rf3      --out-dir folds/rf3
python analysis/oracle_controls/build_control_folds.py --oracle protenix --out-dir folds/protenix

# GPU, via pecli. rf3 needs no extra flags; protenix needs the PAE flag.
pecli prepare rf3      --input folds/rf3/<id>.json      --diffusion-batch-size 1 --seed 42
pecli prepare protenix --input folds/protenix/<id>.json --need-atom-confidence true --sample 1 --seeds 42
pecli submit <run> -y

# fill pae_path in the manifest from the downloaded results, then
python analysis/oracle_controls/compute_control_metrics.py \
    --manifest folds/rf3/folds_manifest.json \
    --out-csv results/oracle_controls/rf3_folds.csv \
    --out-summary results/oracle_controls/rf3_summary.csv
```

## Cost

64 folds/oracle at pecli's realised per-run cost: **rf3 ≈ $4**, **protenix ≈ $8**.
Both run on the default A10G (108–286 tokens per complex, well under the ~700
token limit).
