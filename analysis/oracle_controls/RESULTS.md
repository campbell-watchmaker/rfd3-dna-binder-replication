# Results — ΔminPAE on natural DNA-binding controls

128 folds (8 proteins × 8 DNA targets × 2 oracles), all MSA-free, all on the
default A10G. **128/128 runs SUCCEEDED, 0 failures, 0 scoring errors.**

Realised GPU spend **$4.06** over 3.0 GPU-h — rf3 $1.38 for 64 folds
(**$0.022/fold**), protenix $2.68 (**$0.042/fold**). That is well under the
~$0.06/$0.13 historical per-run averages these were budgeted against, because
these complexes are small (108–286 tokens) and carry no MSA. Worth recording for
sizing the real specificity block: at rf3's measured rate, a 73-design × 46-target
all-by-all (3,358 folds) is roughly **$75**.

Figure: `figures/oracle_controls/oracle_controls.png`
Data: `results/oracle_controls/{rf3,protenix}_{folds,summary}.csv`

## rf3

| protein | class | on-target minPAE | ΔminPAE | spread | argmin | on own site |
|---|---|---|---|---|---|---|
| Zif268 | specific | **1.07** | +1.60 | 8.55 | zif268_site | yes |
| MAX_bHLH | specific | **2.15** | +1.83 | 4.12 | ebox | yes |
| LambdaRep | specific | **3.87** | +1.93 | 6.01 | lambda_OL1 | yes |
| Engrailed | specific | **3.88** | +0.26 | 10.18 | hd_taatta | yes |
| TBP | specific | 15.73 | −2.16 | 2.27 | scramble | no |
| Sac7d | non-specific binder | — | — | 4.54 | scramble (6.60) | n/a |
| GFP | non-binder | — | — | 1.89 | (15.72) | n/a |
| Ubiquitin | non-binder | — | — | **0.94** | (16.23) | n/a |

**argmin hit rate 4/5**, against 0.62 hits expected if the panel-wide minimum
fell at random (8 targets, 5 TFs).

## protenix

| protein | class | on-target minPAE | ΔminPAE | spread | argmin | on own site |
|---|---|---|---|---|---|---|
| Zif268 | specific | 0.52 | +0.72 | 6.63 | zif268_site | yes |
| Engrailed | specific | 0.82 | +0.14 | 3.63 | hd_taatta | yes |
| MAX_bHLH | specific | 0.77 | −0.06 | 3.10 | polygc | no |
| LambdaRep | specific | 7.90 | −3.40 | 5.37 | ebox | no |
| TBP | specific | 8.84 | −0.09 | 4.07 | prnp | no |
| Sac7d | non-specific binder | — | — | 5.01 | ebox (**0.61**) | n/a |
| GFP | non-binder | — | — | 4.34 | (9.75) | n/a |
| Ubiquitin | non-binder | — | — | 2.07 | (**5.78**) | n/a |

**argmin hit rate 2/5.**

## Class separation on absolute minPAE

Best minPAE achieved anywhere on the panel, binders vs non-binders:

| oracle | binders | non-binders | gap | gap excluding TBP |
|---|---|---|---|---|
| rf3 | 1.07 – 13.57 | 15.72 – 16.23 | **+2.15 Å** | **+9.12 Å** (1.07–6.60) |
| protenix | 0.52 – 8.75 | 5.78 – 9.75 | **−2.97 Å** (overlap) | +1.28 Å (0.52–4.50) |

Under protenix, Sac7d — the sequence-*non*-specific binder — attains the lowest
minPAE of any protein in the panel (0.61 Å), and ubiquitin (a non-binder,
5.78 Å) scores better than λ repressor does on its own operator (7.90 Å).

## What this does and does not establish

Reported plainly:

- On rf3, the three protein classes separate on the two axes the design
  predicted: specific TFs reach low on-target minPAE with large panel spread
  (4.1–10.2 Å), non-binders sit high with the smallest spread (0.94, 1.89 Å),
  and the panel-wide argmin lands on the correct cognate site for 4 of 5 TFs.
- On protenix, they do not separate: the binder/non-binder ranges overlap, and
  argmin is 2/5.
- ΔminPAE magnitudes are small in absolute terms even where the ranking is
  correct (rf3: +0.26 to +1.93 Å). The paper's absolute thresholds should not be
  carried over; the identical Zif268 complex scores minPAE 1.07 on rf3 and 0.52
  on protenix, so these scales are not interchangeable.

Two caveats that were predicted in advance and are not post-hoc excuses — both
are recorded in `curated_controls.json` from before any fold was run:

- **TBP was flagged as a stress test, not an easy positive.** It reads the TATA
  box entirely through the *minor* groove and kinks the duplex ~80°, and
  AF3-class predictors are biased toward B-form DNA. It fails on both oracles
  (rf3 on-target minPAE 15.73, i.e. in the non-binder range). This bounds the
  claim: the metric works for major-groove readers, which is the mode rfd3na's
  ori-token + H-bond conditioning targets, and says nothing about minor-groove or
  shape readout.
- **Zif268 is near-certainly in both models' training sets** (it is among the
  most-deposited protein–DNA complexes). It is a strong positive control but a
  weak test of generalisation, so the 4/5 should not be read as a
  generalisation estimate.

Sample size is 5 specific TFs; 4/5 vs 0.62 expected is suggestive, not a
precise hit-rate estimate.

## Follow-up: why are the absolute minPAE values high, and does it matter?

Only Zif268 (1.07 Å) clears the paper's own `minPAE < 1.25` gate; MAX (2.15),
λ repressor (3.87) and Engrailed (3.88) do not, despite all three ranking their
own cognate site first. Three candidate causes were tested or resolved.

**Tested — neither MSA nor best-of-N explains it** (rf3, on-target folds):

| protein | baseline | 5 diffusion samples | Δ | + deep MSA | Δ |
|---|---|---|---|---|---|
| Zif268 | 1.07 | 1.07 | +0.00 | **0.96** | −0.11 |
| MAX_bHLH | 2.15 | 2.15 | +0.00 | — | — |
| LambdaRep | 3.87 | 3.82 | −0.05 | — | — |
| Engrailed | 3.88 | 3.88 | +0.00 | — | — |

baseline = no MSA, 1 diffusion sample (the 64-fold panel config).

Both arms were verified to have actually taken effect rather than silently
no-op'd: the 5-sample runs produced 5 sample directories, and their ranking
scores span 0.8456–0.8458 — the samples are **near-degenerate**, which is why
best-of-N buys nothing. The MSA arm used a genuine 2.9 MB a3m with real UniRef
homologs at ~81% identity, so −0.11 Å is what a deep alignment is worth here,
not an artifact of an empty alignment.

**Resolved from the paper — the actual divergence is templating.** Sehgal et al.
computed their native-TF minPAE benchmark "with the protein templated and run in
**single-sequence mode**", and the specificity block's all-by-all is templated
too (the only place templates are used in their pipeline: "the most recent AF3
prediction before the all-by-all folding was used as the template for the
protein chain").

So our MSA-free choice *matches* the paper. What differs is that we did not
template the protein chain. Templating supplies the protein fold, so minPAE
reflects interface confidence rather than carrying fold uncertainty — which is
the remaining explanation for the offset, by elimination. rf3 supports it
(per-chain CIF as a `path` component plus `template_selection`, mixable with
`seq` components), so this is directly testable and cheap. **Untested as of this
writing.**

**Why the ranking result stands regardless.** The discrimination reported above
was obtained under a config whose absolute values move by ≤0.11 Å under both a
deep MSA and 5× sampling. The 4/5 argmin rate and the +9.1 Å class gap are
therefore not artifacts of an under-powered fold setup. What the offset does
mean is narrower and already stated: **the paper's absolute thresholds must not
be carried over to rf3 outputs** — `minPAE < 1.25` would reject 4 of 5 genuine
TFs here. Note also that `ΔminPAE > 0`, which the paper reports as the criterion
that "enriched for successful designs experimentally", is satisfied by 4/5 TFs on
rf3 and 2/5 on protenix.

Absolute minPAE values for the paper's own native TFs are in Fig. S3, which is
403-locked on bioRxiv, so we cannot yet say whether 1–4 Å is normal under their
templated protocol.

## Templated arm — templating improves ΔminPAE but does NOT close the absolute offset

Run for **LambdaRep and Engrailed only** (the two TFs at ~3.9 Å where the offset
hypothesis actually bites), full 8-target rows so ΔminPAE and argmin stay
computable, plus Zif268 as a plumbing check. Templates are the protein chain of
each protein's own prior untemplated *prediction* — see
`make_predicted_templates.py` for why a prediction and not the crystal chain.

| protein | on-target minPAE | | ΔminPAE | | argmin | in-motif |
|---|---|---|---|---|---|---|
| | untempl. | **templ.** | untempl. | **templ.** | | |
| Zif268 | 1.07 | **1.05** | +1.60 | — | ✓ | yes |
| LambdaRep | 3.87 | **2.94** | +1.93 | **+3.43** | ✓ | yes |
| Engrailed | 3.88 | **3.99** | +0.26 | **+0.72** | ✓ | yes |

- **Absolute minPAE: not explained.** LambdaRep improved 0.93 Å but remains at
  2.94, still well above the paper's `minPAE < 1.25` gate. Engrailed got slightly
  *worse* (+0.11). So templating is not the source of the offset either.
- **ΔminPAE: improved for both.** +1.93 → +3.43 (LambdaRep) and +0.26 → +0.72
  (Engrailed). Discrimination went up, not down.
- **No confidence-for-discrimination trade.** argmin held at 3/3 and all three
  on-target interfaces stayed inside the motif window (100%, baseline 45.8%), so
  the ΔminPAE gain is not templating simply making everything confident.

So the absolute offset has now survived **MSA (−0.11 Å), 5× sampling (≈0), and
templating (−0.93 to +0.11 Å)**. It is a property of rf3's PAE calibration on
this modality, not of our fold configuration. The practical consequence is
unchanged and now well supported: **the paper's absolute minPAE thresholds cannot
be transferred to rf3 and must be recalibrated empirically.** Its
calibration-free `ΔminPAE > 0` criterion transfers fine.

Templating is nonetheless worth adopting in `specs/specificity_block/`, on the
evidence that it raises ΔminPAE — the quantity the pipeline actually ranks on —
without costing discrimination. That is also what the paper does.

## openfold3 arm — worse than chance

40 of 64 folds (the $8.00 cap stopped it; the cut happened to fall after all five
specific TFs, so the argmin test is complete and only the non-binder rows are
missing, meaning no binder/non-binder gap can be computed).

| | argmin on own cognate site | ΔminPAE |
|---|---|---|
| openfold3 | **0/5** (0.62 expected by chance) | negative for all 5 |

openfold3 does emit a full per-token PAE (`pae` key), needs no extra flag, and
folds protein+DNA — it is mechanically usable. It is simply the weakest of the
three at this task.

Cost note worth keeping: openfold3's GPU-hours are **flat at ~0.137 h regardless
of complex size** (fixed ~2.3 GB checkpoint fetch dominates), so unlike protenix
it does not get cheaper for small complexes. $0.31/fold, i.e. ~14× rf3 and ~3×
protenix for the same panel.

## Three-oracle summary

| oracle | argmin on own site | binder/non-binder gap | $/fold |
|---|---|---|---|
| **rf3** | **4/5** | **+9.12 Å** (excl. TBP) | **$0.022** |
| protenix | 2/5 | −2.97 Å (overlap) | $0.042 |
| openfold3 | 0/5 | not measured (cap) | $0.31 |

rf3 is simultaneously the most discriminative and the cheapest, by a wide margin
on both.

## Consequence for the pipeline

`docs/replication_log.md` has been updated: **rf3 becomes the primary
specificity oracle.** It emits a full per-token PAE natively (protenix needs
`--need-atom-confidence true`), costs about half as much per fold, and — on this
panel — is the one whose PAE actually separates DNA binders from non-binders.
protenix is retained as a cross-check rather than the primary.

The PRNP column is included in both heatmaps, so the oracle's behaviour on this
project's actual design target is measured on the same footing as the natural
controls. No control protein has PRNP as its cognate site, so those cells are
off-target readings for every row — they are a baseline for what an
unoptimised protein scores against our target, not a positive result.
