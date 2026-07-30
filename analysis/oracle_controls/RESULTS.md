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
