# Binder block — pipeline runbook (PRNP-site)

The binder-block pipeline generates candidate DNA-binding proteins against the
folded PRNP-site duplex and filters them to self-consistent binders. GPU stages
run on the user's AWS account via **pecli** (`prepare` → review → `submit`);
CPU stages (conditioning geometry, OpenMM relax, filtering) run in Claude Science.

**Scale for the first pass: SMOKE TEST** (~10 designs) to validate the spec
end-to-end before committing budget. Use `sampler_config.json → _smoke_test`.

## Stage 0 — fold the target duplex (GPU, pecli)

Fold the DNA-only duplex to B-form. Input already prepared:
`targets/prnp/prnp_fold_input.json` (both strands, seed 42).

```bash
pecli prepare rf3 --input targets/prnp/prnp_fold_input.json --seed 42
pecli submit <run>
# → prnp_duplex.cif   (the folded target; feeds every downstream stage)
```

> Was `pecli prepare protenix`. Changed to rf3 for consistency with every other fold
> in the pipeline, and because a **protein-free** protenix fold turns out to be
> impossible in pecli at all (issue #189) — the DNA-only duplex is exactly that case.

## Stage 1 — compute conditioning (CPU, here)

Run the geometry driver on the *folded* duplex (not the fold input) to get ori
tokens + major-groove H-bond candidate atoms:

```bash
python scripts/compute_conditioning.py \
    --duplex prnp_duplex.cif --out targets/prnp/conditioning.json
```

This validates the major-groove geometry (purine N7 on major side, ~34°/bp
twist) and warns if N7-major < 90%.

## Stage 2 — generate rfd3na specs (CPU, here)

One spec **per ori placement** (2 for a 12-bp target) — `ori_token` is a single
[x,y,z] per run, so each ori is a separate diffusion job (paper: "~5100 scaffolds
per ori"):

```bash
python scripts/make_rfd3na_specs.py \
    --conditioning targets/prnp/conditioning.json \
    --duplex-cif   prnp_duplex.cif \
    --protein-len  120-150 \
    --design-name  prnp_binder \
    --out-dir      specs/binder_block/rfd3na_specs
```

> **H-bond conditioning: sample it, do not fix it.** For a real campaign use
>
> ```bash
> --hbond-sampling random --hbond-strand either --designs-per-ori 5100 --seed 42
> ```
>
> which emits one spec per design with an independently drawn constraint set: base
> count uniform over `--hbond-bases-min/max` (default 2–5, spanning 4–10 atoms of
> PRNP's 36), positions uniform without replacement inside that ori's window, strand
> drawn per design. Every draw is recorded in `manifest.json` under `hbond_draw`, so
> any design traces back to the exact set that produced it.
>
> This is the shape the paper describes — *"we sample a variety of placements of the
> protein center of mass … and **a diverse set of hydrogen bond (Hbond) condition
> constraints**"*. **The drawing mechanism is ours**; the paper states the diversity,
> not the method. Both atoms of a purine are always kept together (G N7+O6, A N7+N6 are
> the bidentate pairs Arg and Asn form), and terminal base pairs are excluded.
>
> `--hbond-sampling fixed` (the default, purine strand + central 3 bases = 6 of 36
> atoms) is retained so existing arms stay reproducible. It is **not** the paper's
> shape. Consistent with the measurement that 6 vs 8 atoms was indistinguishable over
> 100 refolds — the count was not the operative variable, so no fixed count is the fix.
>
> HBPLUS is **not** needed for conditioning (verified against the rfd3na source
> 2026-09-02). It is used only by the hbond *metrics* and by training-time hbond
> calculation, which is gated behind `TrainingConditionRoute("calculate_hbonds")`.
> Missing HBPLUS costs a warning and the hbond metric, not the conditioning.
>
> **Provenance, corrected 2026-08-05.** This note used to say "the paper conditions
> on a selected subset, e.g. the N7/O6 of the central G/A run". **That is not in the
> paper.** Its only Methods sentence on the topic is *"Hydrogen bond conditioning was
> applied during generation on candidate major groove donor and acceptor atoms
> (Fig. S1)"* — no count, no atom names, no strand; Fig. S1 is in the unreachable
> supplement. The subsetting rule here is ours.
>
> **The paper's actual model is sampling, not a fixed subset.** From Results: *"we
> sample a variety of placements of the protein center of mass … and **a diverse set
> of hydrogen bond (Hbond) condition constraints**"*. So a single fixed subset is the
> wrong shape whatever its size, and the current rule is a stand-in. Consistent with
> the measurement: 6 vs 8 atoms was indistinguishable over 100 refolds
> (`scripts/compare_conditioning_arms.py`).
>
> Two upstream signals suggest our rule is narrower than practice, both
> implementation rather than paper: foundry's training subsamples to roughly a third
> of candidate atoms (~12 for a 36-atom target), and its shipped 14-bp
> `na_binder_design.json` example specifies 16 atoms across **both** strands, mixing
> base *and* phosphate/sugar atoms — where we use purine-strand base edges only.

## Stage 3 — diffuse binders (GPU, pecli, per ori spec)

```bash
python scripts/submit_arm_diffusion.py \
    --spec-dir specs/binder_block/rfd3na_specs \
    --config   specs/binder_block/sampler_config.json \
    --duplex-cif targets/prnp/prnp_duplex.cif \
    --group prnp-binder --designs-per-run 10 --max-spend 5.00 --dry-run
# review, then drop --dry-run
# → per-design <id>.cif (+ <id>.pdb for protein-containing designs) + <id>.json
```

> **Corrected 2026-08-06.** This used to show
> `pecli prepare rfd3na --design-inputs "$spec" --config sampler_config.json:_smoke_test`.
> **Neither flag exists** — `pecli prepare` rejects both with
> `unknown option(s) for rfd3na: config, design_inputs`. The spec goes in via `--input`,
> and every sampler knob is its own flag (`--diffusion-batch-size`, `--n-batches`,
> `--use-classifier-free-guidance`, …). The command as written could never have run.
>
> It also omitted the **target-CIF staging step**, without which a run starts and dies
> immediately: the spec references its duplex by container path (`/workspace/…cif`) and
> `pecli prepare` stages only `config.json` plus the spec. `submit_arm_diffusion.py`
> drives prepare → stage → submit as one sequence so the staging cannot be skipped (it
> has been, twice), translates the config JSON into flags, and enforces a spend cap
> per run rather than only up front.
>
> **`cfg_features` is not a real setting.** `sampler_config.json` carries
> `cfg_features: [active_donor, active_acceptor]` with a paragraph of rationale, but
> pecli's rfd3na tool exposes no such field (`pecli/tools/rfd3na.py`), so it has never
> reached the sampler on any run this project has made. The submitter now reports it as
> skipped rather than passing it and being rejected.

Note the connector chains only the **first** design rfd3na → ligandmpnn; for
sequence design across *all* backbones, run ligandmpnn per design PDB (Stage 5).

## Stage 4 — relax each complex (CPU, here)

Open-source replacement for Rosetta FastRelax; DNA restrained, protein free:

```bash
for pdb in <rfd3na output>/*.pdb; do
    python scripts/relax_openmm.py --complex "$pdb" \
        --out "${pdb%.pdb}_relaxed.pdb" --dna-chains A,B
done
```

## Stage 5 — sequence design (GPU, pecli, per backbone)

```bash
for pdb in <relaxed>/*_relaxed.pdb; do
    pecli prepare ligandmpnn --input "$pdb" \
        --config specs/binder_block/ligandmpnn_config.json
    pecli submit <run>
done
# → <backbone>.fa per backbone: the WT input record, then 5 designs
```

`chains_to_design` is now set to `C` in the config, **read off a real relaxed
backbone**: rfd3na emits the fixed duplex first, so its output is DNA A + DNA B +
protein C. Note this is the *opposite* of the refold layout below (protein A + DNA
B,C) — that mismatch is what caused the Stage 7 RMSD bug where protein Cα found zero
overlap and the DNA aligned on the wrong strand.

## Stage 6 — refold + validate on rf3 (GPU, pecli)

Build per-design complex inputs (protein sequence + both DNA strands) and fold:

```bash
python scripts/build_fold_inputs.py \
    --ligandmpnn-dir <ligandmpnn raw dir> \
    --dna TGAGGAGAGGAG \
    --out-dir folds/refold                # + folds_manifest.json
for cj in folds/refold/*.json; do
    [ "$(basename "$cj")" = folds_manifest.json ] && continue
    pecli prepare rf3 --input "$cj" --diffusion-batch-size 1 --seed 42
    pecli submit <run>
done
```

`--ligandmpnn-dir` reads the `.fa` files directly and drops LigandMPNN's WT input
record by its **missing `id=`** rather than by position (dropping record 0 blindly
deletes a real design from any already-filtered FASTA). No template here,
deliberately: the gate below asks whether the designed sequence *independently* folds
back into its backbone, and a template hands it the answer.

**Oracle: rf3, settled empirically** — 4/5 argmin on the natural-TF control panel vs
protenix 2/5 and openfold3 0/5, at $0.020/fold (~2× cheaper than protenix, ~14×
cheaper than openfold3). See `analysis/oracle_controls/RESULTS.md`. esmfold2 ties rf3
on discrimination and is better calibrated but costs 10×, so it is the spot-check
oracle for top-ranked designs, not the panel-wide one.

**Cost lever at scale:** `--top-n-per-backbone 1` folds one sequence per backbone
instead of five. Refolding is 89% of a backbone's $0.112 cost, and the smoke test
measured per-backbone spread in major-groove H-bonds far exceeding within-backbone
spread (`[0,1,0,0,0]` vs `[13,0,13,0,7]`), i.e. backbone quality dominates sequence
choice. Record it as a deviation — the paper folds 5/backbone at this gate.

## Stage 7 — filter + rank (CPU, here)

```bash
python scripts/build_filter_manifest.py \
    --refold-manifest folds/refold/folds_manifest.json \
    --relaxed-dir relaxed --raw-dir raw/stage6 \
    --out filter_manifest.json
python scripts/filter_binder_block.py \
    --manifest filter_manifest.json \
    --stage pre_resample \
    --out results/binder_block/passers_pre_resample.csv \
    --oracle-comparison results/binder_block/all_designs.csv
```

Gates (paper): DNA-aligned protein Cα-RMSD < 8 Å → LigandMPNN resample → **< 3 Å,
ipTM > 0.7**, high H-bond counts. `--stage` selects which: `pre_resample` uses 8 Å and
**no ipTM gate**, `post_resample` uses 3 Å + ipTM. Passing `--iptm-gate` at
`pre_resample` silently does nothing.

The comparison CSV also carries **`min_pae`**, read from each refold's
`*_confidences.json`, which is what the specificity block's entry gate needs; add
`--min-pae-gate 6.6` to apply it here (that cut is recalibrated for rf3 — see
`specs/specificity_block/PIPELINE.md` Stage 0, *not* the paper's 1.25).

## Hand-off convention

pecli GPU outputs (CIF/PDB/FASTA) come back to Claude Science as artifacts or
into `results/binder_block/`; CPU filtering runs here and commits the ranked
CSVs. Keep the GPU run ids in `docs/replication_log.md` for provenance.
