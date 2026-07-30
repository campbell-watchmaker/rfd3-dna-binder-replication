# Smoke Test Quick Start

## TL;DR

1. **Initialize checklist:** `python scripts/smoke_test_checklist.py --init`
2. **Follow workflow:** See `.claude/plans/pecli-smoke-test-workflow.md`
3. **Track progress:** `python scripts/smoke_test_checklist.py --show`

---

## Before You Start

### Pre-Flight Checklist (manual)

- [ ] **Validate off-target panel** against Sehgal et al. 2026, Table 1
  ```bash
  python scripts/validate_offtargets.py
  ```
  Then manually verify the 9 sequences and report any mismatches.

- [ ] **Decide H-bond conditioning subsetting**
  - Which major-groove atoms will you condition on for rfd3na?
  - Paper uses subset (not all atoms); record decision:
  ```bash
  python scripts/smoke_test_checklist.py \
    --decision hbond_subsetting "N7,O6 of purine core (positions X-Y)"
  ```

- [ ] **Verify pecli auth**
  ```bash
  pecli status  # or equivalent auth check
  ```

- [ ] **Estimate GPU budget for smoke test**
  - ~20-50 designs × 3 oracles = budget estimate
  ```bash
  python scripts/smoke_test_checklist.py \
    --decision gpu_budget_smoke "~$X for 20 designs × 3 oracles"
  ```

---

## Binder Block Stages

### Stage 0: Fold Target Duplex

```bash
# Target already prepped: targets/prnp/prnp_fold_input.json
pecli prepare protenix --input targets/prnp/prnp_fold_input.json --seeds 1
pecli submit <run_id>
# → prnp_duplex.cif

# Record completion
python scripts/smoke_test_checklist.py --stage stage_0_fold_duplex --done
```

### Stage 1: Compute Conditioning

```bash
python scripts/compute_conditioning.py \
    --duplex prnp_duplex.cif \
    --out targets/prnp/conditioning.json

python scripts/smoke_test_checklist.py --stage stage_1_conditioning --done
```

### Stage 2: Generate RFd3na Specs

```bash
python scripts/make_rfd3na_specs.py \
    --conditioning targets/prnp/conditioning.json \
    --duplex-cif prnp_duplex.cif \
    --protein-len 120-150 \
    --design-name prnp_binder \
    --out-dir specs/binder_block/rfd3na_specs

# **IMPORTANT GATE:** Edit the H-bond conditioning in each spec
# to subset major-groove atoms (not all atoms)

python scripts/smoke_test_checklist.py --stage stage_2_rfd3na_specs --done
```

### Stage 3: Diffuse Binders

```bash
for spec in specs/binder_block/rfd3na_specs/prnp_binder_ori*.json; do
    pecli prepare rfd3na --design-inputs "$spec" \
        --config specs/binder_block/sampler_config.json:_smoke_test
    pecli submit <run_id>
    # → ~10 designs per ori, ~20 total
done

# Record how many designs
python scripts/smoke_test_checklist.py \
  --stage stage_3_diffuse --done
# Then manually edit results/smoke_test_checklist.json to add n_designs
```

### Stages 4-7: CPU-Side Processing

```bash
# Stage 4: Relax with OpenMM
for pdb in <rfd3na_output>/*.pdb; do
    python scripts/relax_openmm.py --complex "$pdb" \
        --out "${pdb%.pdb}_relaxed.pdb" --dna-chains A,B
done
python scripts/smoke_test_checklist.py --stage stage_4_relax --done

# Stage 5: Sequence design (GPU, per backbone)
# ... pecli submit ligandmpnn jobs
python scripts/smoke_test_checklist.py --stage stage_5_ligandmpnn --done

# Stage 6: Refold (GPU, 3 oracles or 1 for speed)
# ... pecli submit protenix/openfold3/esmfold2
python scripts/smoke_test_checklist.py --stage stage_6_refold --done

# Stage 7: Filter and validate
python scripts/filter_binder_block.py \
    --rfd3na-designs <rfd3na_output> \
    --refolds specs/binder_block/fold_inputs \
    --pae-dir <oracle_outputs> \
    --out results/binder_block/binder_designs_filtered.csv

python scripts/smoke_test_checklist.py --stage stage_7_filter --done
# Check pass_rate and update checklist
```

---

## Specificity Block Stages

### Stages 1-5: CPU Input Prep + GPU Fold + CPU Ranking

```bash
# Stage 1: Resample sequences (GPU)
# ... pecli submit ligandmpnn resample jobs on binder-block passers

# Stage 2: On-target fold (GPU)
# ... pecli submit protenix fold

# Stage 3: Build all-by-all inputs (CPU, here)
python scripts/build_allbyall_inputs.py \
    --designs <binder_block_passers> \
    --offtargets specs/specificity_block/offtargets.json \
    --on-target-cif prnp_duplex.cif \
    --out specs/specificity_block/allbyall_inputs

# Stage 4: All-by-all fold (GPU, protenix only)
# ... pecli submit protenix jobs for all (design × DNA target) pairs

# Stage 5: Compute ΔminPAE (CPU, here)
python scripts/compute_delta_minpae.py \
    --manifest specs/specificity_block/folds_manifest.json \
    --out results/specificity_block/delta_minpae.csv

python scripts/smoke_test_checklist.py --stage stage_5_delta_minpae --done
```

---

## Post Smoke Test

After Stage 5 completes, verify:

```bash
# Check pass rates and metrics
cat results/binder_block/binder_designs_filtered.csv
cat results/specificity_block/delta_minpae.csv

# Verify sensible trends
python scripts/smoke_test_checklist.py \
  --decision metrics_match_paper_trends "yes - ΔminPAE ranking makes sense"

# Mark post-smoke-test gates
python scripts/smoke_test_checklist.py \
  --stage no_crashes \
  --stage file_handoffs_work \
  --stage cost_estimated \
  --stage decisions_reviewed \
  --done
```

---

## Troubleshooting

- **Stage hangs / no output:** Check pecli submission status, AWS credentials
- **CPU stage crashes:** Run with `--verbose` flag, check input files exist
- **Filter pass rate too low/high:** Adjust thresholds in `scripts/filter_binder_block.py`
- **ΔminPAE all negative:** Check that off-targets are actually worse than on-target

---

## Next: Production Run

Once smoke test is green:

1. Switch `sampler_config.json:_smoke_test` → `_full_run` (~1000 designs)
2. Re-run binder + specificity blocks at scale
3. Launch post-generation analyses:
   - `analysis/delta_minpae/` — validate metric on paper's published data
   - `analysis/tf_embedding/` — map designs onto natural TF sequence space
4. Generate publication figures (`figures/`)
