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
pecli prepare protenix --input targets/prnp/prnp_fold_input.json --seeds 1
pecli submit <run>
# → prnp_duplex.cif   (the folded target; feeds every downstream stage)
```

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

> **Subset the H-bond conditioning before submit.** The generator emits every
> candidate major-groove atom. Conditioning on all of them over-constrains
> diffusion — pick the handful of major-groove acceptors/donors on the
> poly-purine core you actually want the binder to read (the paper conditions on
> a selected subset, e.g. the N7/O6 of the central G/A run). Edit the
> `select_hbond_*` dicts in each spec accordingly. HBPLUS must be installed on
> the GPU side for H-bond conditioning to work.

## Stage 3 — diffuse binders (GPU, pecli, per ori spec)

```bash
for spec in specs/binder_block/rfd3na_specs/prnp_binder_ori*.json; do
    pecli prepare rfd3na --design-inputs "$spec" \
        --config specs/binder_block/sampler_config.json:_smoke_test
    pecli submit <run>
done
# → per-design <id>.cif (+ <id>.pdb for protein-containing designs) + <id>.json
```

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
    # set chains_to_design to the designed protein chain (not the DNA chains)
    pecli submit <run>
done
# → FASTA of 5 sequences/backbone with overall_confidence / ligand_confidence
```

## Stage 6 — refold + validate, THREE oracles (GPU, pecli)

Build per-design complex inputs (protein sequence + both DNA strands) and fold
with each oracle for the comparison (see `fold_config.json`):

```bash
python scripts/build_fold_inputs.py \
    --fasta <ligandmpnn output>.fasta \
    --dna TGAGGAGAGGAG \
    --out-dir specs/binder_block/fold_inputs
for oracle in protenix openfold3 esmfold2; do
    for cj in specs/binder_block/fold_inputs/*.json; do
        pecli prepare $oracle --input "$cj"
        pecli submit <run>
    done
done
```

## Stage 7 — filter + rank (CPU, here)

```bash
python scripts/filter_binder_block.py \
    --designs <refolded cifs, tagged by oracle> \
    --target-dna TGAGGAGAGGAG \
    --out results/binder_block/passers.csv \
    --oracle-comparison results/binder_block/oracle_comparison.csv
```

Gates (paper): DNA-aligned protein Cα-RMSD < 8 Å → resample → **< 3 Å, ipTM >
0.7**, high H-bond counts. The oracle-comparison CSV records RMSD/ipTM/H-bonds
**and** runtime per oracle for the protenix-vs-openfold3-vs-esmfold2 writeup.

## Hand-off convention

pecli GPU outputs (CIF/PDB/FASTA) come back to Claude Science as artifacts or
into `results/binder_block/`; CPU filtering runs here and commits the ranked
CSVs. Keep the GPU run ids in `docs/replication_log.md` for provenance.
