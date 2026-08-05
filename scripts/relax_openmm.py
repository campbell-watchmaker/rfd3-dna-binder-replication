#!/usr/bin/env python3
"""Open-source replacement for Rosetta FastRelax in the binder-block pipeline.

Sehgal et al. 2026 relax the rfd3na-diffused protein-DNA complex with Rosetta
FastRelax before LigandMPNN sequence sampling. Rosetta is free for academic use
but is not permissively licensed, so this replication uses OpenMM (MIT/LGPL)
instead: energy-minimize the complex with Amber ff14SB (protein) + the OL15 DNA
correction (bundled together in OpenMM's `amber14-all.xml`), TIP3P water,
**with DNA atoms held under a positional restraint** and the protein free to
relax. Restraining the DNA matches rfd3na's own treatment of the DNA as fixed
throughout diffusion, and keeps the relax step from moving the target off the
geometry the design was conditioned on.

This is the same class of step as AlphaFold2's post-prediction Amber relax
(clash / stereochemistry cleanup after generation) applied here to the rfd3na
output before sequence design.

Usage:
    python relax_openmm.py --complex rfd3na_output.pdb --out relaxed.pdb \
        --dna-chains B,C --restraint-k 500

`--dna-chains` names the chain IDs to restrain (the DNA duplex); every other
chain is treated as the diffused protein and left free. Runs on CPU for a
single structure in well under a minute at this size; no GPU dispatch needed.
"""
from __future__ import annotations
import argparse

from openmm import app, unit, CustomExternalForce
from openmm import LangevinMiddleIntegrator
from pdbfixer import PDBFixer


def build_system(fixer: PDBFixer, dna_chain_ids: set[str], restraint_k: float):
    ff = app.ForceField("amber14-all.xml", "tip3p.xml")
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addHydrogens(ff)

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,   # single structure, no periodic box needed for a minimization
        constraints=app.HBonds,
    )

    # Positional restraint on DNA heavy atoms: harmonic well centered on the
    # input (fixed) coordinates. Protein atoms are left unrestrained.
    restraint = CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
    restraint.addGlobalParameter("k", restraint_k * unit.kilojoule_per_mole / unit.nanometer**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    n_restrained = 0
    positions = modeller.positions
    for atom in modeller.topology.atoms():
        if atom.residue.chain.id in dna_chain_ids and atom.element is not None and atom.element.symbol != "H":
            p = positions[atom.index]
            restraint.addParticle(atom.index, [p.x, p.y, p.z])
            n_restrained += 1
    system.addForce(restraint)

    return system, modeller, n_restrained


def _strip_5prime_phosphate(complex_path: str, dna_chain_ids: set[str]) -> str:
    """Drop the 5'-terminal phosphate from each DNA chain, and return a temp PDB path.

    WHY. rfd3na emits a full phosphate (P, OP1, OP2) on *every* nucleotide including
    the 5' terminus. Amber's 5'-terminal DNA templates (DT5/DC5/...) begin at O5' and
    have no 5'-phosphate variant, so the residue instead half-matches the internal /
    3' template and createSystem dies with

        No template found for residue 0 (DT). The atoms and bonds in the residue match
        DT3, but the set of externally bonded atoms is missing 1 P atom and has 1 O
        atom too many. Is the chain terminated in a way that is unsupported by the
        force field?

    which is easy to misread as a residue-naming or missing-hydrogen problem. It is
    neither: the terminus simply carries one phosphate the force field cannot type.

    Removing it is chemically harmless here. A 5'-monophosphate on a synthetic duplex
    is arbitrary (rfd3na applies it uniformly rather than modelling a real 5' end), and
    the DNA is positionally restrained throughout the minimisation, so its terminal
    chemistry has negligible influence on how the protein relaxes against it.
    """
    import tempfile

    import numpy as np
    import biotite.structure.io.pdb as pdb

    arr = pdb.PDBFile.read(complex_path).get_structure(model=1)
    drop = np.zeros(arr.array_length(), dtype=bool)
    for ch in sorted(dna_chain_ids):
        sel = arr.chain_id == ch
        if not sel.any():
            continue
        first = arr.res_id[sel].min()
        drop |= sel & (arr.res_id == first) & np.isin(arr.atom_name, ("P", "OP1", "OP2"))
    n = int(drop.sum())
    kept = arr[~drop]

    out = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    f = pdb.PDBFile()
    f.set_structure(kept)
    f.write(out)
    print(f"  stripped {n} 5'-terminal phosphate atom(s) from DNA chains "
          f"{','.join(sorted(dna_chain_ids))}")
    return out


def relax(complex_path: str, out_path: str, dna_chain_ids: set[str], restraint_k: float):
    prepped = _strip_5prime_phosphate(complex_path, dna_chain_ids)
    fixer = PDBFixer(filename=prepped)
    # Don't invent residues. findMissingResidues() would try to fill gaps against the
    # SEQRES record; an rfd3na design has no meaningful SEQRES and any "missing" loop
    # it decided to build would be fabricated backbone, not part of the design.
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    # rfd3na writes heavy atoms only (verified: zero H in its output) while
    # amber14-all.xml is an all-atom force field, so hydrogens have to be built.
    fixer.addMissingHydrogens(7.0)

    system, modeller, n_restrained = build_system(fixer, dna_chain_ids, restraint_k)
    if n_restrained == 0:
        raise ValueError(
            f"No DNA heavy atoms restrained — check --dna-chains against the "
            f"chain IDs actually in {complex_path}."
        )

    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtosecond)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    e0 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(maxIterations=2000)
    e1 = simulation.context.getState(getEnergy=True).getPotentialEnergy()

    state = simulation.context.getState(getPositions=True)
    with open(out_path, "w") as f:
        app.PDBFile.writeFile(modeller.topology, state.getPositions(), f)

    print(f"restrained {n_restrained} DNA heavy atoms (k={restraint_k} kJ/mol/nm^2)")
    print(f"potential energy: {e0} -> {e1}")
    print(f"wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", required=True, help="rfd3na output complex (protein + DNA), PDB")
    ap.add_argument("--out", required=True, help="relaxed complex output PDB")
    ap.add_argument("--dna-chains", required=True, help="comma-separated chain IDs to restrain (the DNA duplex)")
    ap.add_argument("--restraint-k", type=float, default=500.0,
                     help="harmonic restraint constant, kJ/mol/nm^2 (default 500; stiff enough to hold DNA near input geometry)")
    args = ap.parse_args()
    relax(args.complex, args.out, set(args.dna_chains.split(",")), args.restraint_k)


if __name__ == "__main__":
    main()
