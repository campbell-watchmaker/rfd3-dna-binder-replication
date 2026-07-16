"""Target preparation for RFdiffusion3-NA DNA-binder design.

Given a folded double-stranded B-DNA target (mmCIF/PDB), compute the two
sequence-independent conditioning inputs the Sehgal et al. 2026 pipeline feeds
to rfd3na:

  1. ori (center-of-mass) tokens: one per 6 consecutive base pairs, placed 3 Å
     toward the MAJOR groove from the centroid of that 6-bp stretch, perpendicular
     to the local helical axis.
  2. candidate MAJOR-GROOVE hydrogen-bond donor/acceptor atoms on the DNA bases,
     for H-bond conditioning.

The major-groove direction is derived from base-pair frames and validated
chemically (purine N7 must lie on the major-groove side). See
docs/replication_log.md for parameter provenance.

Usage:
    from target_prep import load_duplex, ori_tokens, hbond_candidates
    duplex = load_duplex("target.cif")
    tokens = ori_tokens(duplex, window=6, offset=3.0)
    hbonds = hbond_candidates(duplex)
"""
from __future__ import annotations
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb

# Major-groove edge atoms (Watson-Crick base numbering). Donors/acceptors that
# present into the major groove and are the canonical readout atoms for TFs.
# A: N7 (acc), N6 (don);  G: N7 (acc), O6 (acc);  C: N4 (don);  T: O4 (acc).
MAJOR_GROOVE_ATOMS = {
    "DA": {"N7": "acceptor", "N6": "donor"},
    "DG": {"N7": "acceptor", "O6": "acceptor"},
    "DC": {"N4": "donor"},
    "DT": {"O4": "acceptor"},
}
_PURINES = ("DA", "DG")


def load_duplex(path: str):
    """Read a structure, keep nucleotides, return the AtomArray (model 1)."""
    if path.lower().endswith((".cif", ".mmcif", ".bcif")):
        f = pdbx.CIFFile.read(path)
        arr = pdbx.get_structure(f, model=1)
    else:
        f = pdb.PDBFile.read(path)
        arr = f.get_structure(model=1)
    return arr[struc.filter_nucleotides(arr)]


def _res_atom(dna, chain, res_id, name):
    m = (dna.chain_id == chain) & (dna.res_id == res_id) & (dna.atom_name == name)
    return dna.coord[m][0] if m.sum() else None


def _res_name(dna, chain, res_id):
    m = (dna.chain_id == chain) & (dna.res_id == res_id)
    return dna.res_name[m][0]


def _glycosidic_N(dna, chain, res_id):
    """N9 for purines, N1 for pyrimidines (the base anchor near the helix axis)."""
    rn = _res_name(dna, chain, res_id)
    return _res_atom(dna, chain, res_id, "N9" if rn in _PURINES else "N1")


def base_pairs(dna):
    """Return list of ((chainA,resA),(chainB,resB)) Watson-Crick pairs, 5'->3' on
    the first chain, using biotite's base-pair detection."""
    pairs = struc.base_pairs(dna)  # indices into dna
    out = []
    for i, j in pairs:
        out.append((
            (dna.chain_id[i], int(dna.res_id[i])),
            (dna.chain_id[j], int(dna.res_id[j])),
        ))
    # order by first-chain residue id along the helix
    first_chain = dna.chain_id[0]
    out.sort(key=lambda p: (p[0][0] != first_chain, p[0][1]))
    return out


def bp_frames(dna, pairs):
    """Per-base-pair geometry: center (glycosidic-N midpoint), local helical axis,
    and unit major-groove direction (perpendicular to axis).

    Validated on canonical B-DNA (1BNA): mean interior twist ~34.8°/bp and all
    purine N7 atoms fall on the returned major-groove side.
    """
    centers, c1mid = [], []
    for (ca, ra), (cb, rb) in pairs:
        nA, nB = _glycosidic_N(dna, ca, ra), _glycosidic_N(dna, cb, rb)
        c1A = _res_atom(dna, ca, ra, "C1'")
        c1B = _res_atom(dna, cb, rb, "C1'")
        centers.append((nA + nB) / 2)
        c1mid.append((c1A + c1B) / 2)
    centers = np.array(centers); c1mid = np.array(c1mid)

    n = len(centers)
    axis = np.zeros_like(centers)
    for i in range(n):
        lo, hi = max(0, i - 1), min(n - 1, i + 1)
        v = centers[hi] - centers[lo]
        axis[i] = v / np.linalg.norm(v)

    major = np.zeros_like(centers)
    for i in range(n):
        # C1' midpoint sits on the MINOR-groove side of the base-pair center;
        # major groove is the opposite direction, taken perpendicular to axis.
        d = c1mid[i] - centers[i]
        d = d - np.dot(d, axis[i]) * axis[i]
        d = d / np.linalg.norm(d)
        major[i] = -d
    return centers, axis, major


def ori_tokens(dna, window: int = 6, offset: float = 3.0):
    """One ori token per `window` consecutive base pairs, placed `offset` Å toward
    the major groove from the stretch centroid, perpendicular to the helical axis.

    Returns list of dicts: {bp_start, bp_end, centroid, axis, major_dir, ori_xyz}.
    Uses a sliding step equal to `window` (non-overlapping stretches), plus a final
    trailing window if the sequence length isn't a multiple of `window`.
    """
    pairs = base_pairs(dna)
    centers, axis, major = bp_frames(dna, pairs)
    n = len(pairs)
    starts = list(range(0, max(1, n - window + 1), window))
    if starts and starts[-1] + window < n:
        starts.append(n - window)  # cover the tail
    tokens = []
    for s in starts:
        e = min(s + window, n)
        idx = np.arange(s, e)
        centroid = centers[idx].mean(axis=0)
        ax = axis[idx].mean(axis=0); ax /= np.linalg.norm(ax)
        mg = major[idx].mean(axis=0)
        mg = mg - np.dot(mg, ax) * ax
        mg /= np.linalg.norm(mg)
        ori = centroid + offset * mg
        tokens.append({
            "bp_start": int(s + 1), "bp_end": int(e),
            "centroid": centroid.tolist(),
            "axis": ax.tolist(),
            "major_dir": mg.tolist(),
            "ori_xyz": ori.tolist(),
        })
    return tokens


def hbond_candidates(dna):
    """Candidate major-groove H-bond donor/acceptor atoms per base.

    Returns list of dicts: {chain, res_id, res_name, atom, role, xyz}.
    """
    out = []
    seen = set()
    for chain, res_id, res_name in zip(dna.chain_id, dna.res_id, dna.res_name):
        key = (chain, int(res_id))
        if key in seen:
            continue
        seen.add(key)
        spec = MAJOR_GROOVE_ATOMS.get(res_name)
        if not spec:
            continue
        for atom_name, role in spec.items():
            xyz = _res_atom(dna, chain, int(res_id), atom_name)
            if xyz is not None:
                out.append({
                    "chain": str(chain), "res_id": int(res_id),
                    "res_name": str(res_name), "atom": atom_name,
                    "role": role, "xyz": xyz.tolist(),
                })
    return out


def validate_major_groove(dna, pairs=None):
    """Chemical self-check: fraction of purine N7 atoms on the major-groove side
    and mean interior helical twist. Returns (frac_N7_major, mean_twist_deg)."""
    pairs = pairs or base_pairs(dna)
    centers, axis, major = bp_frames(dna, pairs)
    ok = tot = 0
    for i, ((ca, ra), (cb, rb)) in enumerate(pairs):
        for ch, rid in ((ca, ra), (cb, rb)):
            if _res_name(dna, ch, rid) in _PURINES:
                n7 = _res_atom(dna, ch, rid, "N7")
                if n7 is None:
                    continue
                v = n7 - centers[i]; v = v - np.dot(v, axis[i]) * axis[i]
                v /= np.linalg.norm(v)
                tot += 1; ok += int(np.dot(v, major[i]) > 0)
    ang = [np.degrees(np.arccos(np.clip(np.dot(major[i], major[i + 1]), -1, 1)))
           for i in range(len(major) - 1)]
    interior = ang[2:-2] if len(ang) > 4 else ang
    return ok / max(tot, 1), float(np.mean(interior)) if interior else float("nan")
