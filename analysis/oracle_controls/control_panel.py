#!/usr/bin/env python3
"""Control panel for validating the ΔminPAE specificity metric against a fold oracle.

WHY THIS EXISTS
---------------
The specificity block ranks designs by

    ΔminPAE = min over off-targets of minPAE(off) − minPAE(on)

computed from a fold oracle's PAE matrix. This replication substitutes open
oracles (rf3 / protenix / openfold3) for the paper's AlphaFold3. Before that
metric is trusted on our own designs we need to know it has any discriminative
power *at all* on this oracle. The paper's own binder sequences are not
published, so instead we calibrate against nature: proteins whose DNA-binding
behaviour is already known.

THE DESIGN
----------
Three protein classes, folded against a shared panel of DNA duplexes:

  specific            a natural TF with a known consensus site
  nonspecific_binder  binds duplex DNA avidly but sequence-agnostically
  nonbinder           does not bind DNA at all

Two independent readouts fall out of one all-by-all fold:

  ΔminPAE / spread    separates *specific* from *non-specific + non-binder*
  absolute minPAE     separates *non-binder* from *specific + non-specific*

and one sharper test:

  argmin correctness  for each specific TF, does min-over-panel minPAE land on
                      that TF's own cognate site? This is a hit rate, not a
                      threshold, so it needs no cross-oracle calibration.

Expected pattern if the oracle's PAE carries real specificity information:

  class               ΔminPAE   absolute minPAE   argmin on own site
  specific            large +   low on-target     yes
  nonspecific_binder  ~0        low everywhere    n/a (no cognate site)
  nonbinder           ~0        high everywhere   n/a

A NOTE ON WHY THE DNA IS PADDED TO A FIXED LENGTH
-------------------------------------------------
minPAE is a *minimum* over protein×DNA token pairs, so it is sensitive to how
many pairs exist. Cognate sites differ in length (6 bp E-box vs 17 bp λ
operator), so comparing raw sites would confound specificity with token count.
Every duplex is therefore padded to FIXED_BP with the motif centred in a fixed
neutral flank, making every PAE matrix the same shape.

A NOTE ON MSA
-------------
These folds are run WITHOUT an informative MSA (a single-sequence a3m). That is
deliberate: the de novo designs this metric will ultimately be applied to have
no meaningful MSA either. Giving natural TFs deep alignments while designs get
none would make the controls easier than the real task and inflate the apparent
discriminative power of the metric.
"""
from __future__ import annotations

FIXED_BP = 24  # every duplex is padded to this length (fits the 17-bp λ operator with flanks)

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]


# ---------------------------------------------------------------------------
# Neutral flanking sequence
# ---------------------------------------------------------------------------
# Padding must not itself contain any panel motif, or a "off-target" duplex
# would silently carry an on-target site in its flank. This flank is verified
# against every motif (both strands) by verify_panel() below -- it is not
# assumed to be clean.
NEUTRAL_FLANK = "CTGACTTGCAGTCTGACTTGCAGTCTGACTTGCAGTCTGACTTGCAGT"


# ---------------------------------------------------------------------------
# DNA panel
# ---------------------------------------------------------------------------
# Each entry: id -> (motif, provenance). Motif is the bare cognate site; it gets
# centred in NEUTRAL_FLANK to FIXED_BP by build_duplex().
#
# The five TF sites double as each other's off-targets: with 5 specific TFs in
# the panel, each TF has 1 on-target and 4 biologically-real off-targets, which
# is a far stronger negative set than scrambled sequence alone.
DNA_PANEL = {
    "zif268_site": ("GCGTGGGCGT",
                    "Zif268/EGR1 site; JASPAR MA0162.5 consensus, matches 1AAY crystal"),
    "lambda_OL1":  ("TATCACCGCCAGTGGTA",
                    "bacteriophage λ OL1 operator, 17 bp, as crystallised in 1LMB"),
    "ebox":        ("CACGTG",
                    "canonical E-box, palindromic; MAX dimeric site, matches 1AN2"),
    "hd_taatta":   ("TAATTA",
                    "homeodomain TAAT core; JASPAR MA0220.1, matches 3HDD"),
    "tata":        ("TATAAAA",
                    "TATA box; JASPAR consensus, contained in the 1CDW crystal site"),
    "prnp":        ("TGAGGAGAGGAG",
                    "THIS PROJECT'S REAL TARGET (Sehgal et al. T1). Included so the "
                    "oracle's behaviour on the actual design target is measured on the "
                    "same footing as the natural controls."),
    "scramble":    ("TCGATGCTAGCA",
                    "length-matched neutral control, no known TF motif"),
    "polygc":      ("GCGCGCGCGCGC",
                    "GC-only extreme; minimal base-edge information in the major groove"),
}

# which panel entry is each specific TF's cognate site
ON_TARGET = {
    "Zif268":    "zif268_site",
    "LambdaRep": "lambda_OL1",
    "MAX_bHLH":  "ebox",
    "Engrailed": "hd_taatta",
    "TBP":       "tata",
}


# ---------------------------------------------------------------------------
# Per-protein modelling recipe
# ---------------------------------------------------------------------------
# These are NOT cosmetic. Each entry encodes a confound found during curation
# that would otherwise silently corrupt the fold (see curated_controls.json for
# the full reasoning per entry).
RECIPE = {
    # 3 Zn2+, one per C2H2 finger. The betabetaalpha fold does not exist without
    # them, so omitting Zn would degrade the fold for reasons unrelated to
    # sequence specificity.
    "Zif268":    {"copies": 1, "ligands": ["ZN", "ZN", "ZN"]},

    # bHLH-Zip: the C-terminal leucine zipper is an obligate parallel dimer and
    # specific binding needs both basic regions. A monomer fold is unphysical.
    "MAX_bHLH":  {"copies": 2, "ligands": []},

    # CI binds the 17-bp pseudo-palindrome as a 2-fold dimer, one HTH per
    # half-site. One chain would read only half the operator.
    "LambdaRep": {"copies": 2, "ligands": []},

    # 1:1 monomer on one TAAT core is the biologically relevant complex; the
    # second copy in 3HDD is a lattice artefact.
    "Engrailed": {"copies": 1, "ligands": []},

    # Monomeric saddle. NB minor-groove reader that kinks DNA ~80 deg -- kept as
    # a deliberate stress test, not an easy positive. See notes.
    "TBP":       {"copies": 1, "ligands": []},

    "Sac7d":     {"copies": 1, "ligands": []},
    "Ubiquitin": {"copies": 1, "ligands": []},
    "GFP":       {"copies": 1, "ligands": []},
}

# Sequence repairs required before folding. Applied by load_controls() and
# recorded in the emitted manifest so the run is auditable.
SEQUENCE_FIXES = {
    # 1EMA position 65 is 'X' in the PDB one-letter code: the CRO chromophore, a
    # single modified residue formed by cyclisation of Thr65-Tyr66-Gly67. An
    # AF3-class predictor will reject 'X' or model it as UNK, distorting the
    # barrel. Expanding X -> TYG restores a standard-alphabet 238-mer.
    "GFP": ("X", "TYG"),
}


def build_duplex(motif: str, fixed_bp: int = FIXED_BP, flank: str = NEUTRAL_FLANK):
    """Centre `motif` in neutral flank, padded to exactly fixed_bp.

    Returns (sense, antisense, left_pad, right_pad).
    """
    motif = motif.upper()
    if len(motif) > fixed_bp:
        raise ValueError(f"motif {motif} ({len(motif)} bp) longer than fixed_bp={fixed_bp}")
    total_pad = fixed_bp - len(motif)
    left = total_pad // 2
    right = total_pad - left
    if left > len(flank) or right > len(flank):
        raise ValueError("neutral flank too short for requested padding")
    # take the left pad from the END of the flank so the motif junction differs
    # from the right pad's junction (avoids an accidental palindrome at the seams)
    lpad = flank[len(flank) - left:] if left else ""
    rpad = flank[:right] if right else ""
    sense = lpad + motif + rpad
    assert len(sense) == fixed_bp, (sense, len(sense))
    return sense, revcomp(sense), lpad, rpad


def build_panel(fixed_bp: int = FIXED_BP):
    """Build every DNA target in the panel at a common length."""
    out = {}
    for dna_id, (motif, prov) in DNA_PANEL.items():
        sense, anti, lpad, rpad = build_duplex(motif, fixed_bp)
        out[dna_id] = {
            "id": dna_id,
            "motif": motif,
            "motif_len": len(motif),
            "sense": sense,
            "antisense": anti,
            "left_pad": lpad,
            "right_pad": rpad,
            "provenance": prov,
        }
    return out


def verify_panel(panel: dict) -> list:
    """Check that no duplex carries a motif it is not supposed to.

    Every padded duplex is scanned (both strands) for every panel motif. A hit
    that is not the target's own motif is a contamination: it would make an
    intended off-target a partial on-target and silently compress ΔminPAE.
    Returns a list of problem strings (empty == clean).
    """
    problems = []
    motifs = {k: v[0].upper() for k, v in DNA_PANEL.items()}
    for dna_id, rec in panel.items():
        sense = rec["sense"]
        both = f"{sense}|{revcomp(sense)}"
        for mid, motif in motifs.items():
            present = motif in both or revcomp(motif) in both
            if mid == dna_id:
                if not present:
                    problems.append(
                        f"{dna_id}: own motif {motif} NOT found in built duplex {sense}")
            elif present:
                # A short motif can legitimately be a substring of a longer one
                # (e.g. TAATTA inside a longer AT-rich site). Report it either
                # way and let the caller judge -- silence here would be worse.
                problems.append(
                    f"{dna_id} ({sense}) unexpectedly contains {mid} motif {motif}")
    # the flank itself must be clean
    for mid, motif in motifs.items():
        if motif in NEUTRAL_FLANK or revcomp(motif) in NEUTRAL_FLANK:
            problems.append(f"NEUTRAL_FLANK contains {mid} motif {motif}")
    return problems


def load_controls(path: str) -> list:
    """Read curated_controls.json and apply the recorded sequence repairs."""
    import json
    blob = json.load(open(path))
    controls = blob["controls"]
    for c in controls:
        fix = SEQUENCE_FIXES.get(c["label"])
        if fix:
            bad, good = fix
            if bad in c["protein_sequence"]:
                c["_sequence_fix"] = f"replaced {bad!r} -> {good!r}"
                c["_original_length"] = len(c["protein_sequence"])
                c["protein_sequence"] = c["protein_sequence"].replace(bad, good)
                c["protein_length"] = len(c["protein_sequence"])
        r = RECIPE.get(c["label"], {"copies": 1, "ligands": []})
        c["copies"] = r["copies"]
        c["ligands"] = r["ligands"]
        c["on_target"] = ON_TARGET.get(c["label"])
        # guard: no residual non-standard amino acids
        bad_chars = set(c["protein_sequence"]) - set("ACDEFGHIKLMNPQRSTVWY")
        if bad_chars:
            raise ValueError(
                f"{c['label']} still contains non-standard residues {sorted(bad_chars)} "
                "-- add a SEQUENCE_FIXES entry before folding")
    return controls


if __name__ == "__main__":
    import json
    import os
    panel = build_panel()
    problems = verify_panel(panel)
    print(f"DNA panel: {len(panel)} targets at {FIXED_BP} bp\n")
    for r in panel.values():
        print(f"  {r['id']:12s} {r['sense']}  (motif {r['motif']}, {r['motif_len']} bp)")
    print()
    if problems:
        print(f"VERIFICATION PROBLEMS ({len(problems)}):")
        for p in problems:
            print("  ! " + p)
    else:
        print("verification: clean -- every duplex carries exactly its own motif")

    here = os.path.dirname(os.path.abspath(__file__))
    ctrl_path = os.path.join(here, "curated_controls.json")
    if os.path.exists(ctrl_path):
        controls = load_controls(ctrl_path)
        print(f"\ncontrols: {len(controls)}")
        for c in controls:
            fix = f"  [{c['_sequence_fix']}]" if c.get("_sequence_fix") else ""
            lig = f" +{len(c['ligands'])}Zn" if c["ligands"] else ""
            print(f"  {c['label']:11s} {c['klass']:18s} L={c['protein_length']:4d} "
                  f"x{c['copies']}{lig:5s} on_target={c['on_target'] or '-'}{fix}")
