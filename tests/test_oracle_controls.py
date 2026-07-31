"""Tests for the ΔminPAE oracle control panel.

These are correctness gates, not style checks. Each one guards a failure mode
that would silently produce a plausible-looking but meaningless ΔminPAE:

  - a padded off-target duplex that accidentally contains an on-target motif
    (compresses ΔminPAE toward zero for the wrong reason)
  - unequal duplex lengths (confounds specificity with PAE token count)
  - a DNA strand that rf3 would silently fold as RNA
  - a non-standard residue reaching the folder
  - chain-label mismatch producing an empty protein-DNA PAE block
"""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis", "oracle_controls"))

from compute_delta_minpae import _norm_chain, min_pae, protein_dna_token_masks  # noqa: E402
import control_panel as cp  # noqa: E402
import build_control_folds as bcf  # noqa: E402

CONTROLS_JSON = os.path.join(
    os.path.dirname(__file__), "..", "analysis", "oracle_controls", "curated_controls.json")


# --------------------------------------------------------------------------
# DNA panel
# --------------------------------------------------------------------------

def test_panel_is_length_matched():
    """minPAE is a min over protein x DNA token pairs, so unequal duplex lengths
    would confound specificity with the number of pairs available."""
    panel = cp.build_panel()
    lens = {len(r["sense"]) for r in panel.values()} | {len(r["antisense"]) for r in panel.values()}
    assert lens == {cp.FIXED_BP}, f"duplexes are not all {cp.FIXED_BP} bp: {lens}"


def test_panel_has_no_motif_contamination():
    """An intended off-target must not carry another target's motif."""
    problems = cp.verify_panel(cp.build_panel())
    assert problems == [], "panel contamination:\n" + "\n".join(problems)


def test_every_duplex_is_a_true_reverse_complement():
    for r in cp.build_panel().values():
        assert cp.revcomp(r["sense"]) == r["antisense"], r["id"]


def test_each_specific_tf_has_its_on_target_in_the_panel():
    panel = cp.build_panel()
    for prot, dna_id in cp.ON_TARGET.items():
        assert dna_id in panel, f"{prot} on-target {dna_id} missing from panel"


def test_panel_gives_each_specific_tf_multiple_real_offtargets():
    """The point of a shared panel: every TF's site is an off-target for the others."""
    panel = cp.build_panel()
    for prot, on in cp.ON_TARGET.items():
        offs = [d for d in panel if d != on]
        assert len(offs) >= 4, f"{prot} has only {len(offs)} off-targets"


# --------------------------------------------------------------------------
# Controls / sequence hygiene
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def controls():
    if not os.path.exists(CONTROLS_JSON):
        pytest.skip("curated_controls.json not present")
    return cp.load_controls(CONTROLS_JSON)


def test_no_nonstandard_residues_reach_the_folder(controls):
    """1EMA's position 65 is 'X' (the CRO chromophore). An AF3-class predictor
    rejects X or models it as UNK, distorting the barrel -- load_controls must
    have expanded it to TYG."""
    for c in controls:
        bad = set(c["protein_sequence"]) - set("ACDEFGHIKLMNPQRSTVWY")
        assert not bad, f"{c['label']} has non-standard residues {sorted(bad)}"


def test_gfp_chromophore_was_expanded(controls):
    gfp = next(c for c in controls if c["label"] == "GFP")
    assert gfp.get("_sequence_fix"), "GFP X->TYG repair was not applied"
    assert gfp["protein_length"] == gfp["_original_length"] + 2


def test_obligate_dimers_are_modelled_as_dimers(controls):
    """MAX's leucine zipper and λ repressor's operator half-sites both require two
    chains; a monomer fold would look bad for reasons unrelated to specificity."""
    for label in ("MAX_bHLH", "LambdaRep"):
        c = next(c for c in controls if c["label"] == label)
        assert c["copies"] == 2, f"{label} must be folded as a dimer"


def test_zinc_finger_gets_its_zinc(controls):
    """Without Zn2+ the C2H2 ββα fold does not exist."""
    zif = next(c for c in controls if c["label"] == "Zif268")
    assert zif["ligands"].count("ZN") == 3


def test_classes_cover_all_three_arms(controls):
    got = {c["klass"] for c in controls}
    assert got == {"specific", "nonspecific_binder", "nonbinder"}


# --------------------------------------------------------------------------
# Fold-input emitters
# --------------------------------------------------------------------------

def test_rf3_marks_dna_explicitly():
    """rf3 infers chain_type from the alphabet and tests the all-RNA branch before
    all-DNA, so a T-less strand would silently fold as RNA. chain_type must be
    written for every DNA component."""
    spec, prot, dna = bcf.build_rf3("t", "MKV", 1, [], "GCGCGCGC", "GCGCGCGC")
    comps = spec[0]["components"]
    dna_comps = [c for c in comps if c.get("chain_id") in dna]
    assert len(dna_comps) == 2
    for c in dna_comps:
        assert c["chain_type"] == "polydeoxyribonucleotide", \
            "T-less DNA without explicit chain_type would fold as RNA"


def test_rf3_dimer_uses_distinct_chain_ids():
    spec, prot, dna = bcf.build_rf3("t", "MKV", 2, [], "ACGT", "ACGT")
    assert len(prot) == 2 and len(set(prot)) == 2
    ids = [c["chain_id"] for c in spec[0]["components"]]
    assert len(ids) == len(set(ids)), "duplicate chain_id would be rejected by atomworks"


def test_rf3_ligands_are_separate_components():
    spec, _, _ = bcf.build_rf3("t", "MKV", 1, ["ZN", "ZN", "ZN"], "ACGT", "ACGT")
    zn = [c for c in spec[0]["components"] if c.get("ccd_code") == "ZN"]
    assert len(zn) == 3, "rf3 has no count field; each ion is its own component"
    assert len({c["chain_id"] for c in zn}) == 3


def test_rf3_carries_a_single_sequence_msa():
    """An rf3 input with no MSA auto-routes to a paid msa -> rf3 pipeline."""
    spec, _, _ = bcf.build_rf3("t", "MKV", 1, [], "ACGT", "ACGT")
    p = spec[0]["components"][0]
    assert p["_pecli_rf3_msa_a3m"].startswith(">")
    assert p["_pecli_rf3_msa_a3m"].count(">") == 1, "must be single-sequence (MSA-free)"


def test_protenix_duplex_is_two_entities_not_a_count():
    """count: 2 would duplicate one strand rather than add its complement."""
    spec, _, _ = bcf.build_protenix("t", "MKV", 1, [], "ACGTAA", "TTACGT")
    dna = [s["dnaSequence"] for s in spec[0]["sequences"] if "dnaSequence" in s]
    assert len(dna) == 2
    assert all(d["count"] == 1 for d in dna)
    assert dna[0]["sequence"] != dna[1]["sequence"]


def test_protenix_uses_top_level_list_and_required_count():
    """The engine rejects a top-level object, and reads entity["count"] with no
    default. The previous prnp_fold_input.json failed both."""
    spec, _, _ = bcf.build_protenix("t", "MKV", 2, ["ZN"], "ACGT", "ACGT")
    assert isinstance(spec, list)
    for item in spec[0]["sequences"]:
        assert len(item) == 1, "each sequences item must hold exactly one entity key"
        entity = next(iter(item.values()))
        assert "count" in entity


def test_protenix_carries_unpaired_msa():
    """Without a precomputed MSA, pecli auto-routes to a paid msa -> protenix
    pipeline; --use-msa false does NOT suppress that."""
    spec, _, _ = bcf.build_protenix("t", "MKV", 1, [], "ACGT", "ACGT")
    pc = spec[0]["sequences"][0]["proteinChain"]
    assert pc["unpairedMsa"].count(">") == 1


def test_protenix_entity_order_matches_declared_chains():
    """compute_control_metrics resolves protenix protein-vs-DNA positionally, so
    the emitter's declared chain order must match the entity order."""
    spec, prot, dna = bcf.build_protenix("t", "MKVMKV", 2, [], "ACGT", "ACGT")
    assert prot == ["A", "B"] and dna == ["C", "D"]


# --------------------------------------------------------------------------
# PAE plumbing
# --------------------------------------------------------------------------

def test_chain_label_normalisation_handles_rf3_entity_suffix():
    """rf3 emits token_chain_ids like "A_1". An exact match against "A" yields an
    all-False mask, which used to surface as a confusing "empty PAE block"."""
    assert _norm_chain("A_1") == "A"
    assert _norm_chain("H_12") == "H"
    assert _norm_chain("A") == "A"
    assert _norm_chain("AA_1") == "AA"
    # a non-numeric suffix is part of the name, not an entity index
    assert _norm_chain("chain_x") == "chain_x"


def test_masks_match_across_label_conventions():
    chains_rf3 = ["A_1"] * 3 + ["B_1"] * 2
    prot, dna = protein_dna_token_masks(5, chains_rf3, "A", ["B"], None, None)
    assert prot.tolist() == [True, True, True, False, False]
    assert dna.tolist() == [False, False, False, True, True]


def test_empty_mask_raises_a_diagnostic_error():
    with pytest.raises(ValueError, match="matched no"):
        protein_dna_token_masks(4, ["X_1"] * 4, "A", ["B"], None, None)


def test_min_pae_takes_the_global_interchain_minimum_both_orientations():
    pae = np.array([
        [0.0, 9.0, 9.0, 9.0],
        [9.0, 0.0, 9.0, 2.5],   # protein row 1 -> dna col 3
        [9.0, 9.0, 0.0, 9.0],
        [9.0, 1.5, 9.0, 0.0],   # dna row 3 -> protein col 1 (lower)
    ])
    prot = np.array([True, True, True, False])
    dna = np.array([False, False, False, True])
    assert min_pae(pae, prot, dna) == pytest.approx(1.5)


def test_delta_minpae_ranks_specific_above_promiscuous():
    """Sanity on the metric's direction using the summary arithmetic directly."""
    specific = {"on": 1.0, "off1": 8.0, "off2": 7.0}
    promiscuous = {"on": 1.0, "off1": 1.2, "off2": 1.1}
    d_spec = min(v for k, v in specific.items() if k != "on") - specific["on"]
    d_prom = min(v for k, v in promiscuous.items() if k != "on") - promiscuous["on"]
    assert d_spec > d_prom
    assert d_spec == pytest.approx(6.0)
    assert d_prom == pytest.approx(0.1)


# --------------------------------------------------------------------------
# esmfold2 (complex JSON emitter + two-file PAE)
# --------------------------------------------------------------------------

def test_esmfold2_writes_both_dna_strands():
    """esmfold2 does NOT auto-generate the complementary strand (verified against
    pecli's _validate_esmfold2_complex, which warns when given exactly one), so a
    duplex must be two dna chains."""
    spec, prot, dna = bcf.build_esmfold2("t", "MKV", 1, [], "ACGTAA", "TTACGT")
    dchains = [c for c in spec["chains"] if c["type"] == "dna"]
    assert len(dchains) == 2 and len(dna) == 2
    assert dchains[0]["sequence"] != dchains[1]["sequence"]
    assert isinstance(spec, dict) and spec["id"] == "t"   # object, not a list


def test_esmfold2_puts_all_zinc_in_one_ligand_chain():
    """Two or more ligand chains make the container refuse to emit a token->chain
    map (token_chain_ids = null), which is what slices the PAE. One ligand chain
    with a repeated ccd keeps the map exact."""
    spec, _, _ = bcf.build_esmfold2("t", "MKV", 1, ["ZN", "ZN", "ZN"], "ACGT", "ACGT")
    lig = [c for c in spec["chains"] if c["type"] == "ligand"]
    assert len(lig) == 1, "multiple ligand chains would null out token_chain_ids"
    assert lig[0]["ccd"] == ["ZN", "ZN", "ZN"]


def test_esmfold2_homodimer_is_two_chains_with_distinct_ids():
    spec, prot, dna = bcf.build_esmfold2("t", "MKV", 2, [], "ACGT", "ACGT")
    pchains = [c for c in spec["chains"] if c["type"] == "protein"]
    assert len(pchains) == 2 and len({c["id"] for c in pchains}) == 2
    assert prot == ["A", "B"] and dna == ["C", "D"]
    assert all("count" not in c for c in spec["chains"]), "no count field in the schema"


def test_esmfold2_carries_no_msa_field():
    """esmfold2 declares no requires=["msa"], so there is no carrier to add; an
    unknown extra key would just be ignored (or rejected) by the container."""
    spec, _, _ = bcf.build_esmfold2("t", "MKV", 1, [], "ACGT", "ACGT")
    blob = json.dumps(spec)
    assert "msa" not in blob.lower() and "a3m" not in blob.lower()


def _esmfold2_rec(tmp_path, ids, mapping, n):
    pae = (np.arange(n * n, dtype="float16").reshape(n, n) % 30).astype("float16")
    np.save(tmp_path / "x_pae.npy", pae)
    (tmp_path / "x_pae_tokens.json").write_text(
        json.dumps({"token_chain_ids": ids, "mapping": mapping}))
    return {"fold_id": "t", "oracle": "esmfold2",
            "pae_path": str(tmp_path / "x_pae.npy"),
            "pae_tokens_path": str(tmp_path / "x_pae_tokens.json"),
            "protein_chains": ["A"], "dna_chains": ["B", "C"],
            "protein_len": 5, "protein_copies": 1, "dna_len": 4,
            "ligands": ["ZN", "ZN", "ZN"]}


def test_esmfold2_pae_is_cast_from_float16_and_sliced_by_token_map(tmp_path):
    import compute_control_metrics as ccm
    ids = ["A"] * 5 + ["B"] * 4 + ["C"] * 4 + ["D"] * 3
    rec = _esmfold2_rec(tmp_path, ids, "positional: one token per polymer residue", 16)
    pae, prot, dna, mapping = ccm._load_esmfold2(rec)
    assert pae.dtype == np.float32, "float16 arithmetic would quietly lose precision"
    assert int(prot.sum()) == 5 and int(dna.sum()) == 8
    assert mapping == "exact"


def test_esmfold2_falls_back_positionally_only_when_the_total_reconciles(tmp_path):
    """ZN is monoatomic, so 5 + 4 + 4 + 3x1 = 16 tokens is fully determinate and
    the fallback is sound; any other total means the composition is not what we
    think it is and the fold must yield no number at all."""
    import compute_control_metrics as ccm
    rec = _esmfold2_rec(tmp_path, None, "unavailable: 3 ligand chain(s)", 16)
    pae, prot, dna, mapping = ccm._load_esmfold2(rec)
    assert mapping == "positional_fallback"
    assert int(prot.sum()) == 5 and int(dna.sum()) == 8

    bad = _esmfold2_rec(tmp_path, None, "unavailable: 3 ligand chain(s)", 20)
    with pytest.raises(ValueError, match="Refusing to guess"):
        ccm._load_esmfold2(bad)


def test_esmfold2_refuses_a_matrix_with_no_token_sidecar(tmp_path):
    import compute_control_metrics as ccm
    rec = _esmfold2_rec(tmp_path, ["A"] * 16, "positional: ...", 16)
    rec["pae_tokens_path"] = None
    with pytest.raises(ValueError, match="sidecar"):
        ccm._load_esmfold2(rec)
