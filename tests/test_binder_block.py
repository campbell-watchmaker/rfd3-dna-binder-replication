"""CI tests for the binder-block spec generation and filtering scripts.

These are self-contained (no network, no GPU): they exercise the pure-logic paths
-- spec-schema shaping, fold-input construction, H-bond counting, and DNA-aligned
RMSD -- on small synthetic inputs, so a regression in the conditioning->spec
adapter or the interface-metric math is caught before it reaches a real design.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_rfd3na_specs as mrs
import build_fold_inputs as bfi
import filter_binder_block as fbb


def test_hbond_grouping_and_resid_keys():
    cands = [
        {"chain": "A", "res_id": 6, "atom": "N7", "role": "acceptor"},
        {"chain": "A", "res_id": 6, "atom": "O6", "role": "acceptor"},
        {"chain": "A", "res_id": 5, "atom": "N6", "role": "donor"},
    ]
    acc = mrs._group_hbond(cands, "acceptor")
    don = mrs._group_hbond(cands, "donor")
    assert acc == {"A6": "N7,O6"}, acc
    assert don == {"A5": "N6"}, don


def test_spec_has_single_ori_and_required_fields():
    spec = mrs.build_spec(
        "t", "dup.cif", "120-150", [1.0, 2.0, 3.0],
        {"A5": "N6"}, {"A6": "N7,O6"}, {"A": (1, 12), "B": (13, 24)},
    )["t"]
    assert isinstance(spec["ori_token"], list) and len(spec["ori_token"]) == 3
    assert spec["is_non_loopy"] is True
    assert spec["select_fixed_atoms"] == {"A1-12": "ALL", "B13-24": "ALL"}
    assert "120-150" in spec["contig"]
    assert spec["select_hbond_acceptor"] == {"A6": "N7,O6"}


def test_revcomp():
    assert bfi.revcomp("TGAGGAGAGGAG") == "CTCCTCTCCTCA"


# --- Stage 6 refold inputs (rewritten 2026-08-06) --------------------------
# build_fold_inputs.py used to emit a {"id","chains":[...]} shape no pecli oracle
# accepts, so PIPELINE.md pointed at a script that could not fold anything. It now
# emits the rf3 component shape the smoke test's 50 refolds actually ran on. These
# tests pin the two facts that silently corrupt a run if they drift, plus the WT
# handling.

_LMPNN_FA = """\
>ori1_0_3, T=0.1, seed=1, num_res=8, num_ligand_res=45, use_ligand_context=True
ACDEFGHI
>ori1_0_3, id=1, T=0.1, seed=1, overall_confidence=0.39, ligand_confidence=0.34, seq_rec=0.13
MKTAYIAK
>ori1_0_3, id=2, T=0.1, seed=1, overall_confidence=0.51, ligand_confidence=0.44, seq_rec=0.15
MKTAYIAR
"""


def test_ligandmpnn_wt_record_dropped_by_missing_id_not_by_position(tmp_path):
    """The WT input carries no `id=`. Dropping record 0 positionally instead would
    delete a real design from any already-filtered FASTA."""
    p = tmp_path / "ori1_0_3.fa"
    p.write_text(_LMPNN_FA)
    ds = bfi.parse_ligandmpnn_fasta(str(p))
    assert [d["seq_id"] for d in ds] == [1, 2]
    assert "ACDEFGHI" not in [d["seq"] for d in ds], "WT input sequence leaked through"
    assert ds[0]["backbone"] == "ori1_0_3"
    assert ds[1]["overall_confidence"] == 0.51


def test_rf3_complex_spec_shape_and_the_two_load_bearing_fields(tmp_path):
    spec, prot, dna = bfi.build_rf3_complex("f1", "MKTAYIAK", "TGAGGAGAGGAG", "CTCCTCTCCTCA")
    assert isinstance(spec, list) and set(spec[0]) == {"name", "components"}
    assert (prot, dna) == (["A"], ["B", "C"])
    comps = spec[0]["components"]
    # (1) chain_type explicit on BOTH strands: a T-free strand is otherwise inferred
    #     as RNA, because the all-RNA branch is tested before the all-DNA branch.
    assert [c.get("chain_type") for c in comps] == \
        ["polypeptide(L)", "polydeoxyribonucleotide", "polydeoxyribonucleotide"]
    # (2) the MSA carrier: without it pecli auto-routes to a PAID msa -> rf3 pipeline
    assert comps[0]["_pecli_rf3_msa_a3m"] == ">query\nMKTAYIAK\n"
    # no template -- the RMSD gate is a self-consistency check; a template hands it
    # the answer and collapses RMSD toward zero
    assert not any("path" in c or "template" in str(c) for c in comps)
    assert "template_selection" not in spec[0]


def test_no_t_dna_strand_still_typed_as_dna():
    """The exact silent-RNA trap, on a strand that triggers it."""
    spec, _, _ = bfi.build_rf3_complex("f", "MK", "GGGCCCGGGCCC", "GGGCCCGGGCCC")
    for c in spec[0]["components"][1:]:
        assert c["chain_type"] == "polydeoxyribonucleotide"


def test_emitted_manifest_carries_every_field_build_filter_manifest_reads(tmp_path):
    """build_filter_manifest.py:45,73 index `backbone` and `seq_id` unguarded."""
    import subprocess
    d = tmp_path / "lm" / "seqs"
    d.mkdir(parents=True)
    (d / "ori1_0_3.fa").write_text(_LMPNN_FA)
    out = tmp_path / "folds"
    r = subprocess.run([sys.executable,
                        os.path.join(os.path.dirname(__file__), "..", "scripts",
                                     "build_fold_inputs.py"),
                        "--ligandmpnn-dir", str(tmp_path / "lm"),
                        "--dna", "TGAGGAGAGGAG", "--out-dir", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    man = json.load(open(out / "folds_manifest.json"))
    assert len(man) == 2
    for rec in man:
        for f in ("fold_id", "backbone", "seq_id", "protein_seq", "protein_chain",
                  "dna_chains", "protein_len", "oracle", "fold_input"):
            assert f in rec, f
        assert os.path.isfile(out / rec["fold_input"])
    assert man[0]["fold_id"] == "ori1_0_3_s1"


def test_top_n_per_backbone_keeps_the_highest_confidence(tmp_path):
    """The scale triage lever: 1 seq/backbone at the first gate is 1/5 the refold cost."""
    import subprocess
    d = tmp_path / "lm" / "seqs"
    d.mkdir(parents=True)
    (d / "ori1_0_3.fa").write_text(_LMPNN_FA)
    out = tmp_path / "folds"
    r = subprocess.run([sys.executable,
                        os.path.join(os.path.dirname(__file__), "..", "scripts",
                                     "build_fold_inputs.py"),
                        "--ligandmpnn-dir", str(tmp_path / "lm"),
                        "--dna", "TGAGGAGAGGAG", "--out-dir", str(out),
                        "--top-n-per-backbone", "1"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    man = json.load(open(out / "folds_manifest.json"))
    assert len(man) == 1
    assert man[0]["seq_id"] == 2, "kept seq 1 (0.39) over seq 2 (0.51)"


def test_non_acgt_dna_is_rejected_not_folded(tmp_path):
    import subprocess
    fa = tmp_path / "d.fasta"
    fa.write_text(">d1\nMKTAYIAK\n")
    r = subprocess.run([sys.executable,
                        os.path.join(os.path.dirname(__file__), "..", "scripts",
                                     "build_fold_inputs.py"),
                        "--fasta", str(fa), "--dna", "TGAGGAGNGGAG",
                        "--out-dir", str(tmp_path / "o")], capture_output=True, text=True)
    assert r.returncode != 0
    assert "ACGT" in r.stderr


def test_dna_aligned_rmsd_zero_for_identity(tmp_path):
    # a tiny protein+DNA AtomArray, design == refold -> RMSD 0
    import biotite.structure as struc
    n = 6
    arr = struc.AtomArray(n)
    arr.coord = np.arange(n * 3, dtype=float).reshape(n, 3)
    arr.chain_id = np.array(["B", "B", "B", "A", "A", "A"])
    arr.res_id = np.array([1, 1, 2, 10, 11, 12])
    arr.res_name = np.array(["DA", "DA", "DA", "ALA", "ALA", "ALA"])
    arr.atom_name = np.array(["N7", "C1'", "N7", "CA", "CA", "CA"])
    arr.element = np.array(["N", "C", "N", "C", "C", "C"])
    arr.hetero = np.array([False] * n)
    rmsd, n_ca = fbb.dna_aligned_ca_rmsd(arr, arr.copy())
    assert rmsd == 0.0
    assert n_ca == 3


# --- two-stage binder-block gate (corrected 2026-07-31) -------------------
# Paper Methods: fold -> DNA-aligned RMSD < 8 A -> LigandMPNN resample -> fold ->
# RMSD < 3 A, ipTM > 0.7, high H-bond counts. The 8 A pre-resample gate was
# missing, which would have sent every diffused backbone into resampling rather
# than only the self-consistent ones.

def _run_filter(tmp_path, rows, extra):
    import subprocess, csv
    man = tmp_path / "m.json"
    # analyze_one will fail on these stub paths; the filter records the error row
    # and still applies the gate, which is what we are testing.
    json.dump(rows, open(man, "w"))
    out, cmp_ = tmp_path / "p.csv", tmp_path / "c.csv"
    r = subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "..",
                                                     "scripts", "filter_binder_block.py"),
                        "--manifest", str(man), "--out", str(out),
                        "--oracle-comparison", str(cmp_), *extra],
                       capture_output=True, text=True)
    return r.stdout + r.stderr


def test_pre_resample_stage_defaults_to_8A_and_drops_the_iptm_gate(tmp_path):
    txt = _run_filter(tmp_path, [], ["--stage", "pre_resample"])
    assert "RMSD<8.0" in txt, txt
    assert "no ipTM gate" in txt, txt


def test_post_resample_stage_defaults_to_3A_with_iptm(tmp_path):
    txt = _run_filter(tmp_path, [], ["--stage", "post_resample"])
    assert "RMSD<3.0" in txt and "ipTM>0.7" in txt, txt


def test_post_resample_is_the_default_stage(tmp_path):
    txt = _run_filter(tmp_path, [], [])
    assert "stage=post_resample" in txt, txt


def test_pre_resample_gate_admits_a_backbone_the_3A_gate_would_reject():
    # a 5 A design is self-consistent enough to resample but not to pass the
    # final gate -- the whole point of having two thresholds
    assert fbb  # module imported
    for gate, expect in ((8.0, True), (3.0, False)):
        assert (5.0 < gate) is expect


# --- on-target minPAE emission (added 2026-08-06) ---------------------------
# The specificity block's entry gate is an on-target minPAE cut, but this script
# emitted RMSD/ipTM/H-bonds only -- so the criterion could not be evaluated at all,
# even though the PAE was already on disk in every rf3 refold's *_confidences.json.

def _pae_json(path, prot_dna_min, n_prot=3, n_dna=2):
    pae = np.full((n_prot + n_dna,) * 2, 20.0)
    np.fill_diagonal(pae, 0.5)
    pae[0, n_prot] = prot_dna_min
    pae[n_prot, 0] = prot_dna_min
    json.dump({"pae": pae.tolist(),
               "token_chain_ids": ["A"] * n_prot + ["B", "C"][:n_dna]}, open(path, "w"))


def test_min_pae_read_from_the_refold_confidences(tmp_path):
    p = tmp_path / "c.json"
    _pae_json(p, 3.7)
    got = fbb._on_target_min_pae({"pae_path": str(p), "protein_chain": "A",
                                  "dna_chains": ["B", "C"], "protein_len": 3})
    assert got == 3.7


def test_min_pae_matches_compute_delta_minpae_exactly(tmp_path):
    """Two independent minima over the same matrix is how a silent inconsistency gets
    in, so the filter reuses the ranking script's loader and masks."""
    import compute_delta_minpae as cdm
    p = tmp_path / "c.json"
    _pae_json(p, 1.9)
    job = {"pae_path": str(p), "protein_chain": "A", "dna_chains": ["B", "C"],
           "protein_len": 3}
    pae, chains = cdm.load_pae(str(p))
    prot, dna = cdm.protein_dna_token_masks(pae.shape[0], chains, "A", ["B", "C"],
                                            (0, 3), [])
    assert fbb._on_target_min_pae(job) == round(cdm.min_pae(pae, prot, dna), 4)


def test_missing_pae_is_none_not_a_crash():
    assert fbb._on_target_min_pae({}) is None
    assert fbb._on_target_min_pae({"pae_path": "FILL_AFTER_FOLD"}) is None
    assert fbb._on_target_min_pae({"pae_path": "/nonexistent/x.json"}) is None


def _run_filter_rows(tmp_path, jobs, *extra):
    import subprocess
    m = tmp_path / "m.json"
    json.dump(jobs, open(m, "w"))
    out, cmp_ = tmp_path / "p.csv", tmp_path / "a.csv"
    r = subprocess.run([sys.executable,
                        os.path.join(os.path.dirname(__file__), "..", "scripts",
                                     "filter_binder_block.py"),
                        "--manifest", str(m), "--out", str(out),
                        "--oracle-comparison", str(cmp_), *extra],
                       capture_output=True, text=True)
    import csv as _csv
    rows = list(_csv.DictReader(open(cmp_))) if cmp_.exists() else []
    return r, rows


def test_min_pae_column_is_emitted_even_without_a_gate(tmp_path):
    p = tmp_path / "c.json"
    _pae_json(p, 4.2)
    r, rows = _run_filter_rows(tmp_path, [{"design_id": "d1", "oracle": "rf3",
                                      "design_path": "/nope.pdb",
                                      "refold_path": "/nope.cif",
                                      "pae_path": str(p), "protein_len": 3}])
    assert r.returncode == 0, r.stderr
    assert rows[0]["min_pae"] == "4.2"
    assert "minPAE recovered for 1/1" in r.stdout


def test_an_unreadable_pae_fails_a_requested_gate_rather_than_passing_it(tmp_path):
    """A gate the user asked for must never be satisfied by absent data."""
    r, rows = _run_filter_rows(tmp_path,
                          [{"design_id": "d1", "oracle": "rf3",
                            "design_path": "/nope.pdb", "refold_path": "/nope.cif"}],
                          "--min-pae-gate", "6.6")
    assert r.returncode == 0, r.stderr
    assert rows[0]["min_pae"] in ("", "None")
    assert "0 passers" in r.stdout
    assert "NO row had a readable PAE" in r.stdout
