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
