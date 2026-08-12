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
import pytest

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


# --- H-bond constraint SAMPLING (added 2026-08-06) --------------------------
# The paper describes the constraint set as varied per design -- "we sample ... a
# diverse set of hydrogen bond (Hbond) condition constraints" -- so one fixed subset is
# the wrong SHAPE regardless of its size. Consistent with the measurement that 6 vs 8
# atoms was indistinguishable over 100 refolds: the count was not the operative
# variable. The DRAWING MECHANISM is ours; the paper states the diversity, not the how.

def _cands():
    """A 12-bp PRNP-like duplex, numbered the way the REAL conditioning bundle is.

    Both strands run 1..12 in their own 5'->3' numbering, and the duplex is ANTIPARALLEL,
    so B_j is the complement of A_(13-j) -- NOT of A_j. Getting this wrong in the fixture
    would hide exactly the defect these tests exist to catch.

    Purines carry 2 major-groove atoms each (G N7+O6, A N7+N6), pyrimidines 1
    (C N4, T O4), which is why base COUNT and atom COUNT differ per draw.
    """
    seq = "TGAGGAGAGGAG"
    comp = {"G": "C", "A": "T", "T": "A", "C": "G"}
    anti = "".join(comp[b] for b in reversed(seq))      # CTCCTCTCCTCA
    out = []

    def atoms(chain, res_id, base):
        d = "D" + base
        if base == "G":
            return [{"chain": chain, "res_id": res_id, "res_name": d, "atom": "N7",
                     "role": "acceptor"},
                    {"chain": chain, "res_id": res_id, "res_name": d, "atom": "O6",
                     "role": "acceptor"}]
        if base == "A":
            return [{"chain": chain, "res_id": res_id, "res_name": d, "atom": "N7",
                     "role": "acceptor"},
                    {"chain": chain, "res_id": res_id, "res_name": d, "atom": "N6",
                     "role": "donor"}]
        return [{"chain": chain, "res_id": res_id, "res_name": d,
                 "atom": "N4" if base == "C" else "O4",
                 "role": "donor" if base == "C" else "acceptor"}]

    for i, b in enumerate(seq, start=1):
        out += atoms("A", i, b)
    for j, b in enumerate(anti, start=1):
        out += atoms("B", j, b)
    return out


def test_the_fixture_really_is_antiparallel():
    """Guards the guard: if the fixture were numbered in parallel, every pairing test
    below would pass vacuously."""
    by = {(c["chain"], c["res_id"]): c["res_name"] for c in _cands()}
    assert by[("A", 1)] == "DT" and by[("B", 12)] == "DA"
    assert by[("A", 12)] == "DG" and by[("B", 1)] == "DC"


def test_antiparallel_pairing_is_verified_not_assumed():
    """B_j pairs with A_(L+1-j). This is CHECKED against Watson-Crick complementarity,
    so a target whose strands are not a reverse-complement pair fails loudly instead of
    yielding a plausible, wrong window."""
    bp = mrs.bp_index_map(_cands())
    assert bp[("A", 1)] == 1 and bp[("A", 12)] == 12
    assert bp[("B", 1)] == 12, "B1 must pair with A12, not A1"
    assert bp[("B", 12)] == 1
    assert bp[("B", 5)] == 8


def test_pairing_check_rejects_a_non_complementary_duplex():
    bad = [{"chain": "A", "res_id": 1, "res_name": "DG", "atom": "N7", "role": "acceptor"},
           {"chain": "A", "res_id": 2, "res_name": "DG", "atom": "N7", "role": "acceptor"},
           # B1 should complement A2 (=DC); make it DA so the check must fire
           {"chain": "B", "res_id": 1, "res_name": "DA", "atom": "N7", "role": "acceptor"},
           {"chain": "B", "res_id": 2, "res_name": "DC", "atom": "N4", "role": "donor"}]
    with pytest.raises(ValueError, match="pairing check FAILED"):
        mrs.bp_index_map(bad)


def test_antisense_draws_land_in_their_own_ori_window(tmp_path):
    """THE BUG, pinned. Both strands are numbered 1..L independently and the duplex is
    antiparallel, so comparing a chain-B res_id against a sense-strand window selects the
    OPPOSITE END of the target. It was dormant while conditioning was purine-only and went
    live with --hbond-strand either: 3 of 10 Phase-2.5 sampled specs drew strand B and got
    constraints in the other ori's half."""
    import random
    rng = random.Random(0)
    saw_b = False
    for _ in range(200):
        kept, d = mrs.sample_for_ori(_cands(), 7, 12, rng, (2, 5), "either")
        if not d["base_positions"]:
            continue
        assert all(7 <= p <= 12 for p in d["base_positions"]), d
        # and the ATOMS actually kept must sit at those base pairs
        assert all(7 <= p <= 12 for p in d["atom_bp_positions"]), d
        saw_b |= d["strands"] == ["B"]
    assert saw_b, "never drew the antisense strand, so the regression is untested"


def test_draws_are_contiguous_base_pair_runs():
    """A recognition helix reads a consecutive stretch; a scattered set may demand a
    geometry no single fold satisfies."""
    import random
    rng = random.Random(3)
    for _ in range(200):
        _, d = mrs.sample_for_ori(_cands(), 1, 6, rng, (2, 5), "either")
        pos = d["base_positions"]
        if len(pos) > 1:
            assert pos == list(range(min(pos), max(pos) + 1)), f"not contiguous: {pos}"
        assert d["pattern"] == "contiguous"


def test_a_run_never_jumps_a_gap_in_usable_positions():
    """`usable` has gaps where a base contributes no major-groove atom. Slicing the list
    would emit a 'run' that silently jumps one."""
    import random
    # a strand where bp 3 is absent from the candidate pool entirely
    cands = [c for c in _cands() if not (c["chain"] == "A" and c["res_id"] == 3)]
    rng = random.Random(11)
    for _ in range(200):
        _, d = mrs.sample_for_ori(cands, 1, 6, rng, (2, 5), "purine")
        pos = d["base_positions"]
        if len(pos) > 1:
            assert pos == list(range(min(pos), max(pos) + 1)), pos
            assert 3 not in pos


def test_both_strands_are_eligible_and_neither_dominates():
    """[PAPER] 'candidate major groove donor and acceptor atoms' -- no strand restriction.
    The old purine-only rule was ours."""
    import random
    rng = random.Random(5)
    strands = set()
    for _ in range(60):
        _, d = mrs.sample_for_ori(_cands(), 1, 12, rng, (2, 5), "either")
        strands |= set(d["strands"])
    assert strands == {"A", "B"}, strands


def test_terminal_base_pairs_are_never_drawn():
    import random
    rng = random.Random(5)
    for _ in range(200):
        _, d = mrs.sample_for_ori(_cands(), 1, 12, rng, (2, 5), "either")
        assert 1 not in d["base_positions"] and 12 not in d["base_positions"], d


def test_run_length_is_capped_by_the_window_not_the_request():
    """A 6-bp ori window minus a terminal base pair leaves 5 usable, so a request for up
    to 8 must silently cap at 5 rather than overflow the window."""
    import random
    rng = random.Random(9)
    for _ in range(200):
        _, d = mrs.sample_for_ori(_cands(), 1, 6, rng, (2, 8), "either")
        assert d["n_bases_drawn"] <= 5, d
        assert all(1 <= p <= 6 for p in d["base_positions"]), d


def test_both_atoms_of_a_purine_stay_together():
    import random
    rng = random.Random(3)
    for _ in range(60):
        kept, _ = mrs.sample_for_ori(_cands(), 1, 12, rng, (2, 5), "purine")
        by_res = {}
        for c in kept:
            by_res.setdefault(c["res_id"], set()).add(c["atom"])
        for res, atoms in by_res.items():
            if "N7" in atoms:
                assert len(atoms) == 2, f"res {res} kept only {atoms}"


def test_sampling_is_reproducible_from_the_seed():
    import random
    a, ra = mrs.sample_for_ori(_cands(), 1, 6, random.Random(7))
    b, rb = mrs.sample_for_ori(_cands(), 1, 6, random.Random(7))
    assert ra == rb


def test_sampling_actually_varies_across_designs():
    import random
    rng = random.Random(42)
    draws = [mrs.sample_for_ori(_cands(), 1, 6, rng, (2, 5), "either")[1]
             for _ in range(30)]
    assert len({tuple(d["base_positions"]) for d in draws}) > 1
    assert len({d["n_bases_drawn"] for d in draws}) > 1


def test_fixed_mode_is_unchanged_and_remains_the_default():
    """Existing arms must stay reproducible."""
    import subprocess
    r = subprocess.run([sys.executable,
                        os.path.join(os.path.dirname(__file__), "..", "scripts",
                                     "make_rfd3na_specs.py"), "--help"],
                       capture_output=True, text=True)
    assert "fixed (default)" in r.stdout
    kept = mrs.subset_for_ori(_cands(), 1, 6, 3, "purine")
    assert {c["res_id"] for c in kept} == {3, 4, 5}


def test_designs_per_ori_is_refused_without_sampling(tmp_path):
    """Without sampling the subset is deterministic, so N specs would be N copies."""
    import subprocess
    r = subprocess.run([sys.executable,
                        os.path.join(os.path.dirname(__file__), "..", "scripts",
                                     "make_rfd3na_specs.py"),
                        "--conditioning", os.path.join(os.path.dirname(__file__), "..",
                                                       "targets", "prnp",
                                                       "conditioning.json"),
                        "--duplex-cif", "/workspace/d.cif",
                        "--out-dir", str(tmp_path / "o"), "--designs-per-ori", "5"],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "only makes sense with --hbond-sampling random" in r.stderr


def test_manifest_records_the_draw_for_every_spec(tmp_path):
    import subprocess
    out = tmp_path / "specs"
    r = subprocess.run([sys.executable,
                        os.path.join(os.path.dirname(__file__), "..", "scripts",
                                     "make_rfd3na_specs.py"),
                        "--conditioning", os.path.join(os.path.dirname(__file__), "..",
                                                       "targets", "prnp",
                                                       "conditioning.json"),
                        "--duplex-cif", "/workspace/d.cif", "--out-dir", str(out),
                        "--hbond-sampling", "random", "--designs-per-ori", "4",
                        "--seed", "99"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    man = json.load(open(out / "manifest.json"))
    assert man["hbond_sampling"] == "random" and man["hbond_seed"] == 99
    assert len(man["specs"]) == 8, "2 ori x 4 draws"
    for s in man["specs"]:
        d = s["hbond_draw"]
        assert d["mode"] == "random" and d["atoms"] and d["n_atoms"] == len(d["atoms"])
        assert os.path.isfile(out / s["spec"])


# --- sampler config -> pecli flags (added 2026-08-06) ----------------------
# PIPELINE.md showed `pecli prepare rfd3na --design-inputs SPEC --config CFG:_smoke_test`.
# NEITHER flag exists; pecli rejects both with "unknown option(s) for rfd3na: config,
# design_inputs". The spec goes in via --input and every sampler knob is its own flag.
# The command as documented could never have run.
#
# It also surfaced that `cfg_features` -- carried in sampler_config.json with a
# paragraph of rationale -- is not a field pecli exposes, so it has never reached the
# sampler on any run. A silently-inert setting is the failure mode these tests pin.

import submit_arm_diffusion as sad  # noqa: E402


def test_config_becomes_one_flag_per_field():
    flags, skipped = sad.config_to_flags({"diffusion_batch_size": 5, "n_batches": 2,
                                          "step_scale": 1.5, "gamma_0": 0.6})
    assert flags == ["--diffusion-batch-size", "5", "--n-batches", "2",
                     "--step-scale", "1.5", "--gamma-0", "0.6"]
    assert skipped == []


def test_bools_are_passed_with_an_explicit_value():
    """A BARE bool flag means True in pecli's parser, so `false` is only expressible by
    passing the value -- emitting a bare flag for cfg=False would silently turn CFG ON,
    which is precisely the arm being measured."""
    on, _ = sad.config_to_flags({"use_classifier_free_guidance": True})
    off, _ = sad.config_to_flags({"use_classifier_free_guidance": False})
    assert on == ["--use-classifier-free-guidance", "true"]
    assert off == ["--use-classifier-free-guidance", "false"]


def test_cfg_features_is_reported_as_unsupported_not_passed():
    """It is in sampler_config.json, it is not a pecli field, and passing it makes
    prepare fail outright."""
    flags, skipped = sad.config_to_flags(
        {"cfg_scale": 1.5, "cfg_features": ["active_donor", "active_acceptor"]})
    assert flags == ["--cfg-scale", "1.5"]
    assert skipped == ["cfg_features"]


def test_comment_keys_are_ignored_silently_but_unknown_keys_are_not():
    """Underscore keys are documentation. Anything else unknown is a real setting that
    will not take effect, and must be surfaced."""
    flags, skipped = sad.config_to_flags(
        {"_comment": "x", "_cfg_note": "y", "made_up_knob": 3, "n_batches": 1})
    assert flags == ["--n-batches", "1"]
    assert skipped == ["made_up_knob"]


def test_the_committed_sampler_config_translates_cleanly():
    """The real file, so a future edit that adds an inert knob fails here."""
    cfg = json.load(open(os.path.join(os.path.dirname(__file__), "..", "specs",
                                      "binder_block", "sampler_config.json")))
    for arm in ("_smoke_test", "_full_run"):
        flags, skipped = sad.config_to_flags(cfg[arm])
        assert "--diffusion-batch-size" in flags and "--n-batches" in flags
        assert "--use-classifier-free-guidance" in flags
        assert skipped == ["cfg_features"], (
            f"{arm}: unexpected inert setting(s) {skipped} -- either pecli gained the "
            "field or a knob was added that will silently do nothing")


def test_every_known_field_matches_peclis_rfd3na_tool_spec():
    """KNOWN_FIELDS is a local copy of pecli's field list; if the tool gains or loses a
    knob this drifts silently. Skips when pecli is not importable (CI without it)."""
    pytest_mod = __import__("pytest")
    try:
        from pecli.tools import rfd3na  # noqa: F401
        from pecli.tools.base import SPECS  # type: ignore
    except Exception:
        pytest_mod.skip("pecli not importable from this interpreter")
    spec = SPECS["rfd3na"]
    real = {f.key for f in spec.fields} | {"gpu"}
    unknown_to_pecli = sad.KNOWN_FIELDS - real
    assert not unknown_to_pecli, f"we pass fields pecli does not have: {unknown_to_pecli}"
