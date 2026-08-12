"""CI tests for the specificity-block scripts (self-contained, no network/GPU).

Cover the off-target panel construction and the ΔminPAE math -- the two places a
silent regression would corrupt the specificity ranking.
"""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_offtarget_set as mos
import compute_delta_minpae as cdm


def test_single_base_variants_count_and_content():
    seq = "TGAGGAGAGGAG"  # 12 bp
    variants = mos.single_base_variants(seq)
    assert len(variants) == 3 * len(seq)  # 3 alternatives per position
    # every variant differs from WT at exactly one position
    for name, var in variants:
        assert len(var) == len(seq)
        diffs = [i for i in range(len(seq)) if var[i] != seq[i]]
        assert len(diffs) == 1


def test_minpae_uses_both_orientations_and_global_min():
    # 5x5: protein tokens 0-2 (chain A), DNA tokens 3-4 (B,C). Seed the min in the
    # DNA->protein orientation only, to prove both orientations are checked.
    pae = np.full((5, 5), 20.0)
    np.fill_diagonal(pae, 0.5)
    pae[4, 0] = 1.3  # DNA token 4 vs protein token 0
    chains = np.array(["A", "A", "A", "B", "C"])
    prot = chains == "A"
    dna = np.isin(chains, ["B", "C"])
    assert cdm.min_pae(pae, prot, dna) == 1.3


def test_delta_minpae_ranks_specific_above_promiscuous(tmp_path):
    def write_pae(path, prot_dna_min):
        pae = np.full((5, 5), 20.0)
        np.fill_diagonal(pae, 0.5)
        pae[1, 3] = prot_dna_min
        pae[3, 1] = prot_dna_min
        json.dump({"pae": pae.tolist(), "token_chain_ids": ["A", "A", "A", "B", "C"]}, open(path, "w"))

    jobs = []
    # specific: on low, offs high
    write_pae(tmp_path / "s_on.json", 2.0)
    write_pae(tmp_path / "s_off.json", 15.0)
    jobs += [
        {"design_id": "spec", "dna_id": "on_target", "kind": "on_target", "pae_path": str(tmp_path / "s_on.json"), "oracle": "protenix"},
        {"design_id": "spec", "dna_id": "v1", "kind": "sbs", "pae_path": str(tmp_path / "s_off.json"), "oracle": "protenix"},
    ]
    # promiscuous: on low, an off also low
    write_pae(tmp_path / "p_on.json", 2.0)
    write_pae(tmp_path / "p_off.json", 2.4)
    jobs += [
        {"design_id": "prom", "dna_id": "on_target", "kind": "on_target", "pae_path": str(tmp_path / "p_on.json"), "oracle": "protenix"},
        {"design_id": "prom", "dna_id": "v1", "kind": "sbs", "pae_path": str(tmp_path / "p_off.json"), "oracle": "protenix"},
    ]
    mpath = tmp_path / "m.json"
    json.dump(jobs, open(mpath, "w"))
    out = tmp_path / "delta.csv"
    import subprocess
    script = os.path.join(os.path.dirname(__file__), "..", "scripts", "compute_delta_minpae.py")
    subprocess.run([sys.executable, script, "--manifest", str(mpath), "--out", str(out)], check=True)
    rows = list(csv_dicts(out))
    assert rows[0]["design_id"] == "spec"      # specific ranks first
    assert float(rows[0]["delta_min_pae"]) > float(rows[1]["delta_min_pae"])


def csv_dicts(path):
    import csv
    with open(path) as f:
        yield from csv.DictReader(f)


# --- panel-separation regression (corrected 2026-07-31) --------------------
# The paper's ΔminPAE all-by-all folds against on-target + other Table 1 targets.
# The single-base-variant sweep is wet-lab characterisation of an already-selected
# binder, NOT the ranking panel. Conflating them inflated the all-by-all ~4x AND
# silently changed the metric: ΔminPAE is a MINIMUM over off-targets, so including
# sequences one base from the on-target turns it into a near-worst-case statistic.

def _panel(tmp_path, panel):
    import subprocess
    out = tmp_path / f"off_{panel}.json"
    subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "..",
                                                 "scripts", "make_offtarget_set.py"),
                    "--on-target", "TGAGGAGAGGAG", "--panel", panel, "--out", str(out)],
                   check=True, capture_output=True)
    return json.load(open(out))


def test_ranking_panel_excludes_single_base_variants(tmp_path):
    b = _panel(tmp_path, "ranking")
    kinds = {e["kind"] for e in b["offtargets"]}
    assert "sbs" not in kinds, "single-base variants must not be in the ranking panel"
    assert kinds == {"on_target", "decoy"}
    assert b["n_sbs"] == 0 and b["n_decoys"] > 0


def test_ranking_panel_is_the_default(tmp_path):
    import subprocess
    out = tmp_path / "default.json"
    subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "..",
                                                 "scripts", "make_offtarget_set.py"),
                    "--on-target", "TGAGGAGAGGAG", "--out", str(out)],
                   check=True, capture_output=True)
    assert json.load(open(out))["panel"] == "ranking"


def test_sbs_panel_still_available_and_complete(tmp_path):
    b = _panel(tmp_path, "sbs")
    sbs = [e for e in b["offtargets"] if e["kind"] == "sbs"]
    assert len(sbs) == 3 * 12, "sbs panel must still cover every position x every alt base"
    assert not [e for e in b["offtargets"] if e["kind"] == "decoy"]


def test_ranking_panel_is_much_cheaper_than_the_old_conflated_one(tmp_path):
    rank = _panel(tmp_path, "ranking")
    both = _panel(tmp_path, "both")
    assert len(rank["offtargets"]) < len(both["offtargets"]) / 3, (
        "the corrected ranking panel should be several-fold smaller than the "
        "old decoys+sbs union")


# --- templating is the default for the all-by-all (decided 2026-08-05) -----
# The paper templates the protein chain in exactly one place -- "Templates were not
# used throughout the design campaign with the exception of the all-by-all folding
# step in the specificity block" -- and our own control panel measured it raising
# ΔminPAE (LambdaRep +1.93 -> +3.43, Engrailed +0.26 -> +0.72) with argmin holding
# 3/3 and every on-target interface still inside the motif window. These tests pin
# that default, and pin the two ways it could be got wrong: templating the DNA (which
# would hand the fold the docking geometry it is meant to predict), and carrying
# templating over to the binder block's self-consistency gate (which would hand that
# fold the very backbone it is being asked to independently reproduce).

def _aba(tmp_path, extra, seq="MKTAYIAKQRQISFVKSHFSRQ"):
    import subprocess
    off = tmp_path / "off.json"
    subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "..",
                                                 "scripts", "make_offtarget_set.py"),
                    "--on-target", "TGAGGAGAGGAG", "--panel", "ranking",
                    "--out", str(off)], check=True, capture_output=True)
    fa = tmp_path / "d.fasta"
    fa.write_text(f">d1\n{seq}\n")
    out = tmp_path / ("aba" + str(abs(hash(tuple(extra))) % 9999))
    r = subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "..",
                                                     "scripts", "build_allbyall_inputs.py"),
                        "--design-fasta", str(fa), "--offtargets", str(off),
                        "--out-dir", str(out), *extra],
                       capture_output=True, text=True)
    return r, out


def test_allbyall_requires_a_template_dir_by_default(tmp_path):
    """Templating is not opt-in: omitting it must be an error, not a silent
    fall-through to sequence-only specs."""
    r, _ = _aba(tmp_path, [])
    assert r.returncode != 0
    assert "--template-dir is required" in (r.stdout + r.stderr)


def test_allbyall_templates_the_protein_and_never_the_dna(tmp_path):
    cif = tmp_path / "tmpl"
    cif.mkdir()
    (cif / "d1_template.cif").write_text("data_stub\n")
    r, out = _aba(tmp_path, ["--template-dir", str(cif)])
    assert r.returncode == 0, r.stdout + r.stderr
    spec = json.load(open(out / "d1__on_target.json"))
    entry = spec[0]
    assert entry["template_selection"] == ["A"], "only the protein chain may be templated"
    comps = entry["components"]
    assert "path" in comps[0], "the protein should ride as a template CIF component"
    dna = [c for c in comps if c.get("chain_type") == "polydeoxyribonucleotide"]
    assert len(dna) == 2, "both DNA strands must be present as free sequence components"
    assert all("path" not in c for c in dna), "the DNA must never be templated"


def test_allbyall_targets_rf3_not_protenix(tmp_path):
    """rf3 won the oracle comparison (argmin 4/5 vs protenix 2/5, whose binder /
    non-binder ranges overlap) and is ~2x cheaper."""
    cif = tmp_path / "tmpl2"
    cif.mkdir()
    (cif / "d1_template.cif").write_text("data_stub\n")
    _, out = _aba(tmp_path, ["--template-dir", str(cif)])
    man = json.load(open(out / "folds_manifest.json"))
    assert man, "manifest is empty"
    assert {m["oracle"] for m in man} == {"rf3"}
    assert all(m["templated"] for m in man)


def test_allbyall_untemplated_escape_hatch_is_msa_free(tmp_path):
    r, out = _aba(tmp_path, ["--no-template"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert "NOT the paper's protocol" in r.stdout
    spec = json.load(open(out / "d1__on_target.json"))
    p = spec[0]["components"][0]
    assert p["_pecli_rf3_msa_a3m"].count(">") == 1
    assert "template_selection" not in spec[0]


def test_binder_block_rmsd_warns_against_templating():
    """The binder-block gate must stay untemplated; keep the reason in the source."""
    src = open(os.path.join(os.path.dirname(__file__), "..", "scripts",
                            "filter_binder_block.py")).read()
    assert "MUST BE UNTEMPLATED" in src


# --- Stage 4 -> Stage 5 handoff (added 2026-08-06) --------------------------
# Every other ΔminPAE test hand-builds its manifest, and every one of them omits
# protein_len. That is exactly why a GUARANTEED crash went unnoticed:
# compute_delta_minpae.py did `tuple(j["protein_len"])`, and build_allbyall_inputs.py
# emits protein_len as an int, so `tuple(123)` raised TypeError on record 1 of every
# real manifest. Because Python evaluates call arguments eagerly it fired even when
# chain labels were present and the range was never consulted.
#
# These tests therefore drive the REAL emitter and feed its REAL output to the
# consumer. Any future divergence between the two shapes fails here.

def _emit_real_manifest(tmp_path):
    """Run build_allbyall_inputs.py for real and return its folds_manifest.json."""
    import subprocess
    off = tmp_path / "off.json"
    subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "..",
                                                 "scripts", "make_offtarget_set.py"),
                    "--on-target", "TGAGGAGAGGAG", "--panel", "ranking",
                    "--out", str(off)], check=True, capture_output=True)
    fa = tmp_path / "d.fasta"
    fa.write_text(">d1\nMKTAYIAKQRQISFVKSHFSRQ\n")
    out = tmp_path / "aba_real"
    r = subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "..",
                                                     "scripts", "build_allbyall_inputs.py"),
                        "--design-fasta", str(fa), "--offtargets", str(off),
                        "--out-dir", str(out), "--no-template"],
                       check=True, capture_output=True, text=True)
    return json.load(open(out / "folds_manifest.json")), out


def test_emitted_manifest_protein_len_is_accepted_by_the_consumer(tmp_path):
    """The crash, pinned at its source: the emitter's protein_len must survive the
    consumer's normalisation."""
    man, _ = _emit_real_manifest(tmp_path)
    assert man, "emitter produced no manifest"
    assert isinstance(man[0]["protein_len"], int), \
        "emitter changed protein_len's type; update _protein_len_range too"
    lo, hi = cdm._protein_len_range(man[0])
    assert (lo, hi) == (0, len("MKTAYIAKQRQISFVKSHFSRQ")), (lo, hi)


def test_protein_len_range_accepts_both_repo_conventions():
    """Scalar count (the majority convention) and explicit [lo, hi] range."""
    assert cdm._protein_len_range({"protein_len": 123}) == (0, 123)
    assert cdm._protein_len_range({"protein_len": [0, 123]}) == (0, 123)
    assert cdm._protein_len_range({"protein_len": 86, "protein_copies": 2}) == (0, 172)
    assert cdm._protein_len_range({}) == (0, 0)
    with pytest.raises(ValueError):
        cdm._protein_len_range({"protein_len": [1, 2, 3]})


def test_every_consumer_required_field_is_emitted(tmp_path):
    """Field-name drift between the two scripts is the standing risk here."""
    man, _ = _emit_real_manifest(tmp_path)
    for rec in man:
        for field in ("design_id", "dna_id", "kind", "pae_path", "oracle",
                      "protein_chain", "dna_chains", "protein_len"):
            assert field in rec, f"{field} missing from emitted manifest"
        assert rec["kind"] in ("on_target", "sbs", "decoy"), rec["kind"]
    assert sum(1 for r in man if r["kind"] == "on_target") == 1, \
        "exactly one on-target row per design is required for ΔminPAE"


def test_unrankable_designs_are_reported_not_silently_dropped(tmp_path):
    """A design whose on-target fold is missing must be REPORTED. Silently dropping it
    yields a shorter, entirely credible, wrong ranking on a partially drained batch."""
    import subprocess
    import numpy as np
    # one off-target fold only -- no on-target
    pae = np.full((6, 6), 20.0)
    p = tmp_path / "off.json"
    json.dump({"pae": pae.tolist(), "token_chain_ids": ["A"] * 3 + ["B"] * 3}, open(p, "w"))
    man = [{"design_id": "d1", "dna_id": "decoy_TBP", "kind": "decoy",
            "oracle": "rf3", "pae_path": str(p), "protein_chain": "A",
            "dna_chains": ["B"], "protein_len": 3}]
    mpath = tmp_path / "m.json"
    json.dump(man, open(mpath, "w"))
    r = subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "..",
                                                     "scripts", "compute_delta_minpae.py"),
                        "--manifest", str(mpath), "--out", str(tmp_path / "o.csv")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "UNRANKABLE" in r.stdout, r.stdout
    assert "no on-target fold" in r.stdout, r.stdout


# --- off-target padding (added 2026-08-06) ---------------------------------
# 4 of the 9 Table 1 decoys are 10 bp against a 12-bp on-target. minPAE is a MINIMUM
# over protein x DNA token pairs, so a shorter duplex offers fewer pairs to minimise
# over and loses as an off-target for a reason unrelated to specificity, biasing
# ΔminPAE upward. The old code recorded this as `same_length_as_on` and nothing ever
# read the field.

def _build_panel(tmp_path, *extra):
    import subprocess
    out = tmp_path / "off.json"
    r = subprocess.run([sys.executable,
                        os.path.join(os.path.dirname(__file__), "..", "scripts",
                                     "make_offtarget_set.py"),
                        "--on-target", "TGAGGAGAGGAG", "--out", str(out), *extra],
                       capture_output=True, text=True)
    return r, (json.load(open(out)) if out.exists() else None)


def test_every_duplex_is_one_length_by_default(tmp_path):
    r, b = _build_panel(tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    lens = {len(e["sense"]) for e in b["offtargets"]}
    assert lens == {12}, lens
    assert {len(e["antisense"]) for e in b["offtargets"]} == {12}
    assert b["padded_to_bp"] == 12


def test_padding_preserves_the_motif_and_records_it(tmp_path):
    _, b = _build_panel(tmp_path)
    by_id = {e["id"]: e for e in b["offtargets"]}
    p53 = by_id["decoy_P53"]
    assert p53["motif"] == "AGACATGTCT"
    assert p53["motif"] in p53["sense"], "the actual site must survive padding"
    assert p53["sense"] == p53["left_pad"] + p53["motif"] + p53["right_pad"]
    # the on-target is already at the panel length, so it must be untouched
    assert by_id["on_target"]["sense"] == "TGAGGAGAGGAG"
    assert by_id["on_target"]["left_pad"] == ""


def test_no_pad_is_available_and_leaves_native_lengths(tmp_path):
    _, b = _build_panel(tmp_path, "--no-pad")
    assert b["padded_to_bp"] is None
    assert {len(e["sense"]) for e in b["offtargets"]} == {10, 12}


def test_pad_target_defaults_to_the_longest_so_nothing_is_truncated(tmp_path):
    """build_duplex() raises on a motif longer than fixed_bp; a hardcoded 24 would be
    fine here but would silently break on any panel with a longer site."""
    _, b = _build_panel(tmp_path)
    assert b["padded_to_bp"] == max(len(e["motif"]) for e in b["offtargets"])


def test_padding_is_refused_if_a_flank_introduces_another_motif(tmp_path):
    """The check that makes padding safe: verify_panel() scans every padded duplex,
    both strands, for every other target's motif."""
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis",
                                    "oracle_controls"))
    entries = [
        {"id": "on_target", "kind": "on_target", "sense": "TGAGGAGAGGAG",
         "antisense": "CTCCTCTCCTCA"},
        # a decoy whose padded form will contain the flank-derived motif below
        {"id": "decoy_X", "kind": "decoy", "sense": "AAAA", "antisense": "TTTT"},
        # ... and a "motif" that is a substring of the neutral flank, so padding
        # decoy_X out to 12 bp pulls it in
        {"id": "decoy_flankish", "kind": "decoy", "sense": "CTGACTTG",
         "antisense": "CAAGTCAG"},
    ]
    _, problems = mos.pad_entries(entries, 12)
    assert problems, "a flank-derived motif collision must be reported"
    assert any("decoy_flankish" in p for p in problems), problems
