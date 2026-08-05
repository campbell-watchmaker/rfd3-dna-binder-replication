"""CI tests for the specificity-block scripts (self-contained, no network/GPU).

Cover the off-target panel construction and the ΔminPAE math -- the two places a
silent regression would corrupt the specificity ranking.
"""
import json
import os
import sys

import numpy as np

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
