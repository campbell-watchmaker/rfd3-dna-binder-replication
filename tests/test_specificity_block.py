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
