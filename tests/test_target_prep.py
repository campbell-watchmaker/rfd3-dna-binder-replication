"""CI smoke test for scripts/target_prep.py.

Fetches the canonical B-DNA reference (Drew-Dickerson dodecamer, PDB 1BNA) and
checks the chemical self-consistency of the ori-token / major-groove geometry:
purine N7 atoms must fall on the computed major-groove side, and the mean
interior helical twist must match canonical B-DNA (~34 deg/bp). A regression
in the groove-direction math would silently condition rfd3na on the wrong
side of the DNA, so this is a correctness gate, not a style check.
"""
import os
import sys
import urllib.request

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from target_prep import load_duplex, base_pairs, ori_tokens, hbond_candidates, validate_major_groove

REF_PDB_URL = "https://files.rcsb.org/download/1BNA.cif"


@pytest.fixture(scope="module")
def duplex(tmp_path_factory):
    dst = tmp_path_factory.mktemp("ref") / "1BNA.cif"
    urllib.request.urlretrieve(REF_PDB_URL, dst)
    return load_duplex(str(dst))


def test_major_groove_direction_is_chemically_correct(duplex):
    frac_n7_major, mean_twist = validate_major_groove(duplex)
    assert frac_n7_major >= 0.95, (
        f"only {frac_n7_major:.0%} of purine N7 atoms on the major-groove side "
        "-- groove-direction sign is likely flipped"
    )
    assert 28.0 <= mean_twist <= 40.0, (
        f"mean interior helical twist {mean_twist:.1f} deg/bp is outside the "
        "canonical B-DNA range (~34 deg/bp) -- base-pair frame construction is off"
    )


def test_ori_tokens_cover_all_base_pairs(duplex):
    pairs = base_pairs(duplex)
    tokens = ori_tokens(duplex, window=6, offset=3.0)
    assert len(tokens) == 2, "expected 2 ori tokens for a 12-bp duplex (window=6)"
    covered = set()
    for t in tokens:
        covered.update(range(t["bp_start"], t["bp_end"] + 1))
    assert covered == set(range(1, len(pairs) + 1)), "ori tokens must cover every base pair"


def test_hbond_candidates_nonempty_and_typed(duplex):
    hb = hbond_candidates(duplex)
    assert len(hb) > 0
    roles = {h["role"] for h in hb}
    assert roles <= {"donor", "acceptor"}
