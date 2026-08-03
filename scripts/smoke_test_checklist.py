#!/usr/bin/env python3
"""Track smoke-test execution progress.

Maintains a checklist of smoke-test stages and gates, saving progress to a JSON file.
Use to track which stages have been completed and what decisions were made.

Usage:
    python smoke_test_checklist.py --init                    # Create fresh checklist
    python smoke_test_checklist.py --stage binder_0 --done   # Mark stage complete
    python smoke_test_checklist.py --show                    # Print current status
    python smoke_test_checklist.py --decision h_bond_subset "N7,O6 on purine core"
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from datetime import datetime

CHECKLIST_FILE = Path("results/smoke_test_checklist.json")

DEFAULT_CHECKLIST = {
    "session_start": None,
    "preflight": {
        "offtarget_panel_validated": False,
        "hbond_subsetting_decided": False,
        "pecli_auth_verified": False,
        "gpu_budget_estimated": False,
        "sampler_config_smoke_test": False,
    },
    "binder_block": {
        "stage_0_fold_duplex": {"done": False, "timestamp": None, "notes": ""},
        "stage_1_conditioning": {"done": False, "timestamp": None, "notes": ""},
        "stage_2_rfd3na_specs": {"done": False, "timestamp": None, "notes": ""},
        "stage_3_diffuse": {"done": False, "timestamp": None, "notes": "n_designs: ?"},
        "stage_4_relax": {"done": False, "timestamp": None, "notes": ""},
        "stage_5_ligandmpnn": {"done": False, "timestamp": None, "notes": ""},
        "stage_6_refold": {"done": False, "timestamp": None, "notes": "oracles: ?"},
        "stage_7_filter": {"done": False, "timestamp": None, "notes": "pass_rate: ?"},
    },
    "specificity_block": {
        "stage_1_resample": {"done": False, "timestamp": None, "notes": ""},
        "stage_2_ontarget_fold": {"done": False, "timestamp": None, "notes": ""},
        "stage_3_allbyall_inputs": {"done": False, "timestamp": None, "notes": ""},
        "stage_4_allbyall_fold": {"done": False, "timestamp": None, "notes": ""},
        "stage_5_delta_minpae": {"done": False, "timestamp": None, "notes": ""},
    },
    "post_smoke_test": {
        "metrics_match_paper_trends": False,
        "no_crashes": False,
        "file_handoffs_work": False,
        "cost_estimated": False,
        "decisions_reviewed": False,
    },
    "decisions": {
        "hbond_subsetting": None,
        "oracle_selection": None,
        "filter_thresholds": None,
    },
    "notes": "",
}


def init_checklist():
    """Create a fresh checklist."""
    CHECKLIST_FILE.parent.mkdir(parents=True, exist_ok=True)
    checklist = DEFAULT_CHECKLIST.copy()
    checklist["session_start"] = datetime.now().isoformat()
    with open(CHECKLIST_FILE, "w") as f:
        json.dump(checklist, f, indent=2)
    print(f"✓ Created fresh checklist at {CHECKLIST_FILE}")


def load_checklist():
    """Load existing checklist."""
    if not CHECKLIST_FILE.exists():
        print(f"No checklist found at {CHECKLIST_FILE}")
        print("Run with --init to create one")
        return None
    with open(CHECKLIST_FILE) as f:
        return json.load(f)


def mark_done(stage: str):
    """Mark a stage as complete."""
    checklist = load_checklist()
    if checklist is None:
        return

    # Try to find the stage in any section
    found = False
    for section in ["binder_block", "specificity_block", "preflight", "post_smoke_test"]:
        if stage in checklist.get(section, {}):
            if isinstance(checklist[section][stage], dict):
                checklist[section][stage]["done"] = True
                checklist[section][stage]["timestamp"] = datetime.now().isoformat()
            else:
                checklist[section][stage] = True
            found = True
            break

    if not found:
        print(f"✗ Stage '{stage}' not found in checklist")
        return

    with open(CHECKLIST_FILE, "w") as f:
        json.dump(checklist, f, indent=2)
    print(f"✓ Marked {stage} as complete")


def set_decision(key: str, value: str):
    """Record a decision."""
    checklist = load_checklist()
    if checklist is None:
        return

    if key in checklist.get("decisions", {}):
        checklist["decisions"][key] = value
    else:
        print(f"⚠ Decision key '{key}' not recognized. Adding as custom decision.")
        if "decisions" not in checklist:
            checklist["decisions"] = {}
        checklist["decisions"][key] = value

    with open(CHECKLIST_FILE, "w") as f:
        json.dump(checklist, f, indent=2)
    print(f"✓ Recorded decision: {key} = {value}")


def show_status():
    """Print current checklist status."""
    checklist = load_checklist()
    if checklist is None:
        return

    print("\n" + "=" * 70)
    print(f"SMOKE TEST PROGRESS — Session: {checklist.get('session_start', 'unknown')}")
    print("=" * 70)

    # Preflight
    pf = checklist.get("preflight", {})
    print(f"\nPREFLIGHT:")
    for k, v in pf.items():
        status = "✓" if v else "○"
        print(f"  {status} {k.replace('_', ' ')}")

    # Binder block
    print(f"\nBINDER BLOCK:")
    bb = checklist.get("binder_block", {})
    for stage, info in bb.items():
        if isinstance(info, dict):
            status = "✓" if info.get("done") else "○"
            notes = f" — {info.get('notes', '')}" if info.get('notes') else ""
            print(f"  {status} {stage.replace('_', ' ')}{notes}")
        else:
            print(f"  ○ {stage}")

    # Specificity block
    print(f"\nSPECIFICITY BLOCK:")
    sb = checklist.get("specificity_block", {})
    for stage, info in sb.items():
        if isinstance(info, dict):
            status = "✓" if info.get("done") else "○"
            notes = f" — {info.get('notes', '')}" if info.get('notes') else ""
            print(f"  {status} {stage.replace('_', ' ')}{notes}")
        else:
            print(f"  ○ {stage}")

    # Post smoke test
    print(f"\nPOST SMOKE TEST:")
    pst = checklist.get("post_smoke_test", {})
    for k, v in pst.items():
        status = "✓" if v else "○"
        print(f"  {status} {k.replace('_', ' ')}")

    # Decisions
    print(f"\nDECISIONS:")
    decisions = checklist.get("decisions", {})
    for k, v in decisions.items():
        if v:
            print(f"  • {k}: {v}")
        else:
            print(f"  ○ {k}: (pending)")

    print("\n" + "=" * 70)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--init", action="store_true", help="Create fresh checklist")
    ap.add_argument("--stage", help="Stage name to mark complete (use with --done)")
    ap.add_argument("--done", action="store_true", help="Mark stage as complete")
    ap.add_argument("--decision", nargs=2, metavar=("KEY", "VALUE"), help="Record a decision")
    ap.add_argument("--show", action="store_true", help="Print current status")
    args = ap.parse_args()

    if args.init:
        init_checklist()
    elif args.stage and args.done:
        mark_done(args.stage)
    elif args.decision:
        set_decision(args.decision[0], args.decision[1])
    elif args.show:
        show_status()
    else:
        show_status()


if __name__ == "__main__":
    main()
