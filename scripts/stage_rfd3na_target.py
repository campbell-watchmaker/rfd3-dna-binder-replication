#!/usr/bin/env python3
"""Stage the folded target duplex alongside a prepared rfd3na run.

An rfd3na design spec references its target structure by CONTAINER path, e.g.
`"input": "/workspace/prnp_duplex.cif"`. `pecli prepare` stages only config.json
and the spec itself, so without this step the run starts and immediately fails --
the CIF it points at was never uploaded.

The fix is pecli's `aux_files` companion mechanism: drop the CIF into the staging
directory and list it in the run manifest. `submit_prepared` uploads every
aux_files entry into the run's S3 input prefix, and the container's stage_in
mirrors them into /workspace/, which is exactly where the spec expects it.

Top-level placement is safe here: rfd3na's `find_input_spec()` globs only
`*.json` / `*.yaml` / `*.yml` for the design spec, so a sibling `.cif` cannot be
mistaken for a second spec. (Contrast rf3, whose spec discovery is broader --
`analysis/oracle_controls/stage_templated_run.py` puts its template CIF in a
`templates/` subdirectory for that reason.)

Usage:
    python stage_rfd3na_target.py <staging_dir> <duplex.cif>
"""
from __future__ import annotations
import json
import os
import shutil
import sys


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    sdir, cif = sys.argv[1], sys.argv[2]
    manifest_path = os.path.join(sdir, "manifest.json")
    if not os.path.isfile(manifest_path):
        print(f"no manifest.json in {sdir}")
        return 1
    manifest = json.load(open(manifest_path))
    if manifest.get("kind") == "pipeline":
        print("REFUSING: staged run is a pipeline, not a bare rfd3na run")
        return 1

    spec_name = manifest["input_filename"]
    spec = json.load(open(os.path.join(sdir, spec_name)))
    body = next(iter(spec.values()))
    want = body.get("input", "")
    base = os.path.basename(want)
    if base != os.path.basename(cif):
        print(f"REFUSING: spec points at {want!r} but you passed {cif!r} -- "
              "the container would not find it under that name")
        return 1

    shutil.copyfile(cif, os.path.join(sdir, base))
    aux = list(manifest.get("aux_files", []))
    if base not in aux:
        aux.append(base)
    manifest["aux_files"] = aux
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"staged {base} -> {sdir}")
    print(f"  spec input     : {want}")
    print(f"  aux_files      : {manifest['aux_files']}")
    print(f"  hbond acceptor : {body.get('select_hbond_acceptor')}")
    print(f"  hbond donor    : {body.get('select_hbond_donor')}")
    print(f"  ori_token      : {body.get('ori_token')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
