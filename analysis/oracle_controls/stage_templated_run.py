#!/usr/bin/env python3
"""Turn a BARE staged rf3 run into a TEMPLATED one, then re-run submit-time validation.

Why this two-step dance exists: pecli decides "bare fold vs paid msa->fold pipeline"
at PREPARE time, and its only signal is an in-input MSA -- fold_shapes.has_precomputed_msa
looks for `msa_path` or the `_pecli_rf3_msa_a3m` carrier on a COMPONENT. A templated
input has no `seq` protein component to hang either on (`msa_path` is a
SequenceComponent field; passing it to a CIF `path` component is a TypeError in
atomworks, and hanging it on a DNA component makes rf3 raise
"Unsupported chain type for MSAs: polydeoxyribonucleotide"). So a templated input
prepared directly ALWAYS auto-routes to the paid pipeline.

pecli's documented review/edit flow gets us there anyway:
  1. `pecli prepare rf3 --input <UNTEMPLATED carrier-bearing input>` -> a BARE run;
  2. replace the staged input with the TEMPLATED json (same filename);
  3. drop the template CIF at <sdir>/templates/<x>.cif and list it in the manifest's
     `aux_files` (submit_prepared uploads those into the matching S3 sub-prefix, and
     the container's stage_in mirrors sub-paths into /workspace/, where the rf3
     container's find_input_spec -- a TOP-LEVEL glob -- will not mistake it for a
     second fold spec);
  4. `pecli submit <id>`: it re-checks only the input extension + config.json, and
     does NOT re-run the auto-route decision, so the run stays bare and MSA-free.

This script does 2-3 and then re-runs 4's pre-upload validation locally (free).
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, "/Users/campbell.mcduling/WMG_repos/pecli")
from pecli import fold_shapes  # noqa: E402
from pecli.tools import get_tool  # noqa: E402


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__)
        print("usage: stage_templated.py <staging_dir> <templated_json> <template_cif>")
        return 2
    sdir, tjson, tcif = (Path(a) for a in sys.argv[1:])
    manifest_path = sdir / "run_manifest.json"
    if not manifest_path.is_file():
        # manifest name may differ across versions -- find it
        cands = [p for p in sdir.glob("*.json") if p.name.endswith("manifest.json")]
        if not cands:
            print(f"no manifest in {sdir}: {[p.name for p in sdir.iterdir()]}")
            return 1
        manifest_path = cands[0]
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("kind") == "pipeline":
        print("REFUSING: this staged run is an msa->fold PIPELINE, not a bare fold")
        return 1

    input_name = manifest["input_filename"]
    spec = json.loads(tjson.read_text())
    # the fold name inside the spec must match the staged filename's stem so outputs
    # land under the expected <name>_* files
    spec[0]["name"] = Path(input_name).stem
    (sdir / input_name).write_text(json.dumps(spec, indent=2) + "\n")

    aux_rel = f"templates/{tcif.name}"
    (sdir / "templates").mkdir(exist_ok=True)
    shutil.copyfile(tcif, sdir / aux_rel)
    aux = list(manifest.get("aux_files", []))
    if aux_rel not in aux:
        aux.append(aux_rel)
    manifest["aux_files"] = aux
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    # ---- re-run submit_prepared's pre-upload validation (no upload, no spend) ----
    tool = get_tool(manifest["tool"])
    tool.check_input(input_name)                       # raises on a bad extension
    cfg = tool.validate_config(json.loads((sdir / "config.json").read_text()))
    for rel in manifest["aux_files"]:
        assert (sdir / rel).is_file(), f"missing staged companion {rel}"

    data = json.loads((sdir / input_name).read_text())
    fmt = fold_shapes.detect_format(data)
    comps = data[0]["components"]
    paths = [c["path"] for c in comps if "path" in c]
    print(f"staged  : {sdir}")
    print(f"manifest: kind={manifest.get('kind', 'run')} tool={manifest['tool']} "
          f"aux_files={manifest['aux_files']}")
    print(f"input   : {input_name} detected_shape={fmt} "
          f"template_selection={data[0].get('template_selection')}")
    print(f"          components: {len(comps)}  path={paths}")
    print(f"          seq components: "
          f"{[(c.get('chain_id'), c.get('chain_type'), len(c['seq'])) for c in comps if 'seq' in c]}")
    print(f"          ccd components: "
          f"{[(c.get('chain_id'), c.get('ccd_code')) for c in comps if 'ccd_code' in c]}")
    print(f"config  : diffusion_batch_size={cfg['sampler']['diffusion_batch_size']} "
          f"seed={cfg['sampler']['seed']} n_recycles={cfg['sampler']['n_recycles']}")
    print("PASS: bare rf3 run, templated input + companion CIF staged, submit-time "
          "validation clean (nothing uploaded).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
