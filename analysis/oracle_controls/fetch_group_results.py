#!/usr/bin/env python3
"""Download every run in a pecli run-group, bypassing the CLI's recency-limited resolver.

WHY NOT JUST `pecli results <id>`
---------------------------------
Because it cannot find these runs. `pecli runs --group oracle-controls-rf3` lists all 63
as SUCCEEDED, but `pecli results <that id>` answers `unknown run` for every one of them.

The cause is in pecli, not in our ids (verified in pecli/runner.py `_resolve`):

    find_by_short_id -> runs.list_for_user(user_id, limit=100)   # 100 NEWEST, per user
    fallback         -> runs.list_all(limit=200)                 # 200 NEWEST, team-wide

Both are recency windows. `pecli runs --group` instead queries the sparse group index,
which has no such window. So any run older than the ~100 most recent is listed but
unresolvable, and the two commands disagree about whether it exists. The oracle-control
folds are from 2026-07-30 and are long past both windows.

(This also retires the earlier theory that the run_ids committed in
folds/rf3_manifest.json were corrupted by the manifest race. They are fine. Neither the
committed ids, nor ids scraped from `pecli runs`, nor the full `YYYYMMDD-HHMMSS-xxxxxx`
ids read out of ~/.pecli/jobs.json resolve -- because none of them is the problem.)

The group query returns full Run records, each already carrying `s3_output`, so we skip
resolution entirely and call the same S3 download `fetch_results` would have called.
Free: S3 GETs only, no compute.

Worth filing upstream: `pecli results` should fall back to the group/id index that
`pecli runs` already uses, rather than reporting a listed run as unknown.

Usage:
    python fetch_group_results.py --group oracle-controls-rf3 --out-dir <dir>
    python fetch_group_results.py --group oracle-controls-rf3 --out-dir <dir> \
        --pattern '*confidences.json'      # just the PAE, much less to pull
"""
from __future__ import annotations
import argparse
import os
import re
import sys

DESC_RE = re.compile(r"oracle-controls(?:\s+TEMPLATED)?\s+(\S+)\s+\(")


def fold_id_from_description(desc, fallback):
    """Our submitters embed the fold_id in the run description; fall back to short_id."""
    m = DESC_RE.search(desc or "")
    return m.group(1) if m else fallback


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", required=True)
    ap.add_argument("--out-dir", required=True,
                    help="one subdirectory per fold_id (read from the run description)")
    ap.add_argument("--pattern", default="",
                    help="fnmatch glob restricting what is pulled, e.g. "
                         "'*confidences.json'. Empty fetches everything.")
    ap.add_argument("--limit", type=int, default=500)
    args = ap.parse_args()

    from pecli.config import load
    from pecli.aws import runs as pruns, s3
    from pecli.runner import slugify_group

    cfg = load().resolved()
    # pecli slugifies a --group label at SUBMIT time (runner.slugify_group), so the
    # stored value is e.g. "p25-fixed-cfgon" for a label typed "p25-fixed_cfgon". The
    # CLI slugifies the query too, so `pecli runs --group` matches either way; querying
    # the index directly does not, and returns a confident, wrong 0 runs.
    group = slugify_group(args.group)
    if group != args.group:
        print(f"group {args.group!r} slugified to {group!r} (pecli normalises at submit)")
    try:
        found = pruns.list_by_group(group, limit=args.limit, cfg=cfg)
    except Exception:
        # the sparse group index may not be deployed; the scan fallback is equivalent
        found = pruns.list_all_by_group(group, limit=args.limit, cfg=cfg)
    print(f"group {group}: {len(found)} run(s)")
    if not found:
        return 1

    os.makedirs(args.out_dir, exist_ok=True)
    ok, skipped, failed = 0, [], []
    for run in found:
        fold = fold_id_from_description(run.description, run.short_id)
        dest = os.path.join(args.out_dir, fold)
        if os.path.isdir(dest) and os.listdir(dest):
            ok += 1
            continue
        if run.status != "SUCCEEDED":
            skipped.append((fold, run.status))
            continue
        if not run.s3_output:
            skipped.append((fold, "no s3_output recorded"))
            continue
        try:
            paths = s3.download_outputs(cfg, run.s3_output, dest, pattern=args.pattern)
        except Exception as e:
            failed.append((fold, str(e)[:160]))
            continue
        if not paths:
            failed.append((fold, f"no objects matched pattern {args.pattern!r}"))
            continue
        ok += 1
        print(f"  {fold}: {len(paths)} file(s)")

    print(f"\ndownloaded {ok}/{len(found)} run(s) -> {args.out_dir}")
    for fold, why in skipped:
        print(f"  - skipped {fold}: {why}")
    for fold, why in failed:
        print(f"  ! failed {fold}: {why}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
