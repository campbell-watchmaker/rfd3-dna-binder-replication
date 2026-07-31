#!/usr/bin/env python3
"""Compute minPAE / ΔminPAE / spread / argmin-correctness for the oracle control panel.

Consumes a folds_manifest.json whose `pae_path` fields have been filled in after
the folds returned, and emits a per-fold CSV plus a per-protein summary.

METRICS
-------
minPAE(protein, dna)
    min over protein-token x DNA-token pairs of PAE(i, j), both orientations.
    The paper's specificity primitive.

ΔminPAE(protein)              [specific TFs only -- needs a designated on-target]
    min over off-targets of minPAE(off) − minPAE(on).
    Large positive => the protein is confidently paired only with its own site.

spread(protein)               [defined for every class]
    max(minPAE) − min(minPAE) across the panel. Unlike ΔminPAE this needs no
    on-target, so it is the metric that lets non-specific binders and
    non-binders participate: both should show spread ≈ 0, for opposite reasons.

argmin(protein)               [the sharpest test]
    which DNA target attains the lowest minPAE. For a specific TF this should be
    its own cognate site. Being a rank statistic it needs no cross-oracle
    calibration, so it compares oracles on equal footing.

READING THE RESULT
------------------
  class               ΔminPAE   absolute minPAE   argmin on own site
  specific            large +   low on-target     yes
  nonspecific_binder  n/a       low everywhere    n/a
  nonbinder           n/a       high everywhere   n/a

If specific TFs show spread ≈ 0, or their argmin is essentially random, then the
oracle's PAE does not carry usable sequence-specificity information and ΔminPAE
cannot rank designs on it -- which is exactly what this experiment exists to find
out before the metric is trusted on our own designs.

ORACLE OUTPUT FORMATS (they differ, and the difference matters)
--------------------------------------------------------------
rf3       <name>_confidences.json
          keys: pae [N,N], token_chain_ids [N], token_res_ids [N]
          Chain labels are present, so protein/DNA tokens are selected directly.
          NB labels carry an entity suffix ("A_1"), handled by _norm_chain().

protenix  <name>_full_data_sample_<rank>.json   (only when the run set
          output.need_atom_confidence = true; the default summary file has no
          per-token PAE at all)
          keys: token_pair_pae [N,N], token_asym_id [N], ...
          There are NO chain-id strings -- only integer asym ids. Protein vs DNA
          is therefore resolved positionally, from the order the entities were
          written into the input JSON, which the manifest records.

esmfold2  TWO files: <id>_pae.npy (float16 [N,N], Angstrom on a 0-32 scale --
          cast to float32) and <id>_pae_tokens.json (token_chain_ids + a `mapping`
          provenance string). Only exists for runs submitted with --emit-pae true.
          token_chain_ids is null when the container could not reconcile token
          boundaries (multiple ligand chains); see _load_esmfold2 for the guarded
          positional fallback used in that case.

Usage:
    python compute_control_metrics.py --manifest folds/rf3/folds_manifest.json \
        --out-csv results/oracle_controls/rf3_folds.csv \
        --out-summary results/oracle_controls/rf3_summary.csv
"""
from __future__ import annotations
import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
from compute_delta_minpae import _norm_chain, min_pae  # noqa: E402  (reuse, don't duplicate)


def _load_esmfold2(rec: dict):
    """esmfold2: a float16 .npy matrix + a separate token-map JSON.

    The matrix is float16 Angstrom on a 0-32 scale and must be cast to float32
    before any arithmetic. It is deliberately NOT symmetric (pae[i, j] is the
    error at j when aligned on i, the AlphaFold convention), which is fine here:
    min_pae() takes the minimum over BOTH pae[prot, dna] and pae[dna, prot], and a
    min over a block is direction-agnostic. Do not symmetrise it.

    Note also that PAE saturates near 32 A: a block at the ceiling means "no
    confident relative placement", not 31.8 A of measured error.
    """
    pae = np.load(rec["pae_path"]).astype(np.float32)
    n = pae.shape[0]

    tok_path = rec.get("pae_tokens_path")
    if not tok_path or not os.path.exists(tok_path):
        raise ValueError(
            f"{rec['fold_id']}: esmfold2 PAE matrix without its _pae_tokens.json "
            "sidecar -- the token->chain map is required to slice the matrix.")
    tok = json.load(open(tok_path))
    ids = tok.get("token_chain_ids")

    exp_p = rec["protein_len"] * rec["protein_copies"]
    exp_d = 2 * rec["dna_len"]
    n_lig = len(rec.get("ligands") or [])

    if ids is None:
        # FALLBACK. The guide's advice ("don't assume a positional mapping when the
        # map says unavailable") guards against an UNKNOWN ligand atom count -- a
        # ligand's token count is what the container cannot reconstruct. That does
        # not apply here: every ligand in this panel is ZN, which is MONOATOMIC, so
        # each Zn contributes exactly 1 token and the composition is fully
        # determinate (Zif268 = 90 protein + 24 + 24 DNA + 3x1 Zn = 141 tokens).
        # We therefore resolve positionally -- chains are tokenised one token per
        # residue in spec order, polymers first, then ligands -- but ONLY behind a
        # hard total-token assertion. If the total disagrees the composition is not
        # what we think it is, and we report no number rather than a plausible-
        # looking wrong one.
        expected = exp_p + exp_d + n_lig
        if n != expected:
            raise ValueError(
                f"{rec['fold_id']}: esmfold2 token map unavailable "
                f"({tok.get('mapping')!r}) and the positional fallback does not "
                f"reconcile: PAE is [{n},{n}] but the input implies {exp_p} protein + "
                f"{exp_d} DNA + {n_lig} monoatomic ion = {expected} tokens. "
                "Refusing to guess -- no metric for this fold.")
        prot = np.zeros(n, dtype=bool)
        dna = np.zeros(n, dtype=bool)
        prot[:exp_p] = True
        dna[exp_p:exp_p + exp_d] = True
        return pae, prot, dna, "positional_fallback"

    ids = np.asarray([_norm_chain(c) for c in ids])
    if ids.shape[0] != n:
        raise ValueError(
            f"{rec['fold_id']}: token_chain_ids has {ids.shape[0]} entries but the "
            f"PAE is [{n},{n}].")
    want_p = {_norm_chain(c) for c in rec["protein_chains"]}
    want_d = {_norm_chain(c) for c in rec["dna_chains"]}
    prot = np.isin(ids, list(want_p))
    dna = np.isin(ids, list(want_d))
    if not prot.any() or not dna.any():
        raise ValueError(
            f"{rec['fold_id']}: chain labels {sorted(set(ids.tolist()))} did not match "
            f"protein={sorted(want_p)} / dna={sorted(want_d)}")
    if int(prot.sum()) != exp_p or int(dna.sum()) != exp_d:
        raise ValueError(
            f"{rec['fold_id']}: token map gave protein={int(prot.sum())} (expected "
            f"{exp_p}) and dna={int(dna.sum())} (expected {exp_d}) -- the map and the "
            "submitted sequences disagree; do not trust the metric.")
    # "exact" = the container's own token_chain_ids were used, as opposed to the
    # guarded "positional_fallback" above. (The container spells an exact map
    # "positional: ..." in its `mapping` string; don't reuse that word here, it
    # would read as the fallback.)
    return pae, prot, dna, "exact"


def load_pae_and_masks(rec: dict):
    """Return (pae, prot_mask, dna_mask) for one fold, per-oracle.

    esmfold2 is handled by _load_esmfold2 (its PAE is a .npy + a sidecar, not one
    JSON); use that directly if you need the token-map provenance string.
    """
    if rec["oracle"] == "esmfold2":
        pae, prot, dna, _mapping = _load_esmfold2(rec)
        return pae, prot, dna

    path = rec["pae_path"]
    with open(path) as f:
        d = json.load(f)

    if rec["oracle"] == "rf3":
        pae = np.asarray(d["pae"], dtype=float)
        chains = np.asarray([_norm_chain(c) for c in d["token_chain_ids"]])
        want_p = {_norm_chain(c) for c in rec["protein_chains"]}
        want_d = {_norm_chain(c) for c in rec["dna_chains"]}
        prot = np.isin(chains, list(want_p))
        dna = np.isin(chains, list(want_d))
        if not prot.any() or not dna.any():
            raise ValueError(
                f"{rec['fold_id']}: chain labels {sorted(set(chains.tolist()))} did not "
                f"match protein={sorted(want_p)} / dna={sorted(want_d)}")
        return pae, prot, dna

    if rec["oracle"] == "openfold3":
        # openfold3 <name>_seed_<s>_sample_<k>_confidences.json
        #   keys: pae [N,N], pde [N,N], plddt [n_atoms]
        # There are NO token ids of any kind (not even integer asym ids), so protein
        # vs DNA is resolved purely positionally from the order the chains were
        # written into the query JSON: protein copies, then the 2 DNA strands, then
        # any ligands. The token-count assertion below is the only thing standing
        # between that assumption and a silently wrong number, so it is fatal.
        if "pae" not in d:
            raise ValueError(f"{rec['fold_id']}: no per-token 'pae' in {path}")
        pae = np.asarray(d["pae"], dtype=float)
        n = pae.shape[0]
        exp_p = rec["protein_len"] * rec["protein_copies"]
        exp_d = 2 * rec["dna_len"]
        exp_l = len(rec.get("ligands") or [])
        if n != exp_p + exp_d + exp_l:
            raise ValueError(
                f"{rec['fold_id']}: openfold3 PAE is [{n},{n}] but the input implies "
                f"{exp_p} protein + {exp_d} DNA + {exp_l} ligand = {exp_p + exp_d + exp_l} "
                "tokens. The positional mapping is wrong -- do not trust the metric.")
        prot = np.zeros(n, dtype=bool)
        dna = np.zeros(n, dtype=bool)
        prot[:exp_p] = True
        dna[exp_p:exp_p + exp_d] = True
        return pae, prot, dna

    # protenix: positional resolution via integer asym ids
    key = next((k for k in ("token_pair_pae", "pae", "predicted_aligned_error") if k in d), None)
    if key is None:
        raise ValueError(
            f"{rec['fold_id']}: no per-token PAE in {path}. If this is a protenix run, it "
            "was probably submitted without --need-atom-confidence true, in which case only "
            "the summary confidence exists and the full PAE was never written to disk.")
    pae = np.asarray(d[key], dtype=float)
    asym = np.asarray(d["token_asym_id"])
    order = sorted(set(asym.tolist()))  # ascending; base (0 vs 1) does not matter
    ncopies = rec["protein_copies"]
    if len(order) < ncopies + 2:
        raise ValueError(
            f"{rec['fold_id']}: expected >= {ncopies + 2} asym ids "
            f"(protein x{ncopies} + 2 DNA strands), found {order}")
    prot_ids = order[:ncopies]
    dna_ids = order[ncopies:ncopies + 2]
    prot = np.isin(asym, prot_ids)
    dna = np.isin(asym, dna_ids)
    # sanity: token counts must match the sequences we submitted
    exp_p = rec["protein_len"] * ncopies
    exp_d = 2 * rec["dna_len"]
    if int(prot.sum()) != exp_p or int(dna.sum()) != exp_d:
        raise ValueError(
            f"{rec['fold_id']}: positional asym mapping gave protein={int(prot.sum())} "
            f"(expected {exp_p}) and dna={int(dna.sum())} (expected {exp_d}). The entity "
            "order assumption is wrong for this run -- do not trust the metric.")
    return pae, prot, dna


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    manifest = json.load(open(args.manifest))
    done = [m for m in manifest if m.get("pae_path") not in (None, "", "FILL_AFTER_FOLD")
            and os.path.exists(m["pae_path"])]
    missing = len(manifest) - len(done)

    rows, errors = [], []
    for rec in done:
        try:
            if rec["oracle"] == "esmfold2":
                # keep the token-map provenance in the table: a row scored off the
                # positional fallback must be visibly flagged, not silently equal
                # to one scored off an exact map.
                pae, prot, dna, token_map = _load_esmfold2(rec)
            else:
                pae, prot, dna = load_pae_and_masks(rec)
                token_map = ""
            mp = min_pae(pae, prot, dna)
            rows.append({
                "fold_id": rec["fold_id"], "protein": rec["protein"], "klass": rec["klass"],
                "dna_id": rec["dna_id"], "is_on_target": rec["is_on_target"],
                "oracle": rec["oracle"], "min_pae": round(mp, 4),
                "n_protein_tokens": int(prot.sum()), "n_dna_tokens": int(dna.sum()),
                "iptm": rec.get("iptm"), "ptm": rec.get("ptm"), "plddt": rec.get("plddt"),
                "token_map": token_map,
            })
        except Exception as e:  # a bad fold must not silently vanish from the table
            errors.append(f"{rec['fold_id']}: {e}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)
    if rows:
        with open(args.out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    # ---- per-protein summary ----
    by_prot = {}
    for r in rows:
        by_prot.setdefault(r["protein"], []).append(r)

    on_target_of = {m["protein"]: m.get("on_target_for_protein") for m in manifest}
    summary = []
    for prot, rs in sorted(by_prot.items()):
        vals = {r["dna_id"]: r["min_pae"] for r in rs}
        klass = rs[0]["klass"]
        on = on_target_of.get(prot)
        argmin_dna = min(vals, key=vals.get)
        spread = max(vals.values()) - min(vals.values())
        delta = None
        on_val = None
        if on and on in vals and len(vals) > 1:
            on_val = vals[on]
            offs = [v for k, v in vals.items() if k != on]
            delta = min(offs) - on_val
        summary.append({
            "protein": prot, "klass": klass, "n_dna": len(vals),
            "on_target": on or "", "min_pae_on_target": on_val,
            "delta_min_pae": None if delta is None else round(delta, 4),
            "min_pae_min": round(min(vals.values()), 4),
            "min_pae_max": round(max(vals.values()), 4),
            "min_pae_mean": round(float(np.mean(list(vals.values()))), 4),
            "spread": round(spread, 4),
            "argmin_dna": argmin_dna,
            "argmin_is_on_target": (argmin_dna == on) if on else None,
        })

    with open(args.out_summary, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()) if summary else ["protein"])
        w.writeheader()
        w.writerows(summary)

    # ---- report ----
    print(f"folds scored: {len(rows)}   missing PAE: {missing}   errors: {len(errors)}")
    for e in errors:
        print("  ! " + e)
    print()
    hdr = (f"{'protein':11s} {'class':18s} {'onPAE':>7s} {'ΔminPAE':>8s} "
           f"{'min':>6s} {'max':>6s} {'spread':>7s} {'argmin':>12s} {'ok':>3s}")
    print(hdr)
    print("-" * len(hdr))
    for s in summary:
        d = "" if s["delta_min_pae"] is None else f"{s['delta_min_pae']:+8.3f}"
        o = "" if s["min_pae_on_target"] is None else f"{s['min_pae_on_target']:7.3f}"
        ok = {True: "yes", False: "NO", None: "-"}[s["argmin_is_on_target"]]
        print(f"{s['protein']:11s} {s['klass']:18s} {o:>7s} {d:>8s} "
              f"{s['min_pae_min']:6.3f} {s['min_pae_max']:6.3f} {s['spread']:7.3f} "
              f"{s['argmin_dna']:>12s} {ok:>3s}")

    specific = [s for s in summary if s["klass"] == "specific"]
    hits = [s for s in specific if s["argmin_is_on_target"]]
    if specific:
        print(f"\nargmin hit rate (specific TFs): {len(hits)}/{len(specific)}")
        print(f"  random-chance baseline: ~1/{summary[0]['n_dna']} per TF "
              f"= {len(specific) / summary[0]['n_dna']:.2f} expected hits")
    print(f"\nper-fold table: {args.out_csv}\nsummary: {args.out_summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
