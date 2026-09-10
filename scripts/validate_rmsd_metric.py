"""THE TEST: does dna_aligned_ca_rmsd return a sane number for a real cocrystal
vs its own untemplated rf3 refold?"""
import glob, json, os, sys, numpy as np
sys.path.insert(0,'scripts')
import biotite.structure.io.pdbx as pdbx, biotite.structure as struc
from filter_binder_block import (dna_aligned_ca_rmsd, protein_only_ca_rmsd,
                                 count_protein_dna_hbonds, _match_dna_chains,
                                 _dna_strand_sequences)

CT='/private/tmp/claude-502/-Users-campbell-mcduling-WMG-repos-rfd3-dna-binder-reproduction/c08b913d-7708-4f0e-814c-429673d71b60/scratchpad/ctrlpdb'
RAW='/private/tmp/claude-502/-Users-campbell-mcduling-WMG-repos-rfd3-dna-binder-reproduction/c08b913d-7708-4f0e-814c-429673d71b60/scratchpad/metric_control/raw'
man=json.load(open('/private/tmp/claude-502/-Users-campbell-mcduling-WMG-repos-rfd3-dna-binder-reproduction/c08b913d-7708-4f0e-814c-429673d71b60/scratchpad/metric_control/manifest.json'))

def crystal(pid, prot_chain, dna_chains):
    a=pdbx.get_structure(pdbx.CIFFile.read(f'{CT}/{pid}.cif'), model=1)
    a=a[~struc.filter_solvent(a)]
    keep=np.isin(a.chain_id,[prot_chain]+list(dna_chains)) & \
         (struc.filter_amino_acids(a)|struc.filter_nucleotides(a))
    x=a[keep]
    # Renumber each DNA strand 1..N in res_id order. The metric matches DNA atoms on
    # (chain, res_id, atom_name), and crystal DNA is numbered 201-221 / 51-61 / 101-108
    # while an rf3 refold is always 1..N -- so without this only strands that happen to
    # already be 1..N match, and the fit silently uses a subset.
    nuc=struc.filter_nucleotides(x)
    new=x.res_id.copy()
    for c in set(x.chain_id[nuc]):
        m=nuc&(x.chain_id==c)
        for i,r in enumerate(sorted(set(x.res_id[m].tolist())), start=1):
            new[m&(x.res_id==r)]=i
    x.res_id=new
    return x

print(f"{'pdb':6} {'class':12} {'DNA-aligned':>12} {'protein-only':>13} {'maj-groove H-b':>15}  {'gate':>6}")
print('-'*72)
rows=[]
for rec in man:
    cif=glob.glob(f"{RAW}/{rec['fold_id']}/**/*_model.cif", recursive=True)
    if not cif:
        print(f"{rec['pdb']:6} no refold CIF"); continue
    ref=pdbx.get_structure(pdbx.CIFFile.read(cif[0]), model=1)
    ref=ref[~struc.filter_solvent(ref)]
    ref=ref[struc.filter_amino_acids(ref)|struc.filter_nucleotides(ref)]
    des=crystal(rec['pdb'], rec['prot_chain'], rec['dna_chains'])
    try:
        v,n=dna_aligned_ca_rmsd(des, ref)
        pv=protein_only_ca_rmsd(des, ref)
        _,mg=count_protein_dna_hbonds(ref)
        mapping=_match_dna_chains(des, ref)
        dn=_dna_strand_sequences(des); rn=_dna_strand_sequences(ref)
        assert len(mapping)==len(rn), 'not all refold strands paired'
        tag='PALINDROME' if rec['palindromic'] else ''
        print(f"{rec['pdb']:6} {tag:12} {v:12.2f} {pv:13.2f} {mg:15d}  {'PASS' if v<8 else 'FAIL':>6}")
        rows.append(dict(pdb=rec['pdb'], dna_aligned=round(v,2), protein_only=round(pv,2),
                         major_groove_hbonds=int(mg), palindromic=rec['palindromic'],
                         chain_map={str(k):str(x) for k,x in mapping.items()}))
    except Exception as e:
        print(f"{rec['pdb']:6} RAISED {type(e).__name__}: {str(e)[:70]}")
        rows.append(dict(pdb=rec['pdb'], error=f'{type(e).__name__}: {e}'))
json.dump(rows, open('/tmp/metric_control_results.json','w'), indent=2)
ok=[r for r in rows if 'dna_aligned' in r and not r['palindromic']]
if ok:
    print(f"\nnon-palindromic: {sum(1 for r in ok if r['dna_aligned']<8)}/{len(ok)} under the 8 A gate; "
          f"median {np.median([r['dna_aligned'] for r in ok]):.2f} A")
