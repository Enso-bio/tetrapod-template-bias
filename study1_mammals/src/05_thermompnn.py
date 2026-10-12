"""Per-site ThermoMPNN energies with the reconstructed placental ancestor threaded onto the human AlphaFold backbone.

Also stores, for control analyses: energies with the native human sequence, pLDDT, relative solvent accessibility and
a backbone-dihedral secondary-structure class for every site.
"""
import math

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PPBuilder, ShrakeRupley
from Bio.SeqUtils import seq1

from common import DATA, RESULTS
from thermo import load_model, parse, site_energies

IQ_AA = 'ARNDCQEGHILKMFPSTWYV'
AA20 = 'ACDEFGHIKLMNPQRSTVWY'
TO_AA20 = [IQ_AA.index(a) for a in AA20]
RSA_MAX = {'A': 129, 'R': 274, 'N': 195, 'D': 193, 'C': 167, 'Q': 225, 'E': 223, 'G': 104, 'H': 224, 'I': 197, 'L': 201,
           'K': 236, 'M': 224, 'F': 240, 'P': 159, 'S': 155, 'T': 172, 'W': 285, 'Y': 263, 'V': 174}  # Tien et al. 2013
OUT = DATA / 'energies'
OUT.mkdir(exist_ok=True)


def structure_annotations(pdb_path, length):
    s = PDBParser(QUIET=True).get_structure('m', pdb_path)
    residues = [r for r in s[0]['A'] if 'CA' in r]
    assert len(residues) == length
    ShrakeRupley(n_points=100, probe_radius=1.4).compute(s[0], level='R')
    rsa = np.array([r.sasa / RSA_MAX.get(seq1(r.resname), 167) for r in residues])
    plddt = np.array([r['CA'].bfactor for r in residues])
    ss = np.array(['coil'] * length, dtype=object)
    index = {r.get_id()[1]: k for k, r in enumerate(residues)}
    for pp in PPBuilder().build_peptides(s[0]['A']):
        for res, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
            if phi is None or psi is None:
                continue
            phi, psi = math.degrees(phi), math.degrees(psi)
            k = index[res.get_id()[1]]
            if -160 <= phi <= -20 and -120 <= psi <= 50:
                ss[k] = 'helix_region'
            elif -180 <= phi <= -40 and (psi >= 90 or psi <= -150):
                ss[k] = 'strand_region'
    return plddt, rsa, ss.astype(str)


def main():
    asr = pd.read_csv(RESULTS / 'asr_summary.csv')
    asr = asr[asr.status == 'ok']
    model = load_model()
    rows = []
    for r in asr.itertuples():
        path = OUT / f'{r.symbol}.npz'
        if path.exists():
            continue
        probs = np.load(DATA / 'asr' / r.symbol / 'ancestor.npz')['probs'][:, TO_AA20]    # (L, 20) in AA20 order
        anc = ''.join(AA20[i] for i in probs.argmax(axis=1))
        pdb = DATA / 'alphafold' / f'{r.accession}.pdb'
        d = parse(pdb)
        human = d['seq_chain_A']
        if len(human) != len(anc):
            rows.append({'symbol': r.symbol, 'status': 'length_mismatch'})
            continue
        e_anc = site_energies(model, d, anc)
        e_hum = site_energies(model, d)
        plddt, rsa, ss = structure_annotations(pdb, len(human))
        np.savez_compressed(path, e_anc=e_anc, e_human=e_hum, probs=probs, ancestor=anc, human=human,
                            plddt=plddt, rsa=rsa, ss=ss)
    for r in asr.itertuples():
        path = OUT / f'{r.symbol}.npz'
        if path.exists():
            z = np.load(path)
            anc, human = str(z['ancestor']), str(z['human'])
            rows.append({'symbol': r.symbol, 'status': 'ok', 'length': len(human),
                         'ancestor_identity_to_human': float(np.mean([a == b for a, b in zip(anc, human)]))})
    summary = pd.DataFrame(rows)
    summary.to_csv(RESULTS / 'thermompnn_summary.csv', index=False)
    print(summary.status.value_counts().to_dict())


if __name__ == '__main__':
    main()
