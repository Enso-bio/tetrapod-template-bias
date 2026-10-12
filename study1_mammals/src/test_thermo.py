"""Validate the vectorized ThermoMPNN wrapper against the original per-mutation code path."""
import sys
import time

import numpy as np

from common import DATA
from thermo import AA20, load_model, parse, reference_ddg, site_energies

pdbs = sorted((DATA / 'alphafold').glob('*.pdb'))[:3]
model = load_model()
for pdb in pdbs:
    d = parse(pdb)
    t = time.time()
    E = site_energies(model, d)
    dt = time.time() - t
    seq = d['seq_chain_A']
    L = len(seq)
    pos = list(np.linspace(0, L - 1, 6).astype(int))
    ref = reference_ddg(model, d, pos)
    mine = np.array([[E[i, AA20.index(b)] - E[i, AA20.index(seq[i])] for b in AA20] for i in pos])
    print(pdb.name, 'L', L, f'{dt:.2f}s', 'max abs diff vs original', float(np.abs(ref - mine).max()))
    # sign check: mean ddG of mutations to proline vs to alanine over all sites
    wt = np.array([E[i, AA20.index(seq[i])] for i in range(L) if seq[i] in AA20])
    Ev = np.array([E[i] for i in range(L) if seq[i] in AA20])
    ddg = Ev - wt[:, None]
    print('   mean ddG to P', round(float(ddg[:, AA20.index('P')].mean()), 3), '| to A',
          round(float(ddg[:, AA20.index('A')].mean()), 3), '| mean over all', round(float(ddg.mean()), 3))
    # threading a different sequence changes the energies
    shuffled = ''.join(np.random.default_rng(0).permutation(list(seq)))
    E2 = site_energies(model, d, shuffled)
    print('   corr energies native vs shuffled-sequence context', round(float(np.corrcoef(E.ravel(), E2.ravel())[0, 1]), 3))
