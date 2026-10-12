"""Independent check of the placental-ancestor node: re-root on the outgroup with Biopython and take the MRCA."""
import numpy as np
import pandas as pd
from Bio import Phylo

from common import DATA, RESULTS, SEED

sp = pd.read_csv(DATA / 'species_final.csv')
placental = set(sp[sp.placental].ensembl_name)
asr = pd.read_csv(RESULTS / 'asr_summary.csv')
rng = np.random.default_rng(SEED)
bad = 0
sample = rng.choice(asr.symbol, 60, replace=False)
for g in sample:
    t = Phylo.read(DATA / 'asr' / g / f'{g}.treefile', 'newick')
    tips = [x.name for x in t.get_terminals()]
    out = [x for x in tips if x not in placental]
    pl = [x for x in tips if x in placental]
    t.root_with_outgroup(out if len(out) > 1 else out[0])
    mrca = t.common_ancestor(pl)
    expected = asr.loc[asr.symbol == g, 'node'].iloc[0]
    if mrca.name != expected or {x.name for x in mrca.get_terminals()} != set(pl):
        bad += 1
        print('mismatch', g, mrca.name, expected)
print(f'checked {len(sample)} genes, mismatches: {bad}')
