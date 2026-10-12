"""Independent check of the tetrapod-ancestor node: re-root on the fishes with Biopython and take the MRCA."""
import numpy as np
import pandas as pd
from Bio import Phylo

from common import DATA, RESULTS, SEED

FISH = {'danio_rerio', 'lepisosteus_oculatus', 'latimeria_chalumnae'}
asr = pd.read_csv(RESULTS / 'asr_summary.csv')
rng = np.random.default_rng(SEED)
bad = 0
sample = rng.choice(asr.human_gene, 60, replace=False)
for g in sample:
    t = Phylo.read(DATA / 'asr' / str(g) / f'{g}.treefile', 'newick')
    tips = [x.name for x in t.get_terminals()]
    out = [x for x in tips if x in FISH]
    ing = [x for x in tips if x not in FISH]
    t.root_with_outgroup(out if len(out) > 1 else out[0])
    mrca = t.common_ancestor(ing)
    if mrca.name != asr.loc[asr.human_gene == g, 'node'].iloc[0] or {x.name for x in mrca.get_terminals()} != set(ing):
        bad += 1
        print('mismatch', g)
print(f'checked {len(sample)} genes, mismatches: {bad}')
sp = pd.read_csv(DATA / 'species_final.csv')
print('excluded species:', sp[sp.candidate & ~sp.eligible][['binomial', 'group', 'gene_coverage_fraction']].to_string(index=False))
print('final by group:', sp[sp.eligible].group.value_counts().to_dict())
