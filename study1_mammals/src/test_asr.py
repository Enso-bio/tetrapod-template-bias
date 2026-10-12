"""Smoke test of the ASR step on one gene: IQ-TREE run, placental node detection, state-file columns."""
import importlib.util
import io
import sys

import pandas as pd
from Bio import Phylo

from common import DATA, RAW

spec = importlib.util.spec_from_file_location('asr', DATA.parent / 'src' / '04_asr.py')
asr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(asr)

symbol = sys.argv[1] if len(sys.argv) > 1 else 'SIRT6'
sp = pd.read_csv(DATA / 'species_candidates.csv')
seqs = asr.read_fasta(DATA / 'alignments' / f'{symbol}.fasta')
names = [s for s in ['homo_sapiens'] + sorted(sp.ensembl_name) if s in seqs]
names = list(dict.fromkeys(names))
labels = {k: v for k, v in zip(sp.tree_label, sp.ensembl_name) if v in names}
t = asr.prune_tree(RAW / 'mcc_node.tre', labels)
nwk = asr.newick(t)
print(nwk[:300])
placentals = sorted(sp[sp.placental].ensembl_name)
print(asr.run_gene((symbol, names, placentals, nwk)))
st = (DATA / 'asr' / symbol / f'{symbol}.state').read_text().splitlines()
print([l for l in st if not l.startswith('#')][0][:200])
