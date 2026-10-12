"""Count one-to-one orthologies per partner species in the human genome-specific Compara file."""
import gzip
from collections import Counter

import pandas as pd

from common import DATA, RAW

sp = set(pd.read_csv(DATA / 'species_candidates.csv').ensembl_name)
c = Counter()
first = Counter()
with gzip.open(RAW / 'compara' / 'homo_sapiens.homologies.tsv.gz', 'rt') as f:
    header = f.readline().rstrip('\n').split('\t')
    print(header)
    for line in f:
        x = line.rstrip('\n').split('\t')
        if x[4] != 'ortholog_one2one':
            continue
        a, b = x[2], x[7]
        other = b if a == 'homo_sapiens' else a
        if other in sp:
            c[other] += 1
            first[a == 'homo_sapiens'] += 1
print('human listed first:', dict(first))
print(sorted(c.items(), key=lambda kv: kv[1]))
