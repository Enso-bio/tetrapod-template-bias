"""Feasibility: Ensembl reptile and amphibian genomes with AnAge longevity/mass records."""
import json
import sys

import pandas as pd
import requests

ROOT = sys.argv[1]
ANAGE = sys.argv[2]
UA = {'User-Agent': 'ComparativeLongevityStudy/1.0 (academic)', 'Accept': 'application/json'}


def genomes(taxon):
    r = requests.get(f'https://rest.ensembl.org/info/genomes/taxonomy/{taxon}?content-type=application/json', headers=UA, timeout=300)
    r.raise_for_status()
    return r.json()


rows = []
for group in ['Sauropsida', 'Amphibia', 'Testudines', 'Lepidosauria', 'Crocodylia', 'Aves', 'Actinopterygii', 'Coelacanthiformes', 'Dipnoi']:
    try:
        for g in genomes(group):
            rows.append({'group': group, 'ensembl_name': g['name'], 'display': g['display_name'],
                         'taxon_id': g['taxonomy_id'], 'assembly': g['assembly_name'], 'compara': g.get('has_peptide_compara')})
    except requests.HTTPError as e:
        print(group, 'error', e)
ens = pd.DataFrame(rows).drop_duplicates(['group', 'ensembl_name'])
ens.to_csv(f'{ROOT}/data/raw/ensembl_genomes_by_group.csv', index=False)
print(ens.groupby('group').size().to_string())

an = pd.read_csv(ANAGE)
an['binomial'] = an['binomial'].str.lower()
non_avian = ens[ens.group.isin(['Testudines', 'Lepidosauria', 'Crocodylia', 'Amphibia'])].drop_duplicates('ensembl_name').copy()
non_avian['binomial'] = non_avian.ensembl_name.str.split('_').str[:2].str.join('_')
m = non_avian.merge(an[['binomial', 'Class', 'Order', 'Maximum longevity (yrs)', 'Adult weight (g)', 'Data quality', 'Sample size']],
                    on='binomial', how='left')
pd.set_option('display.width', 220)
pd.set_option('display.max_rows', 200)
print(m[['group', 'ensembl_name', 'compara', 'Order', 'Maximum longevity (yrs)', 'Adult weight (g)', 'Data quality', 'Sample size']]
      .sort_values(['group', 'ensembl_name']).to_string(index=False))
