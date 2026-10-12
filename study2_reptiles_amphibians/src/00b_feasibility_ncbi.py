"""Feasibility: NCBI-annotated (RefSeq) reptile and amphibian reference genomes with AnAge records."""
import sys
import time

import pandas as pd
import requests

ROOT, ANAGE = sys.argv[1], sys.argv[2]
UA = {'User-Agent': 'ComparativeLongevityStudy/1.0 (academic)'}
GROUPS = {'Testudines': 8459, 'Crocodylia': 1294634, 'Lepidosauria': 8504, 'Amphibia': 8292}


def reports(taxon):
    out, token = [], None
    while True:
        params = {'filters.has_annotation': 'true', 'filters.reference_only': 'true', 'page_size': 500}
        if token:
            params['page_token'] = token
        r = requests.get(f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxon}/dataset_report',
                         params=params, headers=UA, timeout=300)
        r.raise_for_status()
        js = r.json()
        out += js.get('reports', [])
        token = js.get('next_page_token')
        if not token:
            return out
        time.sleep(0.4)


rows = []
for group, tax in GROUPS.items():
    for rep in reports(tax):
        org = rep.get('organism', {})
        ann = rep.get('annotation_info', {})
        rows.append({'group': group, 'organism': org.get('organism_name'), 'taxon_id': org.get('tax_id'),
                     'accession': rep.get('accession'), 'annotation': ann.get('provider'),
                     'busco': (ann.get('busco') or {}).get('complete')})
    time.sleep(0.4)
ncbi = pd.DataFrame(rows)
ncbi['binomial'] = ncbi.organism.str.split().str[:2].str.join('_').str.lower()
ncbi.to_csv(f'{ROOT}/data/raw/ncbi_annotated_genomes.csv', index=False)
an = pd.read_csv(ANAGE)
m = ncbi.merge(an[['binomial', 'Class', 'Order', 'Maximum longevity (yrs)', 'Adult weight (g)', 'Data quality', 'Sample size']],
               on='binomial', how='left')
ok_life = m['Maximum longevity (yrs)'].notna() & m['Data quality'].str.lower().isin(['acceptable', 'high']) \
    & m['Sample size'].str.lower().isin(['small', 'medium', 'large', 'huge'])
m['eligible_with_mass'] = ok_life & (m['Adult weight (g)'] > 0)
m['eligible_no_mass'] = ok_life
pd.set_option('display.width', 250)
pd.set_option('display.max_rows', 300)
print(m[['group', 'organism', 'annotation', 'busco', 'Maximum longevity (yrs)', 'Adult weight (g)', 'Data quality', 'Sample size']]
      .sort_values(['group', 'organism']).to_string(index=False))
print('\nannotated reference genomes per group:', ncbi.groupby('group').size().to_dict())
print('eligible WITH mass per group:', m[m.eligible_with_mass].groupby('group').size().to_dict())
print('eligible lifespan (mass missing allowed):', m[m.eligible_no_mass].groupby('group').size().to_dict())
m.to_csv(f'{ROOT}/data/raw/ncbi_anage_feasibility.csv', index=False)
