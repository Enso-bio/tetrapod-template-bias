"""Human gene universe (UniProt Swiss-Prot), GenAge candidate set and the fixed random order for the background."""
import io
import zipfile
from urllib.parse import urlencode

import numpy as np
import pandas as pd

from common import DATA, RAW, SEED, fetch

OLD_FAMILIES = ['SIRT6', 'SOD1', 'SOD2', 'CAT', 'GPX1', 'PRDX1', 'HSPA8', 'OGG1', 'PARK7', 'ACTB', 'GAPDH', 'TPI1']


def main():
    q = {'query': '(organism_id:9606) AND (reviewed:true)', 'format': 'tsv',
         'fields': 'accession,gene_primary,length,sequence,protein_name'}
    uni = pd.read_csv(io.BytesIO(fetch('https://rest.uniprot.org/uniprotkb/stream?' + urlencode(q),
                                       RAW / 'uniprot_human_swissprot.tsv', timeout=600)), sep='\t')
    uni.columns = ['accession', 'symbol', 'length', 'sequence', 'protein_name']
    uni = uni.dropna(subset=['symbol'])
    uni = uni[~uni.symbol.duplicated(keep=False)]           # ambiguous symbols cannot be queried reliably
    uni = uni[(uni.length >= 50) & (uni.length <= 1000)].copy()

    z = zipfile.ZipFile(io.BytesIO(fetch('https://genomics.senescence.info/genes/human_genes.zip', RAW / 'genage_human.zip')))
    name = [n for n in z.namelist() if n.endswith('.csv')][0]
    genage = pd.read_csv(io.BytesIO(z.read(name)))
    genage.to_csv(DATA / 'genage_human.csv', index=False)
    gen_acc = set(genage['uniprot'].dropna().astype(str).str.strip())
    gen_sym = set(genage['symbol'].dropna().astype(str).str.strip())

    uni['genage'] = uni.accession.isin(gen_acc) | uni.symbol.isin(gen_sym)
    uni['old_family'] = uni.symbol.isin(OLD_FAMILIES)
    rng = np.random.default_rng(SEED)
    background = uni[~uni.genage].copy()
    background['background_rank'] = rng.permutation(len(background))
    uni = uni.merge(background[['accession', 'background_rank']], on='accession', how='left')
    uni.to_csv(DATA / 'gene_universe_candidates.csv', index=False)
    print('Swiss-Prot human 50-1000 aa, unique symbols:', len(uni), '| GenAge:', int(uni.genage.sum()),
          'of', len(genage), 'GenAge entries | old families present:', sorted(uni[uni.old_family].symbol))


if __name__ == '__main__':
    main()
