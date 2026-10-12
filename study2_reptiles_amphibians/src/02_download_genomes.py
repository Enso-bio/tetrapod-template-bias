"""Protein FASTA + GFF3 for every candidate and outgroup genome (NCBI Datasets) and NCBI gene_orthologs."""
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from common import DATA, RAW, curl

DATASETS = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download' \
           '?include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_GFF'


def main():
    sp = pd.read_csv(DATA / 'species_candidates.csv')
    out = pd.read_csv(DATA / 'outgroups.csv')
    accs = sorted(set(sp[sp.candidate].accession) | set(out.accession))
    jobs = [(DATASETS.format(acc=a), RAW / 'genomes' / f'{a}.zip') for a in accs]
    jobs.append(('https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz', RAW / 'gene_orthologs.gz'))
    with ThreadPoolExecutor(max_workers=6) as pool:
        for p in pool.map(lambda j: curl(*j), jobs):
            print('done', p.name, round(p.stat().st_size / 1e6, 1), 'MB', flush=True)


if __name__ == '__main__':
    main()
