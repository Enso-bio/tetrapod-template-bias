"""Bulk Ensembl release-116 files: human Compara homologies (+ partner-species files when needed) and peptide FASTA.

Replaces per-gene REST queries (plan amendment E3): same Compara one-to-one orthology calls, bulk download.
"""
import gzip
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

from common import DATA, RAW, USER_AGENT, _record

FTP = 'https://ftp.ensembl.org/pub/release-116'
OUT = RAW / 'compara'
PEP = RAW / 'pep'
OUT.mkdir(exist_ok=True)
PEP.mkdir(exist_ok=True)


def curl(url, path):
    if path.exists():
        return path
    tmp = path.with_suffix(path.suffix + '.part')
    for _ in range(5):
        p = subprocess.run(['curl.exe', '-s', '-S', '-f', '-L', '-A', USER_AGENT, '-o', str(tmp), url],
                           capture_output=True, text=True)
        if p.returncode == 0:
            tmp.replace(path)
            _record(url, path)
            return path
    raise RuntimeError(f'download failed {url}: {p.stderr}')


def partners_in_human_file(species):
    seen = set()
    with gzip.open(OUT / 'homo_sapiens.homologies.tsv.gz', 'rt') as f:
        f.readline()
        for line in f:
            x = line.split('\t', 8)
            if x[4] == 'ortholog_one2one':
                seen.add(x[7])
    return seen & species


def pep_url(species):
    listing = requests.get(f'{FTP}/fasta/{species}/pep/', headers={'User-Agent': USER_AGENT}, timeout=120).text
    names = sorted(set(re.findall(r'href="([^"]+\.pep\.all\.fa\.gz)"', listing)))
    if len(names) != 1:
        raise RuntimeError(f'{species}: pep files {names}')
    return f'{FTP}/fasta/{species}/pep/{names[0]}'


def main():
    sp = pd.read_csv(DATA / 'species_candidates.csv')
    species = set(sp.ensembl_name) - {'homo_sapiens'}
    missing = sorted(species - partners_in_human_file(species))
    print('species absent from the human homology file:', missing, flush=True)
    jobs = [(f'{FTP}/tsv/ensembl-compara/homologies/{s}/Compara.116.protein_default.homologies.tsv.gz',
             OUT / f'{s}.homologies.tsv.gz') for s in missing]
    for s in sorted(species | {'homo_sapiens'}):
        jobs.append((pep_url(s), PEP / f'{s}.pep.all.fa.gz'))
    with ThreadPoolExecutor(max_workers=10) as pool:
        for path in pool.map(lambda j: curl(*j), jobs):
            print('done', path.name, flush=True)


if __name__ == '__main__':
    main()
