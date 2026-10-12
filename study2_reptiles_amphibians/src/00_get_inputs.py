"""Inputs that 01_species_phylo.py expects in data/raw.

AnAge and GenAge are the same archives used in Study 1 (copied from ../study1_mammals/data/raw when that study has
been run; downloaded otherwise). AmphiBIO v1 and ReptTraits v1-2 are downloaded from figshare. Every file is logged in
data/provenance.jsonl.
"""
import shutil

from common import RAW, ROOT, curl, record

STUDY1_RAW = ROOT.parent / 'study1_mammals' / 'data' / 'raw'
FILES = {'anage.zip': 'https://genomics.senescence.info/species/dataset.zip',
         'genage_human.zip': 'https://genomics.senescence.info/genes/human_genes.zip',
         'AmphiBIO_v1.zip': 'https://ndownloader.figshare.com/files/8828578',
         'ReptTraits_v1-2.xlsx': 'https://ndownloader.figshare.com/files/45408133'}


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dst = RAW / name
        if dst.exists():
            continue
        if (STUDY1_RAW / name).exists():
            shutil.copy2(STUDY1_RAW / name, dst)
            record(f'{url} (copied from study1_mammals)', dst)
        else:
            curl(url, dst)
        print('ready', name)


if __name__ == '__main__':
    main()
