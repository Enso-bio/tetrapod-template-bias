"""Print the facts needed for the Supplementary Methods (versions, sources, dates), read from both projects."""
import glob
import json
from collections import Counter
from pathlib import Path

import pandas as pd

P1 = Path(__file__).resolve().parent.parent / 'study1_mammals'
P2 = Path(__file__).resolve().parent.parent / 'study2_reptiles_amphibians'


def provenance(root):
    rows = [json.loads(l) for l in open(root / 'data' / 'provenance.jsonl', encoding='utf-8')]
    df = pd.DataFrame(rows)
    def source(u):
        for key in ['genomics.senescence.info/species', 'genomics.senescence.info/genes', 'MamPhy', 'rest.ensembl.org',
                    'ftp.ensembl.org', 'rest.uniprot.org', 'alphafold.ebi.ac.uk/api', 'alphafold.ebi.ac.uk/files',
                    'eutils.ncbi', 'api.ncbi.nlm.nih.gov/datasets', 'ftp.ncbi.nlm.nih.gov/gene', 'timetree',
                    'figshare', 'github.com']:
            if key in u:
                return key
        return u[:60]
    df['source'] = df.url.map(source)
    return df.groupby('source').agg(n=('url', 'size'), first=('retrieved_utc', 'min'), last=('retrieved_utc', 'max'))


def af_versions(root, accessions):
    c = Counter()
    for a in accessions:
        f = root / 'data' / 'raw' / 'alphafold_api' / f'{a}.json'
        if f.exists():
            info = json.loads(f.read_text())
            c[max(m.get('latestVersion', 0) for m in info)] += 1
    return dict(c)


print('== Project 1 provenance'); print(provenance(P1).to_string())
print('== Project 2 provenance'); print(provenance(P2).to_string())
print('runtime P1:', json.load(open(P1 / 'runtime_versions.json')))
print('ensembl release P1:', json.load(open(P1 / 'data' / 'ensembl_release.json')).get('releases'))
g1 = pd.read_csv(P1 / 'results' / 'asr_summary.csv')
print('AF versions P1 (human templates):', af_versions(P1, g1.accession))
g2 = pd.read_csv(P2 / 'results' / 'asr_summary.csv')
print('AF versions P2 (zebrafish templates):', af_versions(P2, g2.zf_uniprot))
for name, folder in [('opossum', 'alphafold_monodelphis_domestica'), ('platypus', 'alphafold_ornithorhynchus_anatinus')]:
    accs = [Path(p).stem for p in glob.glob(str(P1 / 'data' / folder / '*.pdb'))]
    print(f'AF versions P1 {name}:', af_versions(P1, accs))
sp1 = pd.read_csv(P1 / 'data' / 'species_final.csv')
print('P1 eligible placental:', int((sp1.eligible & sp1.placental).sum()), '| outgroups:', int((sp1.eligible & ~sp1.placental).sum()))
print(sp1[sp1.eligible][['ensembl_name', 'assembly']].head(3).to_string())
sp2 = pd.read_csv(P2 / 'data' / 'species_final.csv')
print(sp2[sp2.eligible][['binomial', 'accession', 'annotation_provider', 'annotation_name']].head(3).to_string())
og = pd.read_csv(P2 / 'data' / 'outgroups.csv')
print(og[['organism', 'accession', 'annotation_name']].to_string())
print('gene_orthologs header file:', (P2 / 'data' / 'raw' / 'gene_orthologs.gz').stat().st_size)
print('timetree summary:', json.load(open(P2 / 'data' / 'timetree_upgma_summary.json')))
print('P2 masses:', sp2[sp2.eligible].mass_source.value_counts().to_dict())
