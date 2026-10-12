"""Supplementary Methods tables (software, data sources, genome assemblies) generated from both projects' files."""
import json
import re
import zipfile
from pathlib import Path

import pandas as pd

P1 = Path(__file__).resolve().parent.parent / 'study1_mammals'
P2 = Path(__file__).resolve().parent.parent / 'study2_reptiles_amphibians'
OUT = Path(__file__).resolve().parent / 'paper' / 'generated'
# UniProtKB release reported by the REST API (X-UniProt-Release header, release date 2 Sep 2026); the data were
# retrieved on 8-9 Oct 2026, within that release.
UNIPROT_RELEASE = '2026\\_03'
THERMOMPNN_SHA256 = 'af449118ccfb4e34971d802321907ed82a8a035ae80f055bf8fc56009a4a838a'


def esc(s):
    return str(s).replace('\\', '').replace('_', '\\_').replace('&', '\\&').replace('%', '\\%')


def macro(name, rows):
    return f'\\newcommand{{\\{name}}}{{%\n' + '\n'.join(rows) + '\n}\n'


def dates(key, roots=(P1, P2)):
    """UTC retrieval date range of all logged downloads whose URL contains key."""
    d = sorted(json.loads(l)['retrieved_utc'][:10] for root in roots
               for l in open(root / 'data' / 'provenance.jsonl', encoding='utf-8') if key in json.loads(l)['url'])
    if not d:
        raise ValueError(f'no provenance entry for {key}')
    a, b = pd.Timestamp(d[0]), pd.Timestamp(d[-1])
    if a == b:
        return f'{a.day} {a:%b %Y}'
    return f'{a.day}--{b.day} {a:%b %Y}' if (a.year, a.month) == (b.year, b.month) else f'{a.day} {a:%b}--{b.day} {b:%b %Y}'


def hagr_build(zip_name):
    """Build number and release date stated in the release notes inside the HAGR archive."""
    with zipfile.ZipFile(P1 / 'data' / 'raw' / zip_name) as z:
        text = re.sub(r'<[^>]+>', ' ', z.read('release.html').decode('utf-8', 'replace'))
    build = re.search(r'Build (\d+)', text).group(1)
    date = re.search(r'Release date:\s*([A-Z][a-z]+ \d+, \d{4})', text).group(1)
    return f'Build {build} (released {pd.Timestamp(date):%d %b %Y})'.replace('(released 0', '(released ')


def main():
    rt = json.load(open(P1 / 'runtime_versions.json'))
    req2 = {l.split('==')[0].lower(): l.split('==')[1].strip()
            for l in open(P2 / 'requirements.txt', encoding='utf-8-sig') if '==' in l}
    software = [
        ('Python', rt['python'], 'All analyses'),
        ('NumPy \\citep{Harris2020}', rt['numpy'], 'Numerical computation, PGLS'),
        ('SciPy \\citep{Virtanen2020}', rt['scipy'], 'Optimisation of $\\lambda$, distributions'),
        ('pandas \\citep{McKinney2010}', rt['pandas'], 'Data handling'),
        ('statsmodels \\citep{Seabold2010}', rt['statsmodels'], 'OLS with HC3 errors (sensitivity)'),
        ('Biopython \\citep{Cock2009}', rt['biopython'], 'Pairwise alignment, structures, trees, accessibility'),
        ('PyTorch \\citep{Paszke2019}', f"{rt['torch'].split('+')[0]} (CUDA {rt['cuda']})", f"ThermoMPNN inference ({rt['gpu']})"),
        ('OmegaConf', rt['omegaconf'], 'ThermoMPNN configuration'),
        ('openpyxl', req2.get('openpyxl', '--'), 'Reading ReptTraits (Study 2)'),
        ('matplotlib \\citep{Hunter2007}', rt['matplotlib'], 'Figures'),
        ('ThermoMPNN \\citep{Dieckhaus2024}', 'Kuhlman-Lab main branch, downloaded 8 Oct 2026',
         'Stability prediction (weights thermoMPNN\\_default.pt)'),
        ('ProteinMPNN \\citep{Dauparas2022}', 'Weights v\\_48\\_020, bundled with ThermoMPNN', 'Frozen backbone and sequence embeddings'),
        ('IQ-TREE \\citep{Minh2020}', rt['iqtree'].split(' ')[0], 'Ancestral sequence reconstruction'),
        ('MAFFT \\citep{Katoh2013}', '7.526 (Windows package)', 'Multiple alignments (Study 2)'),
    ]
    rows = [f'{a} & {b} & {c} \\\\' for a, b, c in software]
    text = macro('SoftwareRows', rows)

    upham_commit = re.search(r'MamPhy_v1/([0-9a-f]{40})/', (P1 / 'src' / '01_phenotype_tree.py').read_text()).group(1)
    data = [
        ('AnAge (HAGR)', hagr_build('anage.zip'), dates('genomics.senescence.info/species'), '1, 2'),
        ('GenAge, human genes (HAGR)', hagr_build('genage_human.zip'), dates('genomics.senescence.info/genes'), '1, 2'),
        ('Ensembl', 'Release 116; REST and FTP (Compara homologies, peptides, genome list)', dates('ensembl.org'), '1'),
        ('UniProtKB', f'Release {UNIPROT_RELEASE}', dates('rest.uniprot.org'), '1, 2'),
        ('AlphaFold DB', 'Model version 6 (all templates used)', dates('alphafold.ebi.ac.uk'), '1, 2'),
        ('MamPhy (Upham et al.)', f'GitHub commit {upham_commit[:7]}; node-dated MCC tree and 100 posterior trees',
         dates('MamPhy'), '1'),
        ('NCBI Taxonomy (E-utilities)', 'Current at retrieval', dates('eutils.ncbi'), '1, 2'),
        ('NCBI Datasets', 'RefSeq genome annotations (Table~\\ref{tab:assembliesR})', dates('api.ncbi.nlm.nih.gov/datasets'), '2'),
        ('NCBI gene\\_orthologs', 'FTP file, Last-Modified 9 Oct 2026', dates('ftp.ncbi.nlm.nih.gov/gene'), '2'),
        ('TimeTree 5', 'Pairwise API (adjusted divergence times)', dates('timetree'), '2'),
        ('AmphiBIO', 'v1 (figshare file 8828578)', dates('files/8828578'), '2'),
        ('ReptTraits', 'v1-2 (figshare file 45408133)', dates('files/45408133'), '2'),
    ]
    rows = [f'{a} & {b} & {c} & {d} \\\\' for a, b, c, d in data]
    text += macro('DataRows', rows)

    def name(b):  # 'bison_bison' -> 'Bison bison'
        return esc(str(b).replace('_', ' ').capitalize())

    def header(label, ncol):
        return f'\\midrule\n\\multicolumn{{{ncol}}}{{@{{}}l}}{{\\textit{{{label}}}}} \\\\'

    # Species named as in AnAge (as in the species tables), with the genome actually used
    sp1 = pd.read_csv(P1 / 'data' / 'species_final.csv')
    sp1 = sp1[sp1.eligible].assign(species=lambda d: d.anage_binomial.fillna(d.binomial).map(name))
    rows = []
    for placental, label in [(True, None), (False, 'Outgroups (marsupials and monotremes)')]:
        if label:
            rows.append(header(label, 3))
        brk = lambda s: esc(s).replace('\\_', '\\_\\allowbreak{}')   # long identifiers may wrap at underscores
        rows += [f'\\emph{{{r.species}}} & \\texttt{{\\footnotesize {brk(r.ensembl_name)}}} & {brk(r.assembly)} \\\\'
                 for r in sp1[sp1.placental == placental].sort_values('species').itertuples()]
    text += macro('AssemblyRowsM', rows)

    def annot(r):  # 'GCF_x-RS_2024_02' -> 'RS_2024_02'; 'NCBI Genus species Annotation Release 102' -> 'Annotation Release 102'
        return esc(str(r.annotation_name).replace(f'{r.accession}-', '').replace(f'NCBI {r.organism} ', ''))
    sp2 = pd.read_csv(P2 / 'data' / 'species_final.csv')
    sp2 = sp2[sp2.eligible].assign(species=lambda d: d.anage_binomial.fillna(d.binomial).map(name)).sort_values('species')
    og = pd.read_csv(P2 / 'data' / 'outgroups.csv').assign(species=lambda d: d.organism.map(esc)).sort_values('species')
    rows = [f'\\emph{{{r.species}}} & \\emph{{{esc(r.organism)}}} & {esc(r.accession)} & {annot(r)} \\\\' for r in sp2.itertuples()]
    rows.append(header('Reference genomes (outgroups, templates and orthology reference)', 4))
    rows += [f'\\emph{{{r.species}}} & \\emph{{{esc(r.organism)}}} & {esc(r.accession)} & {annot(r)} \\\\' for r in og.itertuples()]
    text += macro('AssemblyRowsR', rows)
    sha = '\\allowbreak '.join(THERMOMPNN_SHA256[i:i + 16] for i in range(0, 64, 16))
    text += f'\\newcommand{{\\ThermoSHA}}{{\\texttt{{{sha}}}}}\n'
    text += f"\\newcommand{{\\RunPlatform}}{{{esc(rt['platform'].replace('-SP0', ''))}}}\n"
    (OUT / 'methods_tables.tex').write_text('% Generated by build_methods_tables.py. Do not edit by hand.\n' + text,
                                            encoding='utf-8')
    print('written', OUT / 'methods_tables.tex')


if __name__ == '__main__':
    main()
