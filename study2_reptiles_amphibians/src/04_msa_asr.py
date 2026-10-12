"""MAFFT alignments anchored on the zebrafish template, coverage filters and tetrapod-ancestor reconstruction.

Filters (plan section 4): eligible sites = zebrafish AlphaFold pLDDT >= 70; species-gene pair valid if it covers
>= 80% of eligible sites; gene eligible if valid in >= 60% of analysed species and present in >= 2 of the 3 fishes;
species eligible if valid in >= 70% of eligible genes (genes, then species, then genes again).
"""
import copy
import io
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from Bio import Phylo

from common import DATA, RESULTS, TOOLS

MAFFT = TOOLS / 'mafft-win' / 'mafft.bat'
IQTREE = TOOLS / 'iqtree-2.4.0-Windows' / 'bin' / 'iqtree2.exe'
SEQ, ALN, ASR = DATA / 'sequences', DATA / 'alignments', DATA / 'asr'
ALN.mkdir(exist_ok=True)
ASR.mkdir(exist_ok=True)
FISH = ['danio_rerio', 'lepisosteus_oculatus', 'latimeria_chalumnae']
AA = 'ARNDCQEGHILKMFPSTWYV'


def read_fasta(text):
    seqs, name = {}, None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def plddt(pdb_path):
    values = {}
    for line in open(pdb_path, encoding='utf-8'):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            values[int(line[22:26])] = float(line[60:66])
    return np.array([values[k] for k in sorted(values)])


def align(gene):
    """MAFFT, then keep the columns where zebrafish has a residue."""
    path = ALN / f'{gene}.fasta'
    if path.exists():
        return gene, 'ok'
    p = subprocess.run([str(MAFFT), '--quiet', '--auto', '--thread', '1', str(SEQ / f'{gene}.fasta')],
                       capture_output=True, text=True)
    msa = read_fasta(p.stdout)
    if 'danio_rerio' not in msa:
        return gene, 'mafft_failed'
    cols = [i for i, c in enumerate(msa['danio_rerio']) if c != '-']
    with open(path, 'w', encoding='utf-8') as f:
        for sp, s in msa.items():
            f.write(f'>{sp}\n{"".join(s[i] for i in cols).upper()}\n')
    return gene, 'ok'


def prune(tree, keep):
    t = copy.deepcopy(tree)
    for tip in list(t.get_terminals()):
        if tip.name not in keep:
            t.prune(tip)
    for c in t.find_clades():
        c.comment, c.confidence = None, None
        if not c.is_terminal():
            c.name = None
    buf = io.StringIO()
    Phylo.write(t, buf, 'newick', plain=False)
    return buf.getvalue()


def ancestor_node(treefile, ingroup):
    tree = Phylo.read(treefile, 'newick')
    tips = {t.name for t in tree.get_terminals()}
    ing, out = set(ingroup) & tips, tips - set(ingroup)
    for clade in tree.find_clades(order='preorder'):
        below = {t.name for t in clade.get_terminals()}
        if below == ing:
            return clade.name
        if below == out:
            path = tree.get_path(clade)
            return (path[-2] if len(path) > 1 else tree.root).name
    raise ValueError('ingroup not separable from outgroups')


def run_asr(args):
    gene, names, ingroup, nwk = args
    folder = ASR / str(gene)
    folder.mkdir(exist_ok=True)
    state = folder / f'{gene}.state'
    if not state.exists():
        seqs = read_fasta((ALN / f'{gene}.fasta').read_text(encoding='utf-8'))
        with open(folder / 'aln.fasta', 'w', encoding='utf-8') as f:
            for sp in names:
                f.write(f'>{sp}\n{seqs[sp].replace("U", "C").replace("*", "-")}\n')
        (folder / 'tree.nwk').write_text(nwk)
        p = subprocess.run([str(IQTREE), '-s', 'aln.fasta', '-te', 'tree.nwk', '-m', 'LG+G4', '-asr', '-nt', '1',
                            '-pre', str(gene), '-quiet', '-redo', '-seed', '20261009'], cwd=folder, capture_output=True, text=True)
        if p.returncode != 0 or not state.exists():
            return {'human_gene': gene, 'status': 'iqtree_failed', 'stderr': (p.stderr or p.stdout)[-300:]}
    node = ancestor_node(folder / f'{gene}.treefile', ingroup)
    st = pd.read_csv(state, sep='\t', comment='#')
    probs = st[st.Node == node].sort_values('Site')[[f'p_{a}' for a in AA]].to_numpy()
    np.savez_compressed(folder / 'ancestor.npz', probs=probs, node=node)
    return {'human_gene': gene, 'status': 'ok', 'node': node, 'sites': len(probs),
            'mean_max_posterior': float(probs.max(axis=1).mean())}


def main():
    genes = pd.read_csv(DATA / 'gene_selection_audit.csv')
    genes = genes[genes.gene_set != 'not_used']
    with ThreadPoolExecutor(max_workers=14) as pool:
        st = dict(pool.map(align, genes.human_gene))
    print('MAFFT:', pd.Series(st).value_counts().to_dict(), flush=True)
    sp = pd.read_csv(DATA / 'species_candidates.csv')
    analysed = set(sp[sp.candidate].binomial)
    rows = []
    for g in genes.itertuples():
        if st.get(g.human_gene) != 'ok':
            continue
        conf = plddt(DATA / 'alphafold' / f'{g.zf_uniprot}.pdb')
        seqs = read_fasta((ALN / f'{g.human_gene}.fasta').read_text(encoding='utf-8'))
        if len(conf) != len(seqs['danio_rerio']):
            continue
        el = conf >= 70
        for name, s in seqs.items():
            cov = float(np.mean([c not in '-X' for c, e in zip(s, el) if e])) if el.any() else 0.0
            rows.append({'human_gene': g.human_gene, 'species': name, 'coverage': cov})
    cov = pd.DataFrame(rows)
    cov['valid'] = cov.coverage >= 0.80

    def gene_filter(species_set):
        c = cov[cov.valid]
        frac = c[c.species.isin(species_set)].groupby('human_gene').size() / len(species_set)
        fish = c[c.species.isin(FISH)].groupby('human_gene').size()
        return set(frac[frac >= 0.60].index) & set(fish[fish >= 2].index)

    g1 = gene_filter(analysed)
    sfrac = cov[cov.valid & cov.human_gene.isin(g1)].groupby('species').size() / len(g1)
    final_species = set(sfrac[sfrac >= 0.70].index) & analysed
    final_genes = gene_filter(final_species)
    sp['eligible'] = sp.binomial.isin(final_species)
    sp['gene_coverage_fraction'] = sp.binomial.map(sfrac)
    sp.to_csv(DATA / 'species_final.csv', index=False)
    cov.to_csv(RESULTS / 'pair_coverage.csv.gz', index=False)
    print(f'genes eligible {len(final_genes)} | species {len(final_species)} of {len(analysed)}', flush=True)

    tree = Phylo.read(DATA / 'timetree_upgma.nwk', 'newick')
    ingroup = sorted(final_species | {'homo_sapiens'})
    tasks = []
    for g in sorted(final_genes):
        valid = set(cov[(cov.human_gene == g) & cov.valid].species)
        names = sorted((valid & (final_species | set(FISH) | {'homo_sapiens'})) | {'danio_rerio'})
        tasks.append((g, names, ingroup, prune(tree, set(names))))
    with ThreadPoolExecutor(max_workers=16) as pool:
        res = pd.DataFrame(list(pool.map(run_asr, tasks)))
    res = res.merge(genes[['human_gene', 'symbol', 'zf_uniprot', 'gene_set', 'genage']], on='human_gene', how='left')
    res.to_csv(RESULTS / 'asr_summary.csv', index=False)
    print(res.status.value_counts().to_dict(), '| mean max posterior', round(res.mean_max_posterior.mean(), 3))


if __name__ == '__main__':
    main()
