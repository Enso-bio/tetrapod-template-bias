"""Eligible sites, gene/species filters and ancestral reconstruction of the placental ancestor with IQ-TREE.

Filters (plan section 4): species-gene pair valid if it covers >=80% of eligible sites (human AlphaFold pLDDT>=70);
gene eligible if valid in >=60% of placental species and in >=1 non-placental outgroup; species eligible if valid in
>=70% of eligible genes. Gene filter is applied with all candidate placentals, then species, then genes again with the
final species set.
"""
import copy
import io
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from Bio import Phylo

from common import DATA, RAW, RESULTS, TOOLS

IQTREE = TOOLS / 'iqtree-2.4.0-Windows' / 'bin' / 'iqtree2.exe'
ALN = DATA / 'alignments'
ASR = DATA / 'asr'
ASR.mkdir(exist_ok=True)
AA = 'ARNDCQEGHILKMFPSTWYV'


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path, encoding='utf-8'):
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


def coverage_table(genes, species):
    rows = []
    for g in genes.itertuples():
        conf = plddt(DATA / 'alphafold' / f'{g.accession}.pdb')
        seqs = read_fasta(ALN / f'{g.symbol}.fasta')
        eligible = conf >= 70
        if len(conf) != len(seqs['homo_sapiens']):
            continue
        for sp in species.ensembl_name:
            s = seqs.get(sp)
            cov = 0.0 if s is None else float(np.mean([c not in '-X' for c, e in zip(s, eligible) if e]))
            rows.append({'symbol': g.symbol, 'species': sp, 'coverage': cov})
    return pd.DataFrame(rows)


def prune_tree(tree, labels_to_names):
    """Copy of `tree` (object or nexus path) restricted to the given tips, renamed, annotations removed."""
    tree = Phylo.read(tree, 'nexus') if not hasattr(tree, 'get_terminals') else copy.deepcopy(tree)
    keep = set(labels_to_names)
    for tip in list(tree.get_terminals()):
        if tip.name not in keep:
            tree.prune(tip)
    for clade in tree.find_clades():
        clade.comment, clade.confidence = None, None
        clade.name = labels_to_names[clade.name] if clade.is_terminal() else None
    return tree


def newick(tree):
    buf = io.StringIO()
    Phylo.write(tree, buf, 'newick', plain=False)
    return buf.getvalue()


def placental_ancestor_node(treefile, placentals):
    """Name of the IQ-TREE internal node that is the MRCA of the placental tips (outgroup-rooted)."""
    tree = Phylo.read(treefile, 'newick')
    tips = {t.name for t in tree.get_terminals()}
    pl = set(placentals) & tips
    out = tips - pl
    for clade in tree.find_clades(order='preorder'):
        below = {t.name for t in clade.get_terminals()}
        if below == pl:
            return clade.name
        if below == out:
            parent = tree.get_path(clade)[-2] if len(tree.get_path(clade)) > 1 else tree.root
            return parent.name
    raise ValueError('placentals not separable from outgroup')


def run_gene(args):
    symbol, species_names, placentals, tree_text = args
    folder = ASR / symbol
    folder.mkdir(exist_ok=True)
    state = folder / f'{symbol}.state'
    if not state.exists():
        seqs = read_fasta(ALN / f'{symbol}.fasta')
        with open(folder / 'aln.fasta', 'w', encoding='utf-8') as f:
            for sp in species_names:
                f.write(f'>{sp}\n{seqs[sp].replace("U", "C").replace("*", "-")}\n')
        (folder / 'tree.nwk').write_text(tree_text)
        cmd = [str(IQTREE), '-s', 'aln.fasta', '-te', 'tree.nwk', '-m', 'LG+G4', '-asr', '-nt', '1',
               '-pre', symbol, '-quiet', '-redo', '-seed', '20261008']
        p = subprocess.run(cmd, cwd=folder, capture_output=True, text=True)
        if p.returncode != 0 or not state.exists():
            return {'symbol': symbol, 'status': 'iqtree_failed', 'stderr': (p.stderr or p.stdout)[-300:]}
    node = placental_ancestor_node(folder / f'{symbol}.treefile', placentals)
    st = pd.read_csv(state, sep='\t', comment='#')
    st = st[st.Node == node].sort_values('Site')
    probs = st[[f'p_{a}' for a in AA]].to_numpy()
    np.savez_compressed(folder / 'ancestor.npz', probs=probs, node=node)
    return {'symbol': symbol, 'status': 'ok', 'node': node, 'sites': len(probs),
            'mean_max_posterior': float(probs.max(axis=1).mean())}


def main():
    sp = pd.read_csv(DATA / 'species_candidates.csv')
    genes = pd.read_csv(DATA / 'gene_selection_audit.csv')
    genes = genes[genes.gene_set != 'not_used']
    cov = coverage_table(genes, sp)
    cov['valid'] = cov.coverage >= 0.80
    placental = set(sp[sp.placental].ensembl_name)
    outgroup = set(sp[~sp.placental].ensembl_name)

    def gene_filter(species_set):
        c = cov[cov.valid & cov.species.isin(species_set | outgroup)]
        frac = c[c.species.isin(species_set)].groupby('symbol').size() / len(species_set)
        has_out = c[c.species.isin(outgroup)].groupby('symbol').size()
        return set(frac[frac >= 0.60].index) & set(has_out.index)

    genes1 = gene_filter(placental)
    sfrac = cov[cov.valid & cov.symbol.isin(genes1)].groupby('species').size() / len(genes1)
    final_species = set(sfrac[sfrac >= 0.70].index) & placental
    final_genes = gene_filter(final_species)
    sp['eligible'] = sp.ensembl_name.isin(final_species | outgroup)
    sp['gene_coverage_fraction'] = sp.ensembl_name.map(sfrac)
    sp.to_csv(DATA / 'species_final.csv', index=False)
    cov.to_csv(RESULTS / 'pair_coverage.csv.gz', index=False)
    print(f'genes eligible {len(final_genes)} of {len(genes)} | placental species {len(final_species)} of '
          f'{len(placental)} | outgroups {len(outgroup)}', flush=True)

    labels = dict(zip(sp.tree_label, sp.ensembl_name))
    base_tree = prune_tree(RAW / 'mcc_node.tre', labels)
    back = {v: v for v in labels.values()}
    tasks = []
    for symbol in sorted(final_genes):
        valid = set(cov[(cov.symbol == symbol) & cov.valid].species)
        names = ['homo_sapiens'] + sorted((valid & (final_species | outgroup)) - {'homo_sapiens'})
        t = prune_tree(base_tree, {k: back[k] for k in names})
        tasks.append((symbol, names, sorted(final_species), newick(t)))
    with ThreadPoolExecutor(max_workers=16) as pool:
        res = list(pool.map(run_gene, tasks))
    res = pd.DataFrame(res)
    res = res.merge(genes[['symbol', 'accession', 'gene_set', 'old_family', 'length']], on='symbol', how='left')
    res.to_csv(RESULTS / 'asr_summary.csv', index=False)
    print(res.status.value_counts().to_dict(), '| mean max posterior', round(res.mean_max_posterior.mean(), 3))


if __name__ == '__main__':
    main()
