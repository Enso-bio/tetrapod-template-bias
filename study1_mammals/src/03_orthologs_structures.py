"""Human-anchored one-to-one ortholog alignments (Ensembl 116 FTP) and AlphaFold DB human models.

GenAge and the 12 previous-study families are always processed; background genes are taken in their fixed random
order until 2,000 pass the universe criteria (amendments E2, E3).
"""
import gzip
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices

from common import DATA, RAW, fetch, fetch_json, _record

N_BACKGROUND = 2000
ALN = DATA / 'alignments'
ALN.mkdir(exist_ok=True)
STD = set('ACDEFGHIKLMNPQRSTVWYX')


def read_pep(path):
    """protein id (no version) -> (gene id, symbol, sequence)."""
    out, pid, meta, seq = {}, None, None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if pid:
                    out[pid] = (*meta, ''.join(seq))
                parts = line[1:].split()
                pid = parts[0].split('.')[0]
                kv = dict(p.split(':', 1) for p in parts[1:] if ':' in p)
                meta = (kv.get('gene', '').split('.')[0], kv.get('gene_symbol', ''))
                seq = []
            else:
                seq.append(line.strip())
    if pid:
        out[pid] = (*meta, ''.join(seq))
    return out


def one2one_table(species):
    rows = []
    files = [RAW / 'compara' / 'homo_sapiens.homologies.tsv.gz'] + \
            [RAW / 'compara' / f'{s}.homologies.tsv.gz' for s in species if (RAW / 'compara' / f'{s}.homologies.tsv.gz').exists()]
    for path in files:
        with gzip.open(path, 'rt') as f:
            header = f.readline().rstrip('\n').split('\t')
            for line in f:
                x = dict(zip(header, line.rstrip('\n').split('\t')))
                if x['homology_type'] != 'ortholog_one2one':
                    continue
                if x['species'] == 'homo_sapiens' and x['homology_species'] in species:
                    h, o = ('gene_stable_id', 'protein_stable_id'), ('homology_gene_stable_id', 'homology_protein_stable_id')
                    other = x['homology_species']
                elif x['homology_species'] == 'homo_sapiens' and x['species'] in species:
                    o, h = ('gene_stable_id', 'protein_stable_id'), ('homology_gene_stable_id', 'homology_protein_stable_id')
                    other = x['species']
                else:
                    continue
                rows.append({'human_gene': x[h[0]], 'human_protein': x[h[1]].split('.')[0], 'species': other,
                             'gene': x[o[0]], 'protein': x[o[1]].split('.')[0], 'goc_score': x['goc_score'],
                             'wga_coverage': x['wga_coverage'], 'high_confidence': x['is_high_confidence']})
    return pd.DataFrame(rows).drop_duplicates(['human_gene', 'species'])


_aligner = None


def anchor(args):
    """Semi-global BLOSUM62 alignment; returns partner residues at each human position ('-' if unaligned)."""
    global _aligner
    human, other = args
    if _aligner is None:
        _aligner = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                                   open_gap_score=-10, extend_gap_score=-0.5,
                                   end_insertion_score=0, end_deletion_score=0)
    clean = lambda s: ''.join(c if c in STD else ('C' if c == 'U' else 'X') for c in s.upper())
    aln = _aligner.align(clean(human), clean(other))[0]
    out = ['-'] * len(human)
    for (h0, h1), (o0, o1) in zip(*aln.aligned):
        for i, j in zip(range(h0, h1), range(o0, o1)):
            out[i] = other[j]
    return ''.join(out)


def af_model(row):
    try:
        info = fetch_json(f'https://alphafold.ebi.ac.uk/api/prediction/{row.accession}',
                          RAW / 'alphafold_api' / f'{row.accession}.json', ok_404=True)
        if not info:
            return {'accession': row.accession, 'af_status': 'no_alphafold_model'}
        full = [m for m in info if m.get('sequence', m.get('uniprotSequence')) == row.sequence
                and m.get('sequenceStart', m.get('uniprotStart', 1)) == 1
                and m.get('sequenceEnd', m.get('uniprotEnd')) == row.length]
        if not full:
            return {'accession': row.accession, 'af_status': 'no_full_length_identical_model'}
        m = max(full, key=lambda m: m.get('latestVersion', 0))
        fetch(m['pdbUrl'], DATA / 'alphafold' / f'{row.accession}.pdb')
        return {'accession': row.accession, 'af_status': 'ok', 'af_url': m['pdbUrl'], 'af_version': m.get('latestVersion'),
                'af_global_plddt': m.get('globalMetricValue')}
    except Exception as e:
        return {'accession': row.accession, 'af_status': 'error', 'af_error': repr(e)[:300]}


def main():
    human_file = RAW / 'compara' / 'homo_sapiens.homologies.tsv.gz'
    if not any(human_file.name in l for l in open(DATA / 'provenance.jsonl', encoding='utf-8')):
        _record('https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/homo_sapiens/'
                'Compara.116.protein_default.homologies.tsv.gz', human_file)
    sp = pd.read_csv(DATA / 'species_candidates.csv')
    species = set(sp.ensembl_name) - {'homo_sapiens'}
    o2o = one2one_table(species)
    o2o.to_csv(DATA / 'one2one_orthologs.csv.gz', index=False)
    hpep = read_pep(RAW / 'pep' / 'homo_sapiens.pep.all.fa.gz')
    canon = o2o.drop_duplicates('human_gene')[['human_gene', 'human_protein']]
    canon['symbol'] = canon.human_protein.map(lambda p: hpep.get(p, ('', '', ''))[1])
    canon['human_seq'] = canon.human_protein.map(lambda p: hpep.get(p, ('', '', ''))[2])
    canon = canon[canon.symbol != '']
    canon = canon[~canon.symbol.duplicated(keep=False)]

    uni = pd.read_csv(DATA / 'gene_universe_candidates.csv')
    uni = uni.merge(canon, on='symbol', how='left')
    uni['ensembl_status'] = 'ok'
    uni.loc[uni.human_gene.isna(), 'ensembl_status'] = 'no_compara_one2one_gene'
    uni.loc[uni.human_gene.notna() & (uni.human_seq != uni.sequence), 'ensembl_status'] = 'ensembl_uniprot_sequence_differs'

    # Fixed sets always; background in fixed random order until N_BACKGROUND pass (Ensembl identity + AF model).
    fixed = uni[uni.genage | uni.old_family]
    bg = uni[~uni.genage & ~uni.old_family].sort_values('background_rank')
    bg_ok = bg[bg.ensembl_status == 'ok']
    todo = pd.concat([fixed[fixed.ensembl_status == 'ok'], bg_ok.iloc[:int(N_BACKGROUND * 1.15)]])
    with ThreadPoolExecutor(max_workers=8) as pool:
        af = pd.DataFrame(list(pool.map(af_model, todo.itertuples())))
    uni = uni.merge(af, on='accession', how='left')
    uni['status'] = uni.ensembl_status.where(uni.ensembl_status != 'ok', uni.af_status.fillna('not_processed'))
    passed_bg = uni[~uni.genage & ~uni.old_family & (uni.status == 'ok')].sort_values('background_rank')
    assert len(passed_bg) >= N_BACKGROUND, 'extend the background batch'
    last_rank = passed_bg.background_rank.iloc[N_BACKGROUND - 1]
    uni['gene_set'] = 'not_used'
    uni.loc[uni.genage & (uni.status == 'ok'), 'gene_set'] = 'genage'
    uni.loc[uni.old_family & ~uni.genage & (uni.status == 'ok'), 'gene_set'] = 'old_family'
    uni.loc[~uni.genage & ~uni.old_family & (uni.status == 'ok') & (uni.background_rank <= last_rank), 'gene_set'] = 'background'
    uni.loc[~uni.genage & ~uni.old_family & (uni.background_rank > last_rank), 'status'] = 'beyond_background_quota'
    sel = uni[uni.gene_set != 'not_used']

    # Alignments
    peps = {s: read_pep(RAW / 'pep' / f'{s}.pep.all.fa.gz') for s in sorted(species)}
    o2o_sel = o2o[o2o.human_gene.isin(sel.human_gene)]
    jobs, keys = [], []
    for r in o2o_sel.itertuples():
        seq = peps[r.species].get(r.protein, (None, None, None))[2]
        if seq:
            human = sel.loc[sel.human_gene == r.human_gene, 'sequence'].iloc[0]
            jobs.append((human, seq.rstrip('*')))
            keys.append((r.human_gene, r.species))
    with ProcessPoolExecutor(max_workers=16) as pool:
        anchored = list(pool.map(anchor, jobs, chunksize=200))
    by_gene = {}
    for (g, s), a in zip(keys, anchored):
        by_gene.setdefault(g, {})[s] = a
    n_orth = {}
    for r in sel.itertuples():
        seqs = by_gene.get(r.human_gene, {})
        n_orth[r.symbol] = len(seqs)
        with open(ALN / f'{r.symbol}.fasta', 'w', encoding='utf-8') as f:
            f.write(f'>homo_sapiens\n{r.sequence}\n')
            for s in sorted(seqs):
                f.write(f'>{s}\n{seqs[s]}\n')
    uni['n_one2one_aligned'] = uni.symbol.map(n_orth)
    uni.drop(columns=['sequence', 'human_seq']).to_csv(DATA / 'gene_selection_audit.csv', index=False)
    print(uni.groupby(['gene_set', 'status']).size().to_string())
    print('aligned pairs', len(jobs), '| species pep files', len(peps))


if __name__ == '__main__':
    main()
