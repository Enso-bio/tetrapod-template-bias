"""Gene universe, NCBI one-to-one orthologs, representative proteins and zebrafish AlphaFold templates.

GenAge human genes are always processed; background genes are taken in a fixed random order (seed 20261009) until
2,000 pass the universe criteria (plan section 3).
"""
import gzip
import io
import re
import zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode

import numpy as np
import pandas as pd

from common import DATA, RAW, SEED, fetch, fetch_json

N_BACKGROUND = 2000
SEQ = DATA / 'sequences'
SEQ.mkdir(exist_ok=True)
ZEBRAFISH = 7955


def read_genome(acc):
    """protein accession -> sequence; protein -> GeneID; set of MANE/RefSeq Select proteins."""
    z = zipfile.ZipFile(RAW / 'genomes' / f'{acc}.zip')
    names = z.namelist()
    faa = next(n for n in names if n.endswith('protein.faa'))
    gff = next(n for n in names if n.endswith('.gff'))
    seqs, pid, buf = {}, None, []
    with z.open(faa) as f:
        for line in io.TextIOWrapper(f, encoding='utf-8'):
            if line.startswith('>'):
                if pid:
                    seqs[pid] = ''.join(buf)
                pid, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if pid:
        seqs[pid] = ''.join(buf)
    gene_of, select, symbol = {}, set(), {}
    with z.open(gff) as f:
        for line in io.TextIOWrapper(f, encoding='utf-8'):
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            attr = parts[8]
            p = re.search(r'protein_id=([^;]+)', attr)
            g = re.search(r'GeneID:(\d+)', attr)
            if p and g:
                gene_of[p.group(1)] = int(g.group(1))
                s = re.search(r'(?:^|;)gene=([^;]+)', attr)
                if s:
                    symbol[int(g.group(1))] = s.group(1)
                if 'MANE Select' in attr or 'RefSeq Select' in attr:
                    select.add(p.group(1))
    return seqs, gene_of, select, symbol


def representative(seqs, gene_of, select):
    """GeneID -> (protein accession, sequence): Select isoform if annotated, else the longest."""
    by_gene = defaultdict(list)
    for p, g in gene_of.items():
        if p in seqs:
            by_gene[g].append(p)
    rep = {}
    for g, ps in by_gene.items():
        chosen = [p for p in ps if p in select] or ps
        p = max(chosen, key=lambda x: (len(seqs[x]), x))
        rep[g] = (p, seqs[p])
    return rep, by_gene


def main():
    sp = pd.read_csv(DATA / 'species_candidates.csv')
    sp = sp[sp.candidate]
    out = pd.read_csv(DATA / 'outgroups.csv')
    genomes = pd.concat([sp[['binomial', 'taxon_id', 'accession']], out[['binomial', 'taxon_id', 'accession']]])
    taxa = {int(t): b for t, b in zip(genomes.taxon_id, genomes.binomial)}

    # NCBI orthologs relative to human GeneIDs
    orth = defaultdict(lambda: defaultdict(set))          # human GeneID -> taxid -> {GeneID}
    with gzip.open(RAW / 'gene_orthologs.gz', 'rt') as f:
        next(f)
        for line in f:
            t1, g1, rel, t2, g2 = line.rstrip('\n').split('\t')
            if rel != 'Ortholog':
                continue
            if t1 == '9606' and int(t2) in taxa:
                orth[int(g1)][int(t2)].add(int(g2))
            elif t2 == '9606' and int(t1) in taxa:
                orth[int(g2)][int(t1)].add(int(g1))
    print('human genes with orthologs:', len(orth), flush=True)

    reps, isoforms = {}, {}
    for r in genomes.itertuples():
        seqs, gene_of, select, symbol = read_genome(r.accession)
        if r.binomial == 'homo_sapiens':
            human_symbol = symbol
        reps[r.binomial], iso = representative(seqs, gene_of, select)
        if int(r.taxon_id) == ZEBRAFISH:
            isoforms = {g: [(p, seqs[p]) for p in ps] for g, ps in iso.items()}
        print('genome', r.binomial, len(reps[r.binomial]), 'genes', flush=True)
    human = reps['homo_sapiens']

    # GenAge and universe
    gz = zipfile.ZipFile(RAW / 'genage_human.zip')
    genage = pd.read_csv(io.BytesIO(gz.read([n for n in gz.namelist() if n.endswith('.csv')][0])))
    genage_ids = set(pd.to_numeric(genage['entrez gene id'], errors='coerce').dropna().astype(int))
    uni = []
    for hg, per in orth.items():
        if hg not in human or len(per.get(ZEBRAFISH, ())) != 1:
            continue
        hp, hs = human[hg]
        if not 50 <= len(hs) <= 1000:
            continue
        uni.append({'human_gene': hg, 'symbol': human_symbol.get(hg, ''), 'human_protein': hp, 'human_length': len(hs), 'genage': hg in genage_ids,
                    'zebrafish_gene': next(iter(per[ZEBRAFISH]))})
    uni = pd.DataFrame(uni)
    rng = np.random.default_rng(SEED)
    uni['background_rank'] = np.where(uni.genage, -1, 0)
    bg = uni.index[~uni.genage]
    uni.loc[bg, 'background_rank'] = rng.permutation(len(bg))

    # Zebrafish UniProt entries -> AlphaFold models
    q = {'query': '(organism_id:7955)', 'format': 'tsv', 'fields': 'accession,reviewed,sequence'}
    up = pd.read_csv(io.BytesIO(fetch('https://rest.uniprot.org/uniprotkb/stream?' + urlencode(q),
                                      RAW / 'uniprot_danio_rerio.tsv', timeout=900)), sep='\t')
    up = up.sort_values(['Reviewed', 'Entry'], ascending=[False, True]).drop_duplicates('Sequence')
    seq_to_acc = up.set_index('Sequence').Entry

    def template(row):
        cands = sorted(isoforms.get(row.zebrafish_gene, []), key=lambda x: -len(x[1]))
        for p, s in cands:
            acc = seq_to_acc.get(s)
            if acc is None:
                continue
            info = fetch_json(f'https://alphafold.ebi.ac.uk/api/prediction/{acc}', RAW / 'alphafold_api' / f'{acc}.json', ok_404=True)
            full = [m for m in (info or []) if m.get('sequence', m.get('uniprotSequence')) == s]
            if full:
                m = max(full, key=lambda m: m.get('latestVersion', 0))
                fetch(m['pdbUrl'], DATA / 'alphafold' / f'{acc}.pdb')
                return {'human_gene': row.human_gene, 'zf_protein': p, 'zf_uniprot': acc, 'zf_length': len(s),
                        'af_version': m.get('latestVersion'), 'status': 'ok'}
        return {'human_gene': row.human_gene, 'status': 'no_identical_uniprot_with_alphafold'}

    fixed = uni[uni.genage]
    order = uni[~uni.genage].sort_values('background_rank')
    results, passed, start = [], 0, 0
    with ThreadPoolExecutor(max_workers=4) as pool:
        results += list(pool.map(template, fixed.itertuples()))
        while passed < N_BACKGROUND and start < len(order):
            part = list(pool.map(template, order.iloc[start:start + 500].itertuples()))
            results += part
            passed += sum(r['status'] == 'ok' for r in part)
            start += 500
            print('background processed', start, 'passed', passed, flush=True)
    res = pd.DataFrame(results)
    uni = uni.merge(res, on='human_gene', how='left')
    ok_bg = uni[~uni.genage & (uni.status == 'ok')].sort_values('background_rank')
    last = ok_bg.background_rank.iloc[N_BACKGROUND - 1]
    uni['gene_set'] = 'not_used'
    uni.loc[uni.genage & (uni.status == 'ok'), 'gene_set'] = 'genage'
    uni.loc[~uni.genage & (uni.status == 'ok') & (uni.background_rank <= last), 'gene_set'] = 'background'
    sel = uni[uni.gene_set != 'not_used']

    # Raw ortholog sequences per gene (template isoform for zebrafish, representative isoform for the rest)
    zf_seq = {p: s for ps in isoforms.values() for p, s in ps}
    n_orth = {}
    for r in sel.itertuples():
        recs = [('danio_rerio', zf_seq[r.zf_protein]), ('homo_sapiens', human[r.human_gene][1])]
        for t, b in taxa.items():
            if b in ('danio_rerio', 'homo_sapiens'):
                continue
            ids = orth[r.human_gene].get(t, set())
            if len(ids) == 1 and next(iter(ids)) in reps[b]:
                recs.append((b, reps[b][next(iter(ids))][1]))
        n_orth[r.human_gene] = len(recs)
        with open(SEQ / f'{r.human_gene}.fasta', 'w', encoding='utf-8') as f:
            for b, s in recs:
                f.write(f'>{b}\n{s.rstrip("*")}\n')
    uni['n_sequences'] = uni.human_gene.map(n_orth)
    uni.to_csv(DATA / 'gene_selection_audit.csv', index=False)
    print(uni.groupby(['gene_set', 'status'], dropna=False).size().to_string())


if __name__ == '__main__':
    main()
