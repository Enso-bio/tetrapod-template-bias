"""Species set: NCBI-annotated reptile and amphibian reference genomes, AnAge phenotypes, body mass, TimeTree tree."""
import io
import itertools
import json
import time
import xml.etree.ElementTree as ET
import zipfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.cluster.hierarchy import average, cophenet, to_tree
from scipy.spatial.distance import squareform

from common import DATA, RAW, ROOT, binomial, fetch, fetch_json, record

GROUPS = {'Testudines': 8459, 'Crocodylia': 1294634, 'Lepidosauria': 8504, 'Amphibia': 8292}
OUTGROUPS = {'Danio rerio': 7955, 'Lepisosteus oculatus': 7918, 'Latimeria chalumnae': 7897, 'Homo sapiens': 9606}
DATASETS = 'https://api.ncbi.nlm.nih.gov/datasets/v2'


def genome_reports(taxon, tag):
    out, token, page = [], None, 0
    while True:
        params = {'filters.has_annotation': 'true', 'filters.reference_only': 'true', 'page_size': 500}
        if token:
            params['page_token'] = token
        js = fetch_json(f'{DATASETS}/genome/taxon/{taxon}/dataset_report', RAW / 'ncbi_reports' / f'{tag}_{page}.json',
                        params=params, timeout=300)
        out += js.get('reports', [])
        token, page = js.get('next_page_token'), page + 1
        if not token:
            return out


def taxonomy(taxid):
    """Species-rank taxid and NCBI names (scientific + synonyms) for a taxon."""
    path = RAW / 'ncbi_taxonomy' / f'{taxid}.xml'
    cached = path.exists()
    xml = fetch(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxid}&retmode=xml', path)
    if not cached:
        time.sleep(0.4)
    t = ET.fromstring(xml).find('Taxon')
    species_id, species_name = int(t.findtext('TaxId')), t.findtext('ScientificName')
    if t.findtext('Rank') != 'species':
        for lin in t.find('LineageEx') or []:
            if lin.findtext('Rank') == 'species':
                species_id, species_name = int(lin.findtext('TaxId')), lin.findtext('ScientificName')
    names = [t.findtext('ScientificName')]
    other = t.find('OtherNames')
    if other is not None:
        for tag in ('Synonym', 'EquivalentName', 'GenbankSynonym'):
            names += [e.text for e in other.findall(tag) if e.text]
    return species_id, species_name, names


def timetree_age(a, b):
    """TimeTree 5 pairwise divergence time (Ma) between two NCBI taxids; None if unavailable."""
    lo, hi = sorted((a, b))
    path = RAW / 'timetree' / f'{lo}_{hi}.csv'
    try:
        txt = fetch(f'https://timetree.temple.edu/api/pairwise/{lo}/{hi}', path, ok_404=True, timeout=60)
    except Exception:
        return None
    if not txt:
        return None
    rows = list(pd.read_csv(io.BytesIO(txt)).itertuples())
    if not rows:
        return None
    age = getattr(rows[0], 'adjusted_age', np.nan)
    return float(age) if pd.notna(age) and float(age) > 0 else None


def newick(node, labels):
    def rec(n, parent_h):
        if n.is_leaf():
            return f'{labels[n.id]}:{parent_h:.6f}'
        h = n.dist / 2
        left, right = rec(n.left, h), rec(n.right, h)
        return f'({left},{right}):{parent_h - h:.6f}'
    h = node.dist / 2
    return f'({rec(node.left, h)},{rec(node.right, h)});'


def main():
    # AnAge (same archive as the mammal study, downloaded 8 Oct 2026)
    z = zipfile.ZipFile(RAW / 'anage.zip')
    an = pd.read_csv(io.BytesIO(z.read('anage_data.txt')), sep='\t')
    an['binomial'] = (an['Genus'].str.strip() + '_' + an['Species'].str.strip()).str.lower()
    for col in ['Maximum longevity (yrs)', 'Adult weight (g)']:
        an[col] = pd.to_numeric(an[col], errors='coerce')
    ok = ((an['Maximum longevity (yrs)'] > 0) & an['Data quality'].str.lower().isin(['acceptable', 'high'])
          & an['Sample size'].str.lower().isin(['small', 'medium', 'large', 'huge']))
    pheno = an[ok & an['Class'].isin(['Reptilia', 'Amphibia'])].drop_duplicates('binomial').set_index('binomial')

    # Mass supplements
    ab = pd.read_csv(zipfile.ZipFile(RAW / 'AmphiBIO_v1.zip').open('AmphiBIO_v1.csv'), encoding='latin-1')
    ab['binomial'] = ab.Species.str.lower().str.replace(' ', '_')
    ab_mass = ab.dropna(subset=['Body_mass_g']).drop_duplicates('binomial').set_index('binomial').Body_mass_g
    rt = pd.read_excel(RAW / 'ReptTraits_v1-2.xlsx', sheet_name='Data')
    rt['binomial'] = rt.Species.str.strip().str.lower().str.replace(' ', '_')
    rt = rt.drop_duplicates('binomial').set_index('binomial')
    logged = (DATA / 'provenance.jsonl').read_text(encoding='utf-8') if (DATA / 'provenance.jsonl').exists() else ''
    for url, f in [('https://ndownloader.figshare.com/files/8828578', 'AmphiBIO_v1.zip'),
                   ('https://ndownloader.figshare.com/files/45408133', 'ReptTraits_v1-2.xlsx'),
                   ('https://genomics.senescence.info/species/dataset.zip (copied from PROTEINAS_LONGEVIDAD)', 'anage.zip'),
                   ('https://genomics.senescence.info/genes/human_genes.zip (copied from PROTEINAS_LONGEVIDAD)', 'genage_human.zip')]:
        if f not in logged:
            record(url, RAW / f)

    # NCBI annotated reference genomes
    rows = []
    for group, tax in GROUPS.items():
        for rep in genome_reports(tax, group):
            org = rep.get('organism', {})
            ann = rep.get('annotation_info', {})
            rows.append({'group': group, 'organism': org.get('organism_name'), 'taxon_id': int(org.get('tax_id')),
                         'accession': rep.get('accession'), 'annotation_provider': ann.get('provider'),
                         'annotation_name': ann.get('name'), 'busco_complete': (ann.get('busco') or {}).get('complete')})
    for name, tax in OUTGROUPS.items():
        for rep in genome_reports(tax, name.replace(' ', '_')):
            ann = rep.get('annotation_info', {})
            rows.append({'group': 'outgroup' if tax != 9606 else 'human', 'organism': rep['organism']['organism_name'],
                         'taxon_id': tax, 'accession': rep.get('accession'), 'annotation_provider': ann.get('provider'),
                         'annotation_name': ann.get('name'), 'busco_complete': (ann.get('busco') or {}).get('complete')})
    g = pd.DataFrame(rows)
    # One genome per organism: prefer RefSeq (GCF) annotation, then higher BUSCO.
    g['is_refseq'] = g.accession.str.startswith('GCF_')
    g = g.sort_values(['organism', 'is_refseq', 'busco_complete'], ascending=[True, False, False]).drop_duplicates('organism')
    tx = {t: taxonomy(t) for t in g.taxon_id.unique()}
    g['species_taxid'] = g.taxon_id.map(lambda t: tx[t][0])
    g['species_name'] = g.taxon_id.map(lambda t: tx[t][1])
    g['binomial'] = g.species_name.map(binomial)
    g = g.sort_values(['binomial', 'is_refseq', 'busco_complete'], ascending=[True, False, False]).drop_duplicates('binomial')

    # Phenotype matching (exact binomial, then NCBI synonyms of the same taxon; amendment-E1 rule of the mammal study)
    def match(row):
        if row.binomial in pheno.index:
            return row.binomial, 'exact'
        alts = {binomial(n) for n in tx[row.taxon_id][2] + [row.organism] if len(n.split()) == 2}
        hits = sorted(a for a in alts if a in pheno.index and a != row.binomial)
        return (hits[0], 'ncbi_synonym') if len(hits) == 1 else (None, None)
    m = g[g.group.isin(GROUPS)].apply(match, axis=1, result_type='expand')
    g.loc[m.index, 'anage_binomial'], g.loc[m.index, 'match_type'] = m[0], m[1]
    ing = g[g.group.isin(GROUPS) & g.anage_binomial.notna()].copy()
    p = pheno.loc[ing.anage_binomial]
    ing['order'] = p['Order'].to_numpy()
    ing['max_lifespan_yr'] = p['Maximum longevity (yrs)'].to_numpy()
    ing['data_quality'] = p['Data quality'].to_numpy()
    ing['mass_anage_g'] = p['Adult weight (g)'].to_numpy()
    ing['mass_g'], ing['mass_source'] = ing.mass_anage_g, np.where(ing.mass_anage_g > 0, 'AnAge', None)
    for i, r in ing[~(ing.mass_g > 0)].iterrows():
        keys = [r.binomial, r.anage_binomial]
        if r.group == 'Amphibia':
            v = next((ab_mass[k] for k in keys if k in ab_mass.index), np.nan)
            src = 'AmphiBIO'
        else:
            v = next((rt.loc[k, 'Maximum body mass (g)'] for k in keys if k in rt.index), np.nan)
            src = 'ReptTraits'
        v = pd.to_numeric(v, errors='coerce')
        if v > 0:
            ing.loc[i, ['mass_g', 'mass_source']] = [float(v), src]
    for col, src in [('Mean Tb', 'body_temp_c'), ('Mean Annual Temperature (°C)', 'range_mean_temp_c')]:
        ing[src] = [pd.to_numeric(rt.loc[k, col], errors='coerce') if k in rt.index else np.nan
                    for k in ing.anage_binomial]
    ing['mass_ok'] = ing.mass_g > 0
    out = g[g.group.isin(['outgroup', 'human'])].copy()

    # TimeTree pairwise ages for candidates with phenotype and mass, plus outgroups and human
    cand = pd.concat([ing[ing.mass_ok], out], ignore_index=True)
    ids = cand.species_taxid.tolist()
    pairs = list(itertools.combinations(range(len(ids)), 2))
    with ThreadPoolExecutor(max_workers=4) as pool:
        ages = list(pool.map(lambda ij: timetree_age(ids[ij[0]], ids[ij[1]]), pairs))
    A = np.full((len(ids), len(ids)), np.nan)
    for (i, j), a in zip(pairs, ages):
        A[i, j] = A[j, i] = a if a is not None else np.nan
    np.fill_diagonal(A, 0)
    missing = np.isnan(A).sum(axis=1)
    cand['timetree_missing_pairs'] = missing
    keep = missing == 0
    while not keep.all() and np.isnan(A[np.ix_(keep, keep)]).any():
        worst = np.argmax(np.where(keep, np.isnan(A).sum(axis=1), -1))
        keep[worst] = False
    cand['in_tree'] = keep
    sub = A[np.ix_(keep, keep)]
    labels = cand.loc[keep, 'binomial'].tolist()
    Z = average(squareform(2 * sub, checks=False))
    coph = cophenet(Z) / 2
    dev = np.abs(coph - squareform(sub, checks=False))
    nwk = newick(to_tree(Z), labels)
    (DATA / 'timetree_upgma.nwk').write_text(nwk)
    tree = Phylo.read(io.StringIO(nwk), 'newick')
    assert len(tree.get_terminals()) == len(labels)
    json.dump({'n_taxa': len(labels), 'max_abs_deviation_Ma': float(dev.max()),
               'mean_abs_deviation_Ma': float(dev.mean()), 'root_age_Ma': float(sub.max())},
              open(DATA / 'timetree_upgma_summary.json', 'w'), indent=2)
    pd.DataFrame(A, index=cand.binomial, columns=cand.binomial).to_csv(DATA / 'timetree_pairwise_ages.csv')

    g.to_csv(DATA / 'genome_audit.csv', index=False)
    ing = ing.merge(cand[['binomial', 'in_tree', 'timetree_missing_pairs']], on='binomial', how='left')
    ing['candidate'] = ing.mass_ok & ing.in_tree.fillna(False).astype(bool)
    ing.to_csv(DATA / 'species_candidates.csv', index=False)
    out = out.merge(cand[['binomial', 'in_tree']], on='binomial', how='left')
    out.to_csv(DATA / 'outgroups.csv', index=False)
    print('annotated reference genomes:', g.group.value_counts().to_dict())
    print('with AnAge lifespan:', ing.group.value_counts().to_dict())
    print('with mass:', ing[ing.mass_ok].group.value_counts().to_dict(), '| mass sources:', ing.mass_source.value_counts().to_dict())
    print('candidates (phenotype + mass + TimeTree):', ing[ing.candidate].group.value_counts().to_dict(), '| total', int(ing.candidate.sum()))
    print('outgroups in tree:', out[['organism', 'accession', 'in_tree']].to_dict('records'))
    print('TimeTree UPGMA:', json.load(open(DATA / 'timetree_upgma_summary.json')))
    print(ing[~ing.candidate][['organism', 'group', 'mass_ok', 'in_tree', 'timetree_missing_pairs']].to_string(index=False))


if __name__ == '__main__':
    main()
