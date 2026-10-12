"""AnAge phenotypes, Upham et al. (2019) trees and the Ensembl mammal genome list; species matching."""
import io
import json
import time
import xml.etree.ElementTree as ET
import zipfile

import pandas as pd
from Bio import Phylo

from common import DATA, RAW, fetch, fetch_json, binomial


def ncbi_names(taxon_id):
    """Scientific name plus synonyms/equivalent names recorded by NCBI Taxonomy for one taxon id."""
    path = RAW / 'ncbi_taxonomy' / f'{taxon_id}.xml'
    cached = path.exists()
    xml = fetch(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxon_id}&retmode=xml', path)
    if not cached:
        time.sleep(0.4)  # NCBI limit without API key: 3 requests/s
    taxon = ET.fromstring(xml).find('Taxon')
    names = [taxon.findtext('ScientificName')]
    other = taxon.find('OtherNames')
    if other is not None:
        for tag in ('Synonym', 'EquivalentName', 'GenbankSynonym'):
            names += [e.text for e in other.findall(tag) if e.text]
    return [n.strip() for n in names if n]

UPHAM = ('https://raw.githubusercontent.com/n8upham/MamPhy_v1/740b70189c77fc168be17bd6ade204c755bcb19b/_DATA/'
         'MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_')
TREES = {'mcc_node.tre': 'NDexp_MCC_v2_target.tre',
         'posterior100.trees': 'NDexp_v2_sample100_nexus.trees'}


def main():
    # AnAge
    z = zipfile.ZipFile(io.BytesIO(fetch('https://genomics.senescence.info/species/dataset.zip', RAW / 'anage.zip')))
    an = pd.read_csv(io.BytesIO(z.read('anage_data.txt')), sep='\t')
    an['binomial'] = (an['Genus'].str.strip() + '_' + an['Species'].str.strip()).str.lower()
    for col in ['Maximum longevity (yrs)', 'Adult weight (g)']:
        an[col] = pd.to_numeric(an[col], errors='coerce')
    ok = ((an['Class'] == 'Mammalia') & (an['Maximum longevity (yrs)'] > 0) & (an['Adult weight (g)'] > 0)
          & an['Data quality'].str.lower().isin(['acceptable', 'high'])
          & an['Sample size'].str.lower().isin(['small', 'medium', 'large', 'huge']))
    an['phenotype_eligible'] = ok
    an.to_csv(DATA / 'anage_all.csv', index=False)

    # Trees
    for name, suffix in TREES.items():
        fetch(UPHAM + suffix, RAW / name)
    tree = Phylo.read(RAW / 'mcc_node.tre', 'nexus')
    tips = {binomial(t.name): t.name for t in tree.get_terminals() if not t.name.startswith('_')}

    # Ensembl mammal genomes
    info = fetch_json('https://rest.ensembl.org/info/data?content-type=application/json', RAW / 'ensembl_release.json')
    genomes = fetch_json('https://rest.ensembl.org/info/genomes/taxonomy/Mammalia?content-type=application/json',
                         RAW / 'ensembl_mammal_genomes.json')
    ens = pd.DataFrame([{'ensembl_name': g['name'], 'display_name': g['display_name'], 'taxon_id': int(g['taxonomy_id']),
                         'assembly': g['assembly_name'], 'assembly_level': g.get('assembly_level'),
                         'peptide_compara': g.get('has_peptide_compara'), 'strain': g.get('strain')} for g in genomes])
    ens['binomial'] = ens.ensembl_name.map(binomial)
    ens['exact_binomial'] = ens.ensembl_name.str.lower() == ens.binomial
    # One genome per binomial: exact name if present, otherwise first alphabetically (fixed in the plan).
    ens = ens.sort_values(['binomial', 'exact_binomial', 'ensembl_name'], ascending=[True, False, True])
    ens['chosen_for_binomial'] = ~ens.duplicated('binomial')
    pheno = an[ok].drop_duplicates('binomial').set_index('binomial')
    ens['anage_eligible'] = ens.binomial.isin(pheno.index)
    ens['tree_label'] = ens.binomial.map(tips)
    # Amendment E1: NCBI Taxonomy synonyms of the same taxon id, exact binomials only.
    ens['anage_binomial'] = ens.binomial.where(ens.anage_eligible)
    synonyms = []
    for i, row in ens[ens.chosen_for_binomial & (ens.peptide_compara == 1)
                      & (~ens.anage_eligible | ens.tree_label.isna())].iterrows():
        names = ncbi_names(row.taxon_id)
        alternatives = sorted({n.lower().replace(' ', '_') for n in names if len(n.split()) == 2} - {row.binomial})
        hit_anage = [a for a in alternatives if a in pheno.index]
        hit_tree = [a for a in alternatives if a in tips]
        if not row.anage_eligible and len(hit_anage) == 1:
            ens.loc[i, 'anage_binomial'] = hit_anage[0]
        if pd.isna(row.tree_label) and len(hit_tree) == 1:
            ens.loc[i, 'tree_label'] = tips[hit_tree[0]]
        synonyms.append({'ensembl_name': row.ensembl_name, 'taxon_id': row.taxon_id, 'ncbi_binomials': ';'.join(alternatives),
                         'anage_match': ';'.join(hit_anage), 'tree_match': ';'.join(hit_tree)})
    pd.DataFrame(synonyms).to_csv(DATA / 'species_synonyms_ncbi.csv', index=False)
    ens['anage_eligible'] = ens.anage_binomial.notna()
    ens.to_csv(DATA / 'ensembl_species_audit.csv', index=False)

    sp = ens[ens.chosen_for_binomial & ens.anage_eligible & ens.tree_label.notna() & (ens.peptide_compara == 1)].copy()
    # Two Ensembl genomes must not share one AnAge record or one tree tip.
    sp = sp[~sp.duplicated('anage_binomial', keep=False) & ~sp.duplicated('tree_label', keep=False)]
    sp = sp.join(pheno[['Order', 'Family', 'Maximum longevity (yrs)', 'Adult weight (g)', 'Data quality',
                        'Sample size', 'Specimen origin', 'HAGRID']], on='anage_binomial')
    sp = sp.rename(columns={'Order': 'order', 'Family': 'family', 'Maximum longevity (yrs)': 'max_lifespan_yr',
                            'Adult weight (g)': 'adult_mass_g', 'Data quality': 'data_quality',
                            'Sample size': 'sample_size', 'Specimen origin': 'specimen_origin', 'HAGRID': 'hagrid'})
    non_placental = {'Monotremata', 'Didelphimorphia', 'Dasyuromorphia', 'Diprotodontia', 'Peramelemorphia',
                     'Microbiotheria', 'Paucituberculata', 'Notoryctemorphia'}
    sp['placental'] = ~sp.order.isin(non_placental)
    sp['LQ'] = sp.max_lifespan_yr / (4.88 * sp.adult_mass_g ** 0.153)
    sp.to_csv(DATA / 'species_candidates.csv', index=False)
    (DATA / 'ensembl_release.json').write_text(json.dumps(info))
    print('Ensembl release', info.get('releases'), '| mammal genomes', len(ens), '| unique binomials',
          ens.chosen_for_binomial.sum(), '| matched AnAge+tree+compara', len(sp),
          '| placental', int(sp.placental.sum()))
    print(sp.groupby('order').size().sort_values(ascending=False).to_string())


if __name__ == '__main__':
    main()
