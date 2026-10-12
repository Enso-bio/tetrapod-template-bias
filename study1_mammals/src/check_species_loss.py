"""Audit helper: why Ensembl mammal binomials did not enter the species set."""
import pandas as pd
from common import DATA

ens = pd.read_csv(DATA / 'ensembl_species_audit.csv')
an = pd.read_csv(DATA / 'anage_all.csv')
anb = set(an.binomial)
ch = ens[ens.chosen_for_binomial]
lost = ch[~(ch.anage_eligible & ch.tree_label.notna() & (ch.peptide_compara == 1))].copy()
lost['in_anage_any'] = lost.binomial.isin(anb)
q = an.drop_duplicates('binomial').set_index('binomial')
lost['anage_quality'] = lost.binomial.map(q['Data quality'])
lost['anage_sample'] = lost.binomial.map(q['Sample size'])
pd.set_option('display.width', 200)
print(lost[['ensembl_name', 'display_name', 'peptide_compara', 'in_anage_any', 'anage_quality', 'anage_sample',
            'anage_eligible', 'tree_label']].to_string(index=False))
