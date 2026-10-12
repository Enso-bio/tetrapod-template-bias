"""Run the whole pipeline in order. Downloads are cached in data/raw, so a re-run reuses them.

Usage (from this folder, with the virtual environment active):
    python run_all.py            # every step
    python run_all.py 06 07 09   # only the listed steps
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = [('01_phenotype_tree', []), ('02_genes', []), ('03a_download_ensembl_ftp', []), ('03_orthologs_structures', []),
         ('04_asr', []), ('05_thermompnn', []), ('06_metrics', []), ('07_stats', []), ('08_template_sensitivity', []),
         ('10_outgroup_template', []),                                                   # opossum, post hoc (E4)
         ('10_outgroup_template', ['ornithorhynchus_anatinus', '9258', 'E5_template']),  # platypus (E5)
         ('12_confirm_template', []), ('09_figures', []), ('11_paper_tables', []),
         ('13_msa_platypus', []), ('14_posthoc_tables', [])]                             # post hoc (E6)

wanted = sys.argv[1:]
for step, args in STEPS:
    if wanted and not any(step.startswith(w) for w in wanted):
        continue
    print(f'=== {step} {" ".join(args)}', flush=True)
    subprocess.run([sys.executable, str(ROOT / 'src' / f'{step}.py'), *args], cwd=ROOT / 'src', check=True)
