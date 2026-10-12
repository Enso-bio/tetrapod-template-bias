"""Run the whole reptile/amphibian pipeline in order (downloads are cached in data/raw).

Usage: python run_all.py [step prefixes...]   e.g.  python run_all.py 06 07 09 10
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = ['00_get_inputs', '01_species_phylo', '02_download_genomes', '03_genes_orthologs', '04_msa_asr', '05_thermompnn', '06_metrics',
         '07_stats', '08_template_sensitivity', '09_figures', '10_paper_tables',
         '11_posthoc_human_template', '12_posthoc_pergene_calibration', '13_posthoc_tables', '14_posthoc_figure']  # E2

wanted = sys.argv[1:]
for step in STEPS:
    if wanted and not any(step.startswith(w) for w in wanted):
        continue
    print(f'=== {step}', flush=True)
    subprocess.run([sys.executable, str(ROOT / 'src' / f'{step}.py')], cwd=ROOT / 'src', check=True)
