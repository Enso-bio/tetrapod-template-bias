# Study 2: reptiles and amphibians

40 reptiles and amphibians (turtles, crocodilians, squamates and amphibians) and 1,910 human genes (GenAge genes and
a random background), with zebrafish, spotted gar, coelacanth and human as reference genomes. Proteomes and gene
annotations from NCBI RefSeq (NCBI Datasets), one-to-one orthologs from the NCBI orthologs file, divergence times
from TimeTree 5, lifespan from AnAge and missing body masses from AmphiBIO or ReptTraits.

The design follows Study 1 with the lessons learned there: multiple alignments (MAFFT), a template outside the clade
analysed (zebrafish AlphaFold models) and the tetrapod ancestor as reference.

## Analysis plan

`PLAN_ANALISIS.md` (in Spanish) was written before protein, ortholog or structure data were retrieved;
`plan_hash.json` records its SHA256 at each stage. Dated amendments:

- E1: operational (ReptTraits version and temperature variables).
- E2: post hoc follow-up of the human-template association in reptiles and calibration of the per-gene tests.

## Scripts (`src/`, run in order by `run_all.py`)

| Script | Step |
|---|---|
| `00_get_inputs.py` | AnAge, GenAge, AmphiBIO and ReptTraits files |
| `01_species_phylo.py` | NCBI genomes, phenotypes, body mass and TimeTree-based tree |
| `02_download_genomes.py` | Protein FASTA and GFF3 files (NCBI Datasets) and NCBI orthologs |
| `03_genes_orthologs.py` | Gene universe, one-to-one orthologs and zebrafish AlphaFold templates |
| `04_msa_asr.py` | MAFFT alignments, coverage filters and ancestral reconstruction (IQ-TREE) |
| `05_thermompnn.py` | ThermoMPNN site energies on the zebrafish template |
| `06_metrics.py` | Species × gene metrics and predictor controls |
| `07_stats.py` | Primary tests, secondary and sensitivity analyses |
| `08_template_sensitivity.py` | Coelacanth and human templates |
| `11_posthoc_human_template.py`, `12_posthoc_pergene_calibration.py` | Amendment E2 |
| `09_figures.py`, `10_paper_tables.py`, `13_posthoc_tables.py`, `14_posthoc_figure.py` | Figures and manuscript macros |

`00_feasibility.py` and `00b_feasibility_ncbi.py` are the feasibility checks run before the plan was written.
Helpers: `common.py`, `pgls.py`, `thermo.py`. Check: `check_asr_nodes.py`.

## Other files

- `results/`: all result tables and JSON summaries.
- `data/`: species and gene tables and `provenance.jsonl` (URL, UTC retrieval date and SHA256 of every download).
  TimeTree divergence times are not redistributed; `01_species_phylo.py` queries them.
- `requirements.txt`: Python packages (same environment as Study 1, plus openpyxl).

Tools (IQ-TREE 2.4.0, MAFFT 7.526, ThermoMPNN) go in `tools/`; see the README at the repository root.
