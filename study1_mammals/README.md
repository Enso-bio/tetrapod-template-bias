# Study 1: placental mammals

57 placental mammals and 1,755 human genes (GenAge genes and a random background), with six marsupials and
monotremes as outgroups. Orthologs and proteomes from Ensembl 116 (Compara one-to-one homologies), phylogeny from
the node-dated MamPhy tree (Upham et al. 2019), lifespan and body mass from AnAge.

For each gene, the placental ancestor is reconstructed with IQ-TREE on the fixed species topology and threaded onto
an AlphaFold model; ThermoMPNN predicts the folding-stability effect (ΔΔG) of every amino-acid difference between the
ancestor and each species. Species-level metrics are related to maximum lifespan by phylogenetic generalized least
squares with Pagel's λ.

## Analysis plan

`PLAN_ANALISIS.md` (in Spanish) was written before any protein data were retrieved; `plan_hash.json` records its
SHA256 at each stage. Dated amendments:

- E1–E3: operational (species name reconciliation, background gene order, Ensembl FTP bulk files).
- E4: opossum template (post hoc).
- E5: platypus template (confirmatory, specified before any of its data were retrieved).
- E6: reanalysis with multiple alignments and sites independent of the human sequence (post hoc).

## Scripts (`src/`, run in order by `run_all.py`)

| Script | Step |
|---|---|
| `01_phenotype_tree.py` | AnAge phenotypes, MamPhy trees, Ensembl genome list and species matching |
| `02_genes.py` | GenAge genes and background gene order (UniProt) |
| `03a_download_ensembl_ftp.py` | Ensembl 116 Compara homologies and peptide files |
| `03_orthologs_structures.py` | One-to-one orthologs, alignment to human, AlphaFold DB models |
| `04_asr.py` | Coverage filters and ancestral reconstruction (IQ-TREE) |
| `05_thermompnn.py` | ThermoMPNN site energies on the human template |
| `06_metrics.py` | Species × gene metrics and predictor controls |
| `07_stats.py` | Primary tests (PGLS, bootstrap, competitive test), secondary and sensitivity analyses |
| `08_template_sensitivity.py` | Mouse template (pre-specified) |
| `10_outgroup_template.py` | Opossum (E4) and platypus (E5) templates |
| `12_confirm_template.py` | Confirmatory tests of amendment E5 |
| `13_msa_platypus.py` | Amendment E6 |
| `09_figures.py`, `11_paper_tables.py`, `14_posthoc_tables.py` | Figures and manuscript macros |

Helpers: `common.py` (paths and logged downloads), `pgls.py` (phylogenetic GLS), `thermo.py` (ThermoMPNN wrapper).
Checks: `check_asr_nodes.py`, `check_compara_human.py`, `check_species_loss.py`, `test_asr.py`, `test_thermo.py`.

## Other files

- `results/`: all result tables and JSON summaries.
- `data/`: species and gene tables and `provenance.jsonl` (URL, UTC retrieval date and SHA256 of every download).
- `requirements.txt` and `runtime_versions.json`: Python packages and exact versions.

Tools (IQ-TREE 2.4.0, MAFFT 7.526, ThermoMPNN) go in `tools/`; see the README at the repository root.
