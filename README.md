# Structural templates, predicted protein stability and tetrapod longevity

Code, analysis plans and results for:

> Torres Alegre, E. O. and Mora Jimenez, D. E. *The choice of structural template can create apparent links between
> predicted protein stability and longevity in tetrapods.* (manuscript submitted)

Corresponding author: Enso Onill Torres Alegre (etorresa@nd.edu).

Two comparative studies relate the predicted folding-stability effect (ThermoMPNN ΔΔG) of the amino-acid differences
between a reconstructed ancestor and each species to maximum lifespan, using phylogenetic generalized least squares.
The analyses also show how the choice of structural template can create apparent associations, and how to detect it.

## Contents

| Folder | Contents |
|---|---|
| `study1_mammals/` | 57 placental mammals, 1,755 genes (Ensembl 116); human, mouse, opossum and platypus templates |
| `study2_reptiles_amphibians/` | 40 reptiles and amphibians, 1,910 genes (NCBI RefSeq); zebrafish, coelacanth and human templates |
| `manuscript/` | LaTeX source of the manuscript; every number and table is generated from the result files |

Each study contains:

- `PLAN_ANALISIS.md`: the analysis plan written before protein data were retrieved, with dated amendments (in
  Spanish); `plan_hash.json` records its SHA256 at each stage.
- `src/`: numbered scripts, run in order by `run_all.py`.
- `results/`: all result tables and summaries used in the manuscript.
- `data/`: small species and gene tables, and `provenance.jsonl`, which records the URL, UTC retrieval date and
  SHA256 of every downloaded file.

Large downloads (genomes, proteomes, AlphaFold models) and intermediate files (alignments, ancestral
reconstructions, ThermoMPNN energies) are not included; the scripts download and regenerate them.

## Requirements

- Python 3.12 with the packages in each study's `requirements.txt` (exact versions in
  `study1_mammals/runtime_versions.json`).
- External tools, placed in each study's `tools/` folder:
  - IQ-TREE 2.4.0: `tools/iqtree-2.4.0-Windows/bin/iqtree2.exe`
  - MAFFT 7.526: `tools/mafft-win/mafft.bat`
  - ThermoMPNN (Kuhlman lab, main branch, default weights `thermoMPNN_default.pt`): `tools/ThermoMPNN-main/`
- A CUDA GPU speeds up ThermoMPNN but is not required.
- pdflatex and bibtex for the manuscript.

The pipeline was developed and run on Windows 11. On other systems, adjust the tool paths in `src/common.py`,
`src/thermo.py` and the alignment and reconstruction scripts.

## Reproducing the analyses

```bash
cd study1_mammals
python run_all.py
cd ../study2_reptiles_amphibians
python run_all.py
cd ../manuscript
python build_inputs.py
cd paper
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Downloads are cached in `data/raw`, so steps can be re-run individually (for example `python run_all.py 06 07`).
Each study downloads several gigabytes and writes large intermediate files.

## Data sources and terms

AnAge and GenAge (Human Ageing Genomic Resources), Ensembl, NCBI RefSeq and Datasets, NCBI Taxonomy and gene
orthologs, UniProtKB, AlphaFold DB, MamPhy (Upham et al. 2019), AmphiBIO, ReptTraits and TimeTree 5. These data are
subject to their providers' terms. TimeTree divergence times are queried by the code and not redistributed here.

## Licence and citation

Code is released under the MIT licence (`LICENSE`). Please cite the article and this repository (`CITATION.cff`).
