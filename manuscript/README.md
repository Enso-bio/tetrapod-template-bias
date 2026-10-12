# Manuscript

LaTeX source of the manuscript. Numbers, tables and figures are not typed by hand: they are generated from the result
files of the two studies.

- `build_inputs.py` collects the figures and number macros produced by each study (`study*/paper/figures` and
  `study*/paper/generated`, written by `09_figures.py` and the `*_paper_tables.py` scripts), prefixes the macros
  (`M` = mammals, `R` = reptiles and amphibians), computes the percentage and power macros and calls
  `build_methods_tables.py`.
- `build_methods_tables.py` builds the Supplementary Methods tables (software versions, data sources with retrieval
  dates, genome assemblies) from the studies' version records, provenance logs and species tables.
- `collect_methods_facts.py` prints the same facts for checking.
- `paper/` contains `main.tex`, the bibliography and the generated inputs as used for the submitted version, so the
  manuscript compiles without re-running the pipeline:

```bash
cd paper
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```
