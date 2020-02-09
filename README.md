# Hughes et al

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/j-berg/hughes_rnaseq_2019/master)

The following contains the necessary code for reproducing select figures from the associated manuscript   
The [interactive notebook](https://mybinder.org/v2/gh/j-berg/hughes_rnaseq_2019/master) contains the relevant code for reproducing the normalization and processing of the conA RNAseq time course data, as well as the code to reproduce figures relating to this data. This notebook is hosted on BinderHub, so it can be run interactively within your browser.   

### Associated citation:
```
Hughes CE, Coody TK, Jeong M, Berg JA, Winge DR, Hughes AL. Amino acid toxicity drives age-related 
mitochondrial decline by altering iron metabolism. (2020) Cell. DOI: 10.1016/j.cell.2019.12.035
```

### Navigation:
- docs: Contains markdown files for associated [website](https://j-berg.github.io/hughes_rnaseq_2019/)
- go_final: A deeper look into the genes enriched within oxidative reductive GO terms at the 6hr timepoint
- metadata: Metadata matrices and lists for plotting within the Jupyter notebook
- plots: Plots output during analysis (code contained in the Jupyter notebook)
  - volcano_plots: Volcano plots for each time-point
  - go: GO term enrichment plots
    - data: Raw threshold data from each time-point
- processed_data: Normalized count tables for conA RNA-seq samples used in study
- raw_data: Raw data files output by the University of Utah Sequencing Core and Bioinformatics core

