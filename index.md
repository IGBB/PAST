# PAST: Pathway Association Study Tool

The Pathway Association Study Tool (PAST) was developed to facilitate easier and more efficient GWAS-based metabolic pathway analysis. PAST was designed for use with maize but is usable for other species as well. It tracks all SNP marker - trait associations, regardless of significance or magnitude. PAST groups SNPs into linkage blocks based on linkage disequilibrium (LD) data and identifies a tagSNP from each block. PAST then identifies genes within a user-defined distance of the tagSNPs, and transfers the attributes of the tagSNP to the gene(s), including the allele effect, R2 and p-value of the original SNP-trait association found from the GWAS analysis.  Finally, PAST uses the gene effect values to calculate an enrichment score (ES) and p-value for each pathway. PAST is easy to use as an online tool, standalone R script, or as a downloadable R Shiny application. It uses as input TASSEL files that are generated as output from the General Linear or Mixed Linear Models (GLM and MLM), or files from any association analysis that has been similarly formatted, as well as genome annotations in GFF format, and a metabolic pathways file.

# Installation

PAST is available on [Bioconductor](https://bioconductor.org/packages/release/bioc/html/PAST.html).

To install PAST from GitHub, run the following commands in R:

```R
install.packages("devtools")
library("devtools")
install_github("IGBB/PAST")
```

# Running PAST

## Functions Overview

* `load_GWAS_data`: takes TASSEL-formatted association and statistics files and the names of the necessary columns in the files as input and loads all of the GWAS data into a single dataframe.
* `load_LD`: takes a linkage disequilibrium file, cleans it up, and splits it by chromosome
* `assign_SNPs_to_genes`: takes loaded GWAS data, loaded LD data, a GFF annotations file, a window size, a r^2 cut-off value, and a number of cores and assigned SNPs to genes
* `find_pathway_significance`: takes genes, a pathways database file, a minimum number of genes that must be in a pathway to retain the pathway for analysis, the mode of the analysis, and a number of cores and calculates the enrichment score, the p-value, and the q-value for each pathway
* `plot_pathways`: takes the pathways data from `find_pathway_significance`, a filter type (p-value or q-value), a filter cut-off, the mode used to generate the data, and a directory and plots the pathways that pass the filter, saving them to the given directory

# Citation

Thrash A, Tang JD, DeOrnellis M, Peterson DG, Warburton ML (2020). “PAST: The Pathway Association Studies Tool to Infer Biological Meaning from GWAS Datasets.” Plants, 9(1), 58. [https://doi.org/10.3390/plants9010058](https://doi.org/10.3390/plants9010058). 