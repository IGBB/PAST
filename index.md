---
layout: page
permalink: /
title: PAST - Pathway Association Study Tool
---

The Pathway Association Study Tool (PAST) was developed to facilitate easier and more efficient GWAS-based metabolic pathway analysis. PAST was designed for use with maize but is usable for other species as well. It tracks all SNP marker - trait associations, regardless of significance or magnitude. PAST groups SNPs into linkage blocks based on linkage disequilibrium (LD) data and identifies a tagSNP from each block. PAST then identifies genes within a user-defined distance of the tagSNPs, and transfers the attributes of the tagSNP to the gene(s), including the allele effect, R2 and p-value of the original SNP-trait association found from the GWAS analysis.  Finally, PAST uses the gene effect values to calculate an enrichment score (ES) and p-value for each pathway. PAST is easy to use as an online tool, standalone R script, or as a downloadable R Shiny application. It uses as input TASSEL files that are generated as output from the General Linear or Mixed Linear Models (GLM and MLM), or files from any association analysis that has been similarly formatted, as well as genome annotations in GFF format, and a metabolic pathways file.

## Installation

PAST requires R > 3.5. R can be downloaded from [here](https://www.r-project.org/). R can often be installed via package manager on Linux. Users who are unfamiliar with R are encouraged to use [RStudio](https://rstudio.com/products/rstudio/) to run the following commands. Once R is installed, PAST can be installed from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/PAST.html) with the following lines of code in an R Console.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PAST")
```

## Example Data

Example data is included in the PAST package (see [Running PAST through R](running_r)). The data can also be downloaded [here](/assets/example_data.zip). This example data only contains two pathways and is meant to illustrate how PAST works.

## Citation

Thrash A, Tang JD, DeOrnellis M, Peterson DG, Warburton ML (2020). “PAST: The Pathway Association Studies Tool to Infer Biological Meaning from GWAS Datasets.” Plants, 9(1), 58. [https://doi.org/10.3390/plants9010058](https://doi.org/10.3390/plants9010058). 