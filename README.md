# BayesSUR-RE

This repository is for the results in the paper *Structured Bayesian variable selection for multiple correlated response variables and high-dimensional predictors* by 'Zhao Z., Banterle M., Lewin A. and Zucknick M. (2021)', [arXiv:2101.05899](https://arxiv.org/abs/2101.05899).

# Data
## Abstract

The main data is the Genomics of Drug Sensitivity in Cancer (Garnett et al., 2012, Nature) which is publicly available. The gene set of MAPK is also is publicly available from Kyoto Encyclopedia of Genes and Genomes (KEGG).

## Availability

No restrictions.

## Description

The pharmacogenomic datasets in Garnett et al. (2012, Nature) can be downloaded from 
```diff
ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-5.0/
```
There are three data files needed: `gdsc_en_input_w5.csv`, `gdsc_drug_sensitivity_fitted_data_w5.csv` and `gdsc_tissue_output_w5.csv`.

The gene set of MAPK can be downloaded from KEGG gene set enrichment analysis
```diff
http://software.broadinstitute.org/gsea/msigdb/cards/KEGG_MAPK_SIGNALING_PATHWAY
```
Download the `.txt` version and delete the head content, i.e. the first two lines.

# Code
## Abstract

We are including all of the code that will enable reproducing our results.

## Description
We developed an R package BayesSUR which is available on CRAN for implementing our approach. All of the code are also available on the [first author’s GitHub](https://github.com/zhizuio/BayesSUR-RE/).

## Optional Information 
The real data analysis is computationally intensive due to high-dimensional genomic predictors. Running with one thread of CPU takes a few hours, but producing the same results in this article. It can be run faster on a cluster with parallelization. However, since the core code of our approach is in C++ for computational efficiency, Rcpp with OpenMP for parallelization is difficult reproduce the results on different type of machines.

# Instructions for Use
## Reproducibility

All tables, Figures 3-5, and Figures 7-10 can be reproduced through the provided code. The general steps are:

1. Load simulation functions through files `sim.ssur.R` and `sim.ssur.re.R`.
2. Run script `simulation_results.R` line by line to reproduce all simulation results in the article.
3. Download real data by following the above Data section.
4. Run script `GDSC_preprocess1.R` to get the real data ready for modelling (need to uncomment some lines to obtain datasets corresponding 'Feature sets II and III' in the article); run script `GDSC_preprocess2.R` to get the independent validation data.
5. Run script `GDSC_results.R` line by line to reproduce all results of real data analysis in the article.
