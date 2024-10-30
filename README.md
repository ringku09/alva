
<!-- README.md is generated from README.Rmd. Please edit that file -->

# toxassay <a href="https://dplyr.tidyverse.org"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/dplyr)](https://cran.r-project.org/package=dplyr)
[![R-CMD-check](https://github.com/tidyverse/dplyr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tidyverse/dplyr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/tidyverse/dplyr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/tidyverse/dplyr?branch=main)
<!-- badges: end -->

## Overview

**ToxAssay** is an R-based software package designed for the
comprehensive evaluation of drug-induced toxicity utilizing complex
toxicogenomics databases. The package offers a comprehensive suite of
functions, encompassing:

- Automated acquisition, preprocessing, and annotation of data by
  compound name for efficient management.
- Identification of molecular markers for targeted toxicity, such as
  differentially expressed genes (DEGs),adverse outcome pathways (AOPs),
  functional pathways, and Protein-Protein Interaction (PPI) networks,
  etc.
- Development of optimized machine-learning classifiers for predicting
  the targeted toxicity in test samples of compounds.

Workflow of ToxAssay in five main steps:
<img src="man/figures/workflow.png" width="890"/>

## Installation

``` r
# Install released version from CRAN
install.packages("toxassay")
```

To get a bug fix or to use a feature from the development version, you
can install the development version of toxassay from GitHub.

``` r
# install.packages("pak")
pak::pak("FanLiuLab/toxassay")
```

## Usages

### Download and preprocess analysis-ready data

In the current version of **ToxAssay**, users can download and process
raw compound perturbation data from the [Open
TG-GATEs](https://dbarchive.biosciencedbc.jp/en/open-tggates/download.html)
and [DrugMatrix](https://ntp.niehs.nih.gov/data/drugmatrix) database s,
as well as curated relational data from the [CTD](https://ctdbase.org).

Available compounds in databases for different parameters:

``` r
# Compounds in Open TG-GATEs data base
# Data available for compounds in Rat in-vivo liver Single data
tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")

# Compounds in Open DrugMatrix database
# Data available for compounds in `Liver`
drugmatrix_compounds(tissue = "Liver")
```

Automatically download and process data by compound name from the
databases:

``` r
# Get data from Open TG-GATEs
get_tggates(compounds = c("2-nitrofluorene", "N-nitrosomorpholine"), species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")

# Get data from Open DrugMatrix
get_drugmatrix(compounds = c("Ethanol", "Estriol"), tissue = "Liver")

# Get data from CTD
get_ctd(compounds = c("Ethanol", "Estriol"), gene_parm = "cgixns", disease_parm = "diseases")
```

<img src="man/figures/downloading.png" width="890"/>

### Differential expression analysis

``` r
library(toxassay)
set.seed(100)
sim_data <- simulate_tgxdata(n_de = 10, n_ee = 100, n_com = c(5,5), d=1)
com_group <- list(non_toxic = paste0("Compound", 1:5), toxic = paste0("Compound", 6:10))
tgx_degs(com_group, 
         ge_matrix = sim_data$expression, 
         metadata = sim_data$metadata, 
         gr_diff = TRUE, log10p = TRUE)
#> # A tibble: 110 × 8
#>   probe_id gene_symbol non_toxic toxic `non_toxic-toxic`  p_value log10p
#> * <chr>    <chr>           <dbl> <dbl>             <dbl>    <dbl>  <dbl>
#> 1 DE5      DE5            -0.507 0.541            -1.05  2.21e-35   34.7
#> 2 DE1      DE1            -0.423 0.673            -1.10  2.24e-34   33.6
#> 3 DE3      DE3            -0.509 0.415            -0.924 5.92e-30   29.2
#> 4 DE6      DE6            -0.494 0.557            -1.05  7.80e-30   29.1
#> 5 DE4      DE4            -0.800 0.293            -1.09  2.18e-29   28.7
#> # ℹ 105 more rows
#> # ℹ 1 more variable: sig_type <chr>
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/ringku09/toxassay/issues). For questions and
other discussion, please send an [E-mail](mailto:ringku_740@yahoo.com).
