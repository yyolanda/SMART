
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SMART

<!-- badges: start -->
<!-- badges: end -->

SMART performs deconvolution on spatial transcriptomics data using
marker-gene-assisted topic models.

## Installation

You can install the development version of SMART from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("yyolanda/SMART")
```

## Example

The base model of SMART takes the spatial transcriptomics data matrix
(`stMat`) and a list of marker gene symbols for each cell type
(`markers`) as inputs. Run the base model with following:

``` r
library(SMART)
data(MPOA)
res <- SMART_base(stData=MPOA$stMat, markerGs=MPOA$markers, 
                  noMarkerCts=1,outDir='SMART_results')
```
