
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

## SMART-base example

The base model of SMART takes the spatial transcriptomics data matrix
(`stMat`) and a list of marker gene symbols for each cell type
(`markers`) as inputs. 
Run the base model with following:

``` r
library(SMART)
data(visium_mouse_brain)
res <- SMART_base(stData=visium_mouse_brain$stMat,
                  markerGs=visium_mouse_brain$markers, 
                  noMarkerCts=1,outDir='SMART_results',seed=8001)
```

![alt text](https://github.com/yyolanda/SMART/blob/main/figures/visium_mousebrain_scatterpie.png?raw=true)

## SMART-covariate example

The covariate model of SMART takes the spatial transcriptomics data matrix
(`stMat`), a list of marker gene symbols for each cell type
(`markers`), and a data.frame of covariate levels (`covDat`) as inputs. 
Run the base model with following:

``` r
library(SMART)
data(MPOA)
mod <- SMART_covariate(stData=MPOA$stMat, markerGs=MPOA$markers,       
                       covarDat=MPOA$covDat,covars='sex',
                       noMarkerCts=1,outDir='SMART_results', seed=1)
res <- get_est(model=mod$model, stData=MPOA$stMat, 
               covarDat=MPOA$covDat, covar='sex')
```


