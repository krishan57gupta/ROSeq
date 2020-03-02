
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ROSeq - A rank based approach to modeling gene expression

Author: Krishan Gupta

## Introduction

ROSeq - A rank based approach to modeling gene expression with filtered
and normalized read count matrix. Takes in the complete filtered and
normalized read count matrix, the location of the two sub-populations
and the number of cores to be used.

## Installation

The developer version of the R package can be installed with the
following R commands:

``` r
if (\!requireNamespace(“BiocManager”, quietly = TRUE))
install.packages(“BiocManager”) BiocManager::install(‘ROSeq’)
```
or can be installed with the following R commands:
``` r
library(devtools) 
install\_github(‘krishan57gupta/ROSeq’)
```
## Vignette tutorial

This vignette uses a tung dataset already inbuilt in same package, to
demonstrate a standard pipeline. This vignette can be used as a tutorial
as well. Ref: Tung, P.-Y.et al.Batch effects and the effective design of
single-cell geneexpression studies.Scientific reports7, 39921 (2017).
\#\# Example

Libraries need to be loaded before running.

``` r
library(ROSeq)
library(compcodeR)
#> Loading required package: sm
#> Warning in fun(libname, pkgname): couldn't connect to display ":0"
#> Package 'sm', version 2.2-5.6: type help(sm) for summary information
library(edgeR)
#> Loading required package: limma
```

``` r
samples<-list()
samples$count<-ROSeq::L_Tung_single$NA19098_NA19101_count
samples$group<-ROSeq::L_Tung_single$NA19098_NA19101_group
head(samples$count[1:5,1:5])
#>                 NA19098.r1.A01 NA19098.r1.A02 NA19098.r1.A03 NA19098.r1.A04
#> ENSG00000237683              0              0              0              1
#> ENSG00000187634              0              0              0              0
#> ENSG00000188976              3              6              1              3
#> ENSG00000187961              0              0              0              0
#> ENSG00000187583              0              0              0              0
#>                 NA19098.r1.A05
#> ENSG00000237683              0
#> ENSG00000187634              0
#> ENSG00000188976              4
#> ENSG00000187961              0
#> ENSG00000187583              0
```

``` r
samples$count=apply(samples$count,2,function(x) as.numeric(x))
gkeep <- apply(samples$count,1,function(x) sum(x>0)>5)
samples$count<-samples$count[gkeep,]
samples$count<-edgeR::cpm(samples$count)
```

``` r
output<-ROSeq(countData=samples$count, condition = samples$group, nbits=0, numCores=1)
```

``` r
output[1:5,]
#>           pVals      pAdj       log2FC
#> [1,] 0.75009454 0.8766822  0.003179705
#> [2,] 0.05000403 0.1732159 -1.882145352
#> [3,] 0.08819975 0.2440436  0.043243212
#> [4,] 0.56798417 0.7405707 -0.622076958
#> [5,] 0.34364800 0.5498973  0.602577144
```
