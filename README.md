
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

install\_github(“krishan57gupta/ROSeq”)

## Vignette tutorial

This vignette uses a small simulated data, to demonstrate a standard
pipeline. This vignette can be used as a tutorial as well.

## Example

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
n_samples=8
n_genes=80
n_de_genes=40
mydata <- generateSyntheticData(dataset = "mydata", n.vars = n_genes, samples.per.cond = n_samples, n.diffexp = n_de_genes)
samples=list()
samples$count=mydata@count.matrix
samples$group=mydata@sample.annotations$condition
head(samples$count)
#>    sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9
#> g1 2056045 2290491 2718801 1333862 1886137 1390042 1632543 1942499 2004335
#> g2   73609  215422  187118  149236  169996  237892  140990  262417  129362
#> g3     381    1824    1219     650    1831    5364    2546     733    4302
#> g4  156480  123282   90208   91322  155018   47613   69752   81955  105999
#> g5   39736   23485   49667   18461   30440   21621   14239   21104   71991
#> g6   39912   56223   71500   35357   41198   80787   72443   68313   76313
#>    sample10 sample11 sample12 sample13 sample14 sample15 sample16
#> g1  2051525  2376852  3259187  1871437  1013470  2039306  1873336
#> g2   233315   101101   283907   341021   237257   500810   298621
#> g3     3355      840     1364     3909     2953     3733     8560
#> g4   109535    48966   217667    79984    64135   179715    63412
#> g5    52022    57039    94626   105057    50264    18529    70295
#> g6    47912    61340    74679    66145    61980    87384   100585
```

``` r
samples$count=apply(samples$count,2,function(x) as.numeric(x))
gkeep <- apply(samples$count,1,function(x) sum(x>0)>5)
samples$count=samples$count[gkeep,]
samples$count=edgeR::cpm(samples$count)
```

``` r
output=ROSeq(countData=samples$count, condition = samples$group, numCores=1)
```

``` r
output
#>              pVals         pAdj
#>  [1,] 8.810070e-01 1.0000000000
#>  [2,] 1.109624e-02 0.0355079824
#>  [3,] 1.109624e-02 0.0355079824
#>  [4,] 5.903587e-04 0.0118071738
#>  [5,] 5.903587e-04 0.0118071738
#>  [6,] 1.000000e+00 1.0000000000
#>  [7,] 1.000000e+00 1.0000000000
#>  [8,] 4.824637e-01 1.0000000000
#>  [9,] 1.000000e+00 1.0000000000
#> [10,] 1.109624e-02 0.0355079824
#> [11,] 1.000000e+00 1.0000000000
#> [12,] 1.000000e+00 1.0000000000
#> [13,] 8.810070e-01 1.0000000000
#> [14,] 1.109624e-02 0.0355079824
#> [15,] 1.109624e-02 0.0355079824
#> [16,] 1.000000e+00 1.0000000000
#> [17,] 4.824637e-01 1.0000000000
#> [18,] 3.978029e-01 1.0000000000
#> [19,] 1.000000e+00 1.0000000000
#> [20,] 1.000000e+00 1.0000000000
#> [21,] 8.810070e-01 1.0000000000
#> [22,] 4.172948e-01 1.0000000000
#> [23,] 1.000000e+00 1.0000000000
#> [24,] 1.000000e+00 1.0000000000
#> [25,] 1.000000e+00 1.0000000000
#> [26,] 1.000000e+00 1.0000000000
#> [27,] 1.000000e+00 1.0000000000
#> [28,] 1.000000e+00 1.0000000000
#> [29,] 1.000000e+00 1.0000000000
#> [30,] 1.000000e+00 1.0000000000
#> [31,] 1.000000e+00 1.0000000000
#> [32,] 1.000000e+00 1.0000000000
#> [33,] 1.000000e+00 1.0000000000
#> [34,] 1.109624e-02 0.0355079824
#> [35,] 1.000000e+00 1.0000000000
#> [36,] 1.000000e+00 1.0000000000
#> [37,] 1.109624e-02 0.0355079824
#> [38,] 5.915839e-03 0.0355079824
#> [39,] 1.109624e-02 0.0355079824
#> [40,] 1.109624e-02 0.0355079824
#> [41,] 7.836796e-02 0.2161874790
#> [42,] 1.109624e-02 0.0355079824
#> [43,] 1.000000e+00 1.0000000000
#> [44,] 7.836796e-02 0.2161874790
#> [45,] 8.810070e-01 1.0000000000
#> [46,] 7.836796e-02 0.2161874790
#> [47,] 2.944461e-02 0.0905987967
#> [48,] 1.000000e+00 1.0000000000
#> [49,] 1.109624e-02 0.0355079824
#> [50,] 8.810070e-01 1.0000000000
#> [51,] 1.000000e+00 1.0000000000
#> [52,] 1.000000e+00 1.0000000000
#> [53,] 1.109624e-02 0.0355079824
#> [54,] 1.109624e-02 0.0355079824
#> [55,] 1.000000e+00 1.0000000000
#> [56,] 1.000000e+00 1.0000000000
#> [57,] 1.109624e-02 0.0355079824
#> [58,] 1.000000e+00 1.0000000000
#> [59,] 1.109624e-02 0.0355079824
#> [60,] 1.000000e+00 1.0000000000
#> [61,] 1.000000e+00 1.0000000000
#> [62,] 1.217270e-06 0.0000973816
#> [63,] 1.109624e-02 0.0355079824
#> [64,] 1.000000e+00 1.0000000000
#> [65,] 1.000000e+00 1.0000000000
#> [66,] 1.109624e-02 0.0355079824
#> [67,] 1.000000e+00 1.0000000000
#> [68,] 1.000000e+00 1.0000000000
#> [69,] 8.810070e-01 1.0000000000
#> [70,] 1.109624e-02 0.0355079824
#> [71,] 1.000000e+00 1.0000000000
#> [72,] 1.000000e+00 1.0000000000
#> [73,] 1.000000e+00 1.0000000000
#> [74,] 1.000000e+00 1.0000000000
#> [75,] 8.810070e-01 1.0000000000
#> [76,] 5.903587e-04 0.0118071738
#> [77,] 1.000000e+00 1.0000000000
#> [78,] 1.000000e+00 1.0000000000
#> [79,] 1.109624e-02 0.0355079824
#> [80,] 1.109624e-02 0.0355079824
```
