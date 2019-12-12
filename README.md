
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
n_samples<-8
n_genes<-80
n_de_genes<-40
mydata <- generateSyntheticData(dataset = "mydata", n.vars = n_genes, samples.per.cond = n_samples, n.diffexp = n_de_genes)
samples<-list()
samples$count<-mydata@count.matrix
samples$group<-mydata@sample.annotations$condition
head(samples$count)
#>    sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9
#> g1      87      72     754     755    5345    1662    3985    5558    1440
#> g2    3509    3457    3176    2268    3779    1560     942     841    4555
#> g3    1963    2832    4233    3146    2047    4262     823    9706    1595
#> g4  243746  307664  237762  265680  183348  207577  190751  219311  291206
#> g5  813285 1171247  892096  830688  924787  574875  557931  801904 1620191
#> g6  361412  433473  280553  341094  253030  219771  136059  187882  228770
#>    sample10 sample11 sample12 sample13 sample14 sample15 sample16
#> g1     6282    25172     7026     4505     6071     7030     4391
#> g2     1420     2603     2274     2090     1242     2344     3280
#> g3     2584     1115     2706     7021     2191     3651     3879
#> g4   196232   221162   242308   272768   141316   149760   317407
#> g5  2119321  2295027  1382800  2186911   517134  1721600  2175866
#> g6   366267   471256   249531   561744   218908   383814   100271
```

``` r
samples$count=apply(samples$count,2,function(x) as.numeric(x))
gkeep <- apply(samples$count,1,function(x) sum(x>0)>5)
samples$count<-samples$count[gkeep,]
samples$count<-edgeR::cpm(samples$count)
```

``` r
output<-ROSeq(countData=samples$count, condition = samples$group, nbits=10, numCores=1)
```

``` r
output
#>              pVals         pAdj
#>  [1,] 1.109624e-02 0.0467210294
#>  [2,] 8.810070e-01 1.0000000000
#>  [3,] 3.463372e-01 1.0000000000
#>  [4,] 1.000000e+00 1.0000000000
#>  [5,] 8.810070e-01 1.0000000000
#>  [6,] 1.000000e+00 1.0000000000
#>  [7,] 1.000000e+00 1.0000000000
#>  [8,] 1.000000e+00 1.0000000000
#>  [9,] 1.109624e-02 0.0467210294
#> [10,] 1.000000e+00 1.0000000000
#> [11,] 3.612172e-01 1.0000000000
#> [12,] 1.109624e-02 0.0467210294
#> [13,] 1.000000e+00 1.0000000000
#> [14,] 1.109624e-02 0.0467210294
#> [15,] 1.000000e+00 1.0000000000
#> [16,] 1.109624e-02 0.0467210294
#> [17,] 8.810070e-01 1.0000000000
#> [18,] 1.109624e-02 0.0467210294
#> [19,] 1.000000e+00 1.0000000000
#> [20,] 5.853661e-04 0.0078048819
#> [21,] 5.853661e-04 0.0078048819
#> [22,] 1.000000e+00 1.0000000000
#> [23,] 1.000000e+00 1.0000000000
#> [24,] 1.000000e+00 1.0000000000
#> [25,] 1.109624e-02 0.0467210294
#> [26,] 1.000000e+00 1.0000000000
#> [27,] 1.000000e+00 1.0000000000
#> [28,] 1.000000e+00 1.0000000000
#> [29,] 1.000000e+00 1.0000000000
#> [30,] 1.000000e+00 1.0000000000
#> [31,] 3.978029e-01 1.0000000000
#> [32,] 3.612172e-01 1.0000000000
#> [33,] 1.109624e-02 0.0467210294
#> [34,] 1.000000e+00 1.0000000000
#> [35,] 1.109624e-02 0.0467210294
#> [36,] 1.000000e+00 1.0000000000
#> [37,] 8.810070e-01 1.0000000000
#> [38,] 8.810070e-01 1.0000000000
#> [39,] 8.810070e-01 1.0000000000
#> [40,] 1.000000e+00 1.0000000000
#> [41,] 8.810070e-01 1.0000000000
#> [42,] 1.000000e+00 1.0000000000
#> [43,] 8.810070e-01 1.0000000000
#> [44,] 4.824637e-01 1.0000000000
#> [45,] 1.000000e+00 1.0000000000
#> [46,] 1.000000e+00 1.0000000000
#> [47,] 5.853661e-04 0.0078048819
#> [48,] 8.810070e-01 1.0000000000
#> [49,] 8.810070e-01 1.0000000000
#> [50,] 8.810070e-01 1.0000000000
#> [51,] 1.000000e+00 1.0000000000
#> [52,] 1.000000e+00 1.0000000000
#> [53,] 1.000000e+00 1.0000000000
#> [54,] 1.217270e-06 0.0000973816
#> [55,] 1.000000e+00 1.0000000000
#> [56,] 8.810070e-01 1.0000000000
#> [57,] 1.109624e-02 0.0467210294
#> [58,] 1.109624e-02 0.0467210294
#> [59,] 4.824637e-01 1.0000000000
#> [60,] 1.000000e+00 1.0000000000
#> [61,] 1.109624e-02 0.0467210294
#> [62,] 1.000000e+00 1.0000000000
#> [63,] 5.853661e-04 0.0078048819
#> [64,] 1.109624e-02 0.0467210294
#> [65,] 1.000000e+00 1.0000000000
#> [66,] 1.000000e+00 1.0000000000
#> [67,] 8.810070e-01 1.0000000000
#> [68,] 8.810070e-01 1.0000000000
#> [69,] 8.810070e-01 1.0000000000
#> [70,] 7.218966e-05 0.0028875863
#> [71,] 8.810070e-01 1.0000000000
#> [72,] 1.000000e+00 1.0000000000
#> [73,] 1.000000e+00 1.0000000000
#> [74,] 1.000000e+00 1.0000000000
#> [75,] 8.810070e-01 1.0000000000
#> [76,] 4.824637e-01 1.0000000000
#> [77,] 1.000000e+00 1.0000000000
#> [78,] 1.000000e+00 1.0000000000
#> [79,] 1.000000e+00 1.0000000000
#> [80,] 3.612172e-01 1.0000000000
```
