
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
TMMnormalization <- function(countTable){
  nf=calcNormFactors(countTable ,method= "TMM")
  nf= colSums(countTable)*nf
  scalingFactors = nf/mean(nf)
  countTableTMM <- t(t(countTable)/scalingFactors)
  return(countTableTMM)
}
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
#> g1  887032 1677034 1372273  332156  803457  607740 1142627 1904313 1819862
#> g2   83559  200218  123773   51849   64747  113590   77668  126622   79677
#> g3     982     662    1031     237     955     811     625     741     879
#> g4   11276   18398    2602    3139    1151    5972    5603    3559    2575
#> g5  144956  239432  409210  159390  236250  270163  238493  186287  340569
#> g6      10     341    1159      59       0     956    1051     340     144
#>    sample10 sample11 sample12 sample13 sample14 sample15 sample16
#> g1  1551774   861660  2545325  1153328  1282872  2421862  2212632
#> g2    86632    37582    88417   103230    84350    83530   113202
#> g3      677      460     1112     2059      396     1080      585
#> g4     3564     3500     1094     3265     1746     1192     4046
#> g5   183528   165918   437415   171600   109330   144774   229697
#> g6      168        0      933        1       72     3688       80
```

``` r
samples$count=apply(samples$count,2,function(x) as.numeric(x))
gkeep <- apply(samples$count,1,function(x) sum(x>0)>5)
samples$count=samples$count[gkeep,]
samples$count=TMMnormalization(samples$count)
```

``` r
output=ROSeq(countData=samples$count, condition = samples$group, numCores=1)
```

``` r
output
#>              pVals         pAdj
#>  [1,] 1.109624e-02 0.0467210294
#>  [2,] 1.000000e+00 1.0000000000
#>  [3,] 1.217270e-06 0.0000486908
#>  [4,] 1.109624e-02 0.0467210294
#>  [5,] 1.000000e+00 1.0000000000
#>  [6,] 8.970041e-01 1.0000000000
#>  [7,] 1.200996e-03 0.0192159283
#>  [8,] 1.000000e+00 1.0000000000
#>  [9,] 1.109624e-02 0.0467210294
#> [10,] 1.000000e+00 1.0000000000
#> [11,] 1.000000e+00 1.0000000000
#> [12,] 5.415617e-01 1.0000000000
#> [13,] 1.000000e+00 1.0000000000
#> [14,] 1.109624e-02 0.0467210294
#> [15,] 1.000000e+00 1.0000000000
#> [16,] 1.000000e+00 1.0000000000
#> [17,] 4.172948e-01 1.0000000000
#> [18,] 1.109624e-02 0.0467210294
#> [19,] 1.109624e-02 0.0467210294
#> [20,] 1.000000e+00 1.0000000000
#> [21,] 5.853661e-04 0.0117073229
#> [22,] 1.000000e+00 1.0000000000
#> [23,] 1.000000e+00 1.0000000000
#> [24,] 1.000000e+00 1.0000000000
#> [25,] 8.810070e-01 1.0000000000
#> [26,] 1.000000e+00 1.0000000000
#> [27,] 1.000000e+00 1.0000000000
#> [28,] 1.109624e-02 0.0467210294
#> [29,] 8.810070e-01 1.0000000000
#> [30,] 1.109624e-02 0.0467210294
#> [31,] 8.810070e-01 1.0000000000
#> [32,] 1.000000e+00 1.0000000000
#> [33,] 1.000000e+00 1.0000000000
#> [34,] 1.217270e-06 0.0000486908
#> [35,] 8.810070e-01 1.0000000000
#> [36,] 1.000000e+00 1.0000000000
#> [37,] 1.000000e+00 1.0000000000
#> [38,] 1.000000e+00 1.0000000000
#> [39,] 1.000000e+00 1.0000000000
#> [40,] 8.810070e-01 1.0000000000
#> [41,] 8.810070e-01 1.0000000000
#> [42,] 1.109624e-02 0.0467210294
#> [43,] 3.612172e-01 1.0000000000
#> [44,] 1.000000e+00 1.0000000000
#> [45,] 8.810070e-01 1.0000000000
#> [46,] 1.000000e+00 1.0000000000
#> [47,] 1.109624e-02 0.0467210294
#> [48,] 1.109624e-02 0.0467210294
#> [49,] 1.000000e+00 1.0000000000
#> [50,] 1.109624e-02 0.0467210294
#> [51,] 1.000000e+00 1.0000000000
#> [52,] 1.000000e+00 1.0000000000
#> [53,] 8.810070e-01 1.0000000000
#> [54,] 1.000000e+00 1.0000000000
#> [55,] 8.810070e-01 1.0000000000
#> [56,] 4.824637e-01 1.0000000000
#> [57,] 8.810070e-01 1.0000000000
#> [58,] 1.000000e+00 1.0000000000
#> [59,] 2.944461e-02 0.1177784357
#> [60,] 1.000000e+00 1.0000000000
#> [61,] 1.000000e+00 1.0000000000
#> [62,] 1.000000e+00 1.0000000000
#> [63,] 7.836796e-02 0.2985446139
#> [64,] 4.950205e-01 1.0000000000
#> [65,] 1.000000e+00 1.0000000000
#> [66,] 8.810070e-01 1.0000000000
#> [67,] 1.000000e+00 1.0000000000
#> [68,] 1.000000e+00 1.0000000000
#> [69,] 1.000000e+00 1.0000000000
#> [70,] 8.810070e-01 1.0000000000
#> [71,] 1.000000e+00 1.0000000000
#> [72,] 8.810070e-01 1.0000000000
#> [73,] 3.978029e-01 1.0000000000
#> [74,] 4.824637e-01 1.0000000000
#> [75,] 1.109624e-02 0.0467210294
#> [76,] 1.109624e-02 0.0467210294
#> [77,] 1.000000e+00 1.0000000000
#> [78,] 1.000000e+00 1.0000000000
#> [79,] 5.853661e-04 0.0117073229
#> [80,] 1.000000e+00 1.0000000000
```
