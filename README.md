
<!-- README.md is generated from README.Rmd. Please edit that file -->

ROSeq - A rank based approach to modeling gene expression Author:
Krishan Gupta

date: 2/12/2019

## Introduction

ROSeq - A rank based approach to modeling gene expression with filtered
and normalized read count matrix. Takes in the complete filtered and
normalized read count matrix, the location of the two sub-populations
and the number of cores to be used.

## Installation

The developer version of the R package can be installed with the
following R commands:

library(devtools) install\_github(“krishan57gupta/ROSeq”)

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
#> g1   27756   27856   14996   15170   17536   17596   10680   25327   11115
#> g2  426294  263846  179803  205621  263160  257101  245900  295316  354956
#> g3   12066   23135    8180    4403   19442    4094    7123   16410   54821
#> g4   21052   15864    9887   30811   18720    7015   21955   21186   14550
#> g5  259737  335046  322967  258223  615820  387975  430683  423946  253189
#> g6   79393   25852   28569      12    7582     110   22632    3871   16691
#>    sample10 sample11 sample12 sample13 sample14 sample15 sample16
#> g1    24807    15339    29194    17566    16155    21736    10238
#> g2   503347   123557   340325   255690   307608   245292   202565
#> g3    32799    35674    68551    46630    20091    54231    35089
#> g4    25444    10804    19855    15379    29884     5807    16013
#> g5   483363   335090   503896   387804   385752   291252   407910
#> g6      251    21448    34916     6030   205805    18976        7
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
#>              pVals        pAdj
#>  [1,] 0.8810070374 1.000000000
#>  [2,] 1.0000000000 1.000000000
#>  [3,] 1.0000000000 1.000000000
#>  [4,] 1.0000000000 1.000000000
#>  [5,] 0.8810070374 1.000000000
#>  [6,] 1.0000000000 1.000000000
#>  [7,] 0.0783679611 0.208981230
#>  [8,] 0.0783679611 0.208981230
#>  [9,] 0.0110962445 0.032877761
#> [10,] 0.8810070374 1.000000000
#> [11,] 0.0110962445 0.032877761
#> [12,] 0.0110962445 0.032877761
#> [13,] 1.0000000000 1.000000000
#> [14,] 1.0000000000 1.000000000
#> [15,] 1.0000000000 1.000000000
#> [16,] 1.0000000000 1.000000000
#> [17,] 1.0000000000 1.000000000
#> [18,] 0.0783679611 0.208981230
#> [19,] 0.8810070374 1.000000000
#> [20,] 1.0000000000 1.000000000
#> [21,] 1.0000000000 1.000000000
#> [22,] 0.0110962445 0.032877761
#> [23,] 0.8810070374 1.000000000
#> [24,] 1.0000000000 1.000000000
#> [25,] 0.8810070374 1.000000000
#> [26,] 0.8810070374 1.000000000
#> [27,] 1.0000000000 1.000000000
#> [28,] 0.0110962445 0.032877761
#> [29,] 0.0110962445 0.032877761
#> [30,] 0.0110962445 0.032877761
#> [31,] 0.4950205039 1.000000000
#> [32,] 0.0110962445 0.032877761
#> [33,] 0.0110962445 0.032877761
#> [34,] 0.0110962445 0.032877761
#> [35,] 0.0110962445 0.032877761
#> [36,] 0.8810070374 1.000000000
#> [37,] 0.0110962445 0.032877761
#> [38,] 0.4172948490 1.000000000
#> [39,] 1.0000000000 1.000000000
#> [40,] 0.0110962445 0.032877761
#> [41,] 1.0000000000 1.000000000
#> [42,] 1.0000000000 1.000000000
#> [43,] 1.0000000000 1.000000000
#> [44,] 0.8810070374 1.000000000
#> [45,] 0.0110962445 0.032877761
#> [46,] 0.0110962445 0.032877761
#> [47,] 1.0000000000 1.000000000
#> [48,] 0.8810070374 1.000000000
#> [49,] 0.0110962445 0.032877761
#> [50,] 0.8810070374 1.000000000
#> [51,] 0.8810070374 1.000000000
#> [52,] 1.0000000000 1.000000000
#> [53,] 0.0005903587 0.009445739
#> [54,] 0.0110962445 0.032877761
#> [55,] 1.0000000000 1.000000000
#> [56,] 0.0005853661 0.009445739
#> [57,] 0.8810070374 1.000000000
#> [58,] 1.0000000000 1.000000000
#> [59,] 0.0110962445 0.032877761
#> [60,] 1.0000000000 1.000000000
#> [61,] 0.0005853661 0.009445739
#> [62,] 0.0110962445 0.032877761
#> [63,] 0.8810070374 1.000000000
#> [64,] 1.0000000000 1.000000000
#> [65,] 0.0005903587 0.009445739
#> [66,] 1.0000000000 1.000000000
#> [67,] 1.0000000000 1.000000000
#> [68,] 0.0110962445 0.032877761
#> [69,] 0.0005853661 0.009445739
#> [70,] 1.0000000000 1.000000000
#> [71,] 1.0000000000 1.000000000
#> [72,] 0.8810070374 1.000000000
#> [73,] 1.0000000000 1.000000000
#> [74,] 0.0110962445 0.032877761
#> [75,] 1.0000000000 1.000000000
#> [76,] 0.0110962445 0.032877761
#> [77,] 1.0000000000 1.000000000
#> [78,] 0.4824636857 1.000000000
#> [79,] 0.8810070374 1.000000000
#> [80,] 1.0000000000 1.000000000
```
