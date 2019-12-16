
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
#> g1  512133 1181751  315219 1188471  474266  550543  536917  874092 1451815
#> g2    6821    3202    4921    4111   17350    9786   61689   10509   31487
#> g3   76303   59987   54713   46357   63074  127375   45984   49261  102099
#> g4  167465  107886  122828   99240  110297  119792  142645  132061  156902
#> g5   21035   18048   13706   17034    7615   11723   16735   18567   64021
#> g6     177     332     341     826     296     158    1019     454     585
#>    sample10 sample11 sample12 sample13 sample14 sample15 sample16
#> g1   411762  2241466   406504   809686  1186339   266116  1184855
#> g2    16683      944     6289      125    48179     5912    19635
#> g3    47417    69272    52844    57943    57269    29546    40994
#> g4   122008   237716   178541   173678   240362   137913   240891
#> g5    17362    39367    27804    36173    62737    42718    47700
#> g6      172      347      249      273      470      368      559
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
output
#>              pVals         pAdj
#>  [1,] 2.988297e-02 7.711734e-02
#>  [2,] 1.000000e+00 1.000000e+00
#>  [3,] 1.000000e+00 1.000000e+00
#>  [4,] 1.000000e+00 1.000000e+00
#>  [5,] 1.000000e+00 1.000000e+00
#>  [6,] 1.000000e+00 1.000000e+00
#>  [7,] 2.988297e-02 7.711734e-02
#>  [8,] 7.540960e-04 4.309120e-03
#>  [9,] 4.637550e-02 8.627999e-02
#> [10,] 8.558477e-02 1.521507e-01
#> [11,] 1.000000e+00 1.000000e+00
#> [12,] 2.553986e-04 1.857444e-03
#> [13,] 4.637550e-02 8.627999e-02
#> [14,] 2.988297e-02 7.711734e-02
#> [15,] 7.540960e-04 4.309120e-03
#> [16,] 1.000000e+00 1.000000e+00
#> [17,] 5.399796e-05 4.849143e-04
#> [18,] 4.637550e-02 8.627999e-02
#> [19,] 1.000000e+00 1.000000e+00
#> [20,] 1.000000e+00 1.000000e+00
#> [21,] 5.455286e-05 4.849143e-04
#> [22,] 5.399796e-05 4.849143e-04
#> [23,] 1.000000e+00 1.000000e+00
#> [24,] 1.000000e+00 1.000000e+00
#> [25,] 2.988297e-02 7.711734e-02
#> [26,] 5.399796e-05 4.849143e-04
#> [27,] 1.000000e+00 1.000000e+00
#> [28,] 2.988297e-02 7.711734e-02
#> [29,] 1.000000e+00 1.000000e+00
#> [30,] 5.911492e-01 1.000000e+00
#> [31,] 4.637550e-02 8.627999e-02
#> [32,] 2.988297e-02 7.711734e-02
#> [33,] 2.988297e-02 7.711734e-02
#> [34,] 1.000000e+00 1.000000e+00
#> [35,] 1.000000e+00 1.000000e+00
#> [36,] 1.000000e+00 1.000000e+00
#> [37,] 1.000000e+00 1.000000e+00
#> [38,] 2.884585e-01 5.016670e-01
#> [39,] 1.000000e+00 1.000000e+00
#> [40,] 1.000000e+00 1.000000e+00
#> [41,] 2.988297e-02 7.711734e-02
#> [42,] 2.277580e-02 7.711734e-02
#> [43,] 1.000000e+00 1.000000e+00
#> [44,] 7.540960e-04 4.309120e-03
#> [45,] 1.000000e+00 1.000000e+00
#> [46,] 1.000000e+00 1.000000e+00
#> [47,] 2.553986e-04 1.857444e-03
#> [48,] 5.399796e-05 4.849143e-04
#> [49,] 3.591405e-08 2.873124e-06
#> [50,] 1.000000e+00 1.000000e+00
#> [51,] 4.637550e-02 8.627999e-02
#> [52,] 2.988297e-02 7.711734e-02
#> [53,] 1.000000e+00 1.000000e+00
#> [54,] 1.000000e+00 1.000000e+00
#> [55,] 2.988297e-02 7.711734e-02
#> [56,] 2.988297e-02 7.711734e-02
#> [57,] 2.988297e-02 7.711734e-02
#> [58,] 4.637550e-02 8.627999e-02
#> [59,] 1.000000e+00 1.000000e+00
#> [60,] 4.637550e-02 8.627999e-02
#> [61,] 2.988297e-02 7.711734e-02
#> [62,] 1.000000e+00 1.000000e+00
#> [63,] 1.036419e-05 2.763785e-04
#> [64,] 8.178125e-02 1.486932e-01
#> [65,] 4.637550e-02 8.627999e-02
#> [66,] 2.988297e-02 7.711734e-02
#> [67,] 1.000000e+00 1.000000e+00
#> [68,] 1.000000e+00 1.000000e+00
#> [69,] 1.000000e+00 1.000000e+00
#> [70,] 1.000000e+00 1.000000e+00
#> [71,] 4.637550e-02 8.627999e-02
#> [72,] 4.637550e-02 8.627999e-02
#> [73,] 2.277580e-02 7.711734e-02
#> [74,] 5.455286e-05 4.849143e-04
#> [75,] 5.847200e-07 2.338880e-05
#> [76,] 1.000000e+00 1.000000e+00
#> [77,] 1.000000e+00 1.000000e+00
#> [78,] 4.637550e-02 8.627999e-02
#> [79,] 4.637550e-02 8.627999e-02
#> [80,] 2.988297e-02 7.711734e-02
```
