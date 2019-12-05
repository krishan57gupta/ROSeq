
<!-- README.md is generated from README.Rmd. Please edit that file -->

ROSeq - A rank based approach to modeling gene expression Author:
Krishan Gupta

date: 2/12/2019

## Introduction

ROSeq - A rank based approach to modeling gene expression with filtered
and normalized read count matrix. Takes in the complete filtered and
normalized read count matrix, the location of the two sub-populations
and the number of cores to be used.

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
#> g1   13388    5354    2912    9852    6690    3106    7636    4905    9355
#> g2   54976  114733  272759  185489  129384  270818  147838  231502  181089
#> g3  119381   57112   88210   59239   81405  126708   96010   45177   71887
#> g4    5886    6247    5723    9829    6338   13715    9939    8859    9190
#> g5   63909   51210   77790   39713   39736   63126   92695  192066   78064
#> g6   74868   58381   21519   47324   45029   30004   54013   71637  135862
#>    sample10 sample11 sample12 sample13 sample14 sample15 sample16
#> g1    14611    23262     8729    11512    15253     7098     7660
#> g2   116589   183259   152345   280733   461248   414182  1043036
#> g3   191901    85463    59813   142915   107260    73348   119712
#> g4    12645    16255     6732    14732    10426     8290    11595
#> g5    85685   145385    57978    95292    82625    59602   107131
#> g6    51141   285782    69757    36207    85190   131683   114503
```

``` r
gkeep <- apply(samples$count,1,function(x) sum(x>0)>5)
samples$count=samples$count[gkeep,]
samples$count=TMMnormalization(samples$count)
```

``` r
output=ROSeq(countData=samples$count, condition = samples$group, numCores=1)
```

``` r
output
#>            pVals         pAdj
#> g1  1.000000e+00 1.0000000000
#> g2  8.810070e-01 1.0000000000
#> g3  1.000000e+00 1.0000000000
#> g4  1.000000e+00 1.0000000000
#> g5  8.810070e-01 1.0000000000
#> g6  3.612172e-01 1.0000000000
#> g7  6.571580e-04 0.0075103774
#> g8  4.824637e-01 1.0000000000
#> g9  1.000000e+00 1.0000000000
#> g10 8.810070e-01 1.0000000000
#> g11 1.000000e+00 1.0000000000
#> g12 4.824637e-01 1.0000000000
#> g13 5.853661e-04 0.0075103774
#> g14 1.000000e+00 1.0000000000
#> g15 2.535364e-01 0.8113166378
#> g16 1.000000e+00 1.0000000000
#> g17 1.000000e+00 1.0000000000
#> g18 1.109624e-02 0.0467210294
#> g19 8.810070e-01 1.0000000000
#> g20 1.109624e-02 0.0467210294
#> g21 1.000000e+00 1.0000000000
#> g22 1.000000e+00 1.0000000000
#> g23 1.000000e+00 1.0000000000
#> g24 1.109624e-02 0.0467210294
#> g25 5.853661e-04 0.0075103774
#> g26 2.245272e-02 0.0898108804
#> g27 1.109624e-02 0.0467210294
#> g28 1.000000e+00 1.0000000000
#> g29 1.000000e+00 1.0000000000
#> g30 1.000000e+00 1.0000000000
#> g31 1.000000e+00 1.0000000000
#> g32 1.200996e-03 0.0120099552
#> g33 1.000000e+00 1.0000000000
#> g34 1.000000e+00 1.0000000000
#> g35 1.109624e-02 0.0467210294
#> g36 1.376253e-04 0.0055050120
#> g37 1.000000e+00 1.0000000000
#> g38 1.000000e+00 1.0000000000
#> g39 8.810070e-01 1.0000000000
#> g40 1.000000e+00 1.0000000000
#> g41 8.810070e-01 1.0000000000
#> g42 2.944461e-02 0.1121699388
#> g43 1.109624e-02 0.0467210294
#> g44 1.217270e-06 0.0000973816
#> g45 1.000000e+00 1.0000000000
#> g46 2.535364e-01 0.8113166378
#> g47 1.000000e+00 1.0000000000
#> g48 1.000000e+00 1.0000000000
#> g49 1.000000e+00 1.0000000000
#> g50 1.000000e+00 1.0000000000
#> g51 1.000000e+00 1.0000000000
#> g52 1.000000e+00 1.0000000000
#> g53 1.109624e-02 0.0467210294
#> g54 1.000000e+00 1.0000000000
#> g55 1.000000e+00 1.0000000000
#> g56 1.000000e+00 1.0000000000
#> g57 1.000000e+00 1.0000000000
#> g58 5.903587e-04 0.0075103774
#> g59 1.109624e-02 0.0467210294
#> g60 1.109624e-02 0.0467210294
#> g61 1.000000e+00 1.0000000000
#> g62 1.000000e+00 1.0000000000
#> g63 7.836796e-02 0.2725842126
#> g64 8.810070e-01 1.0000000000
#> g65 7.836796e-02 0.2725842126
#> g66 1.109624e-02 0.0467210294
#> g67 3.612172e-01 1.0000000000
#> g68 1.000000e+00 1.0000000000
#> g69 1.000000e+00 1.0000000000
#> g70 8.810070e-01 1.0000000000
#> g71 1.000000e+00 1.0000000000
#> g72 1.109624e-02 0.0467210294
#> g73 2.801196e-01 0.8619065721
#> g74 1.000000e+00 1.0000000000
#> g75 1.000000e+00 1.0000000000
#> g76 8.810070e-01 1.0000000000
#> g77 5.853661e-04 0.0075103774
#> g78 1.000000e+00 1.0000000000
#> g79 1.000000e+00 1.0000000000
#> g80 8.810070e-01 1.0000000000
```
