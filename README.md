
<!-- README.md is generated from README.Rmd. Please edit that file -->

ROSeq - A rank based approach to modeling gene expression Author:
Krishan Gupta

date: 2/12/2019

Desktop Installation
===============

The developer version of the R package can be installed with the following R commands:

``` r
library(devtools)
install_github("krishan57gupta/ROSeq")
```

Vignette tutorial
------------------
This vignette uses a small data set of simulated data to demonstrate a standard pipeline. This vignette can be used as a tutorial as well.

Setting up directories
----------------------

``` r
library(ROSeq)
```


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
n_samples=10
n_genes=100
n_de_genes=50
mydata <- generateSyntheticData(dataset = "mydata", n.vars = n_genes, samples.per.cond = n_samples, n.diffexp = n_de_genes)
samples=list()
samples$count=mydata@count.matrix
samples$group=mydata@sample.annotations$condition
head(samples$count)
#>    sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9
#> g1    1356    1245     803    1616    2999    1205    1055    2130    3289
#> g2     158    4767   23868   14287   24858   20848    3381    2797   36909
#> g3  162384  225911  114693  245087  197813  155566  176720  194967  221882
#> g4  356391  445022  348340  345196  476913  315737  412461  415997  389053
#> g5  156234  224224   76847  140907  177125  189574  105330  163920   96003
#> g6  126276  272071  194994  142560  246846  192265  166681  253189  375595
#>    sample10 sample11 sample12 sample13 sample14 sample15 sample16 sample17
#> g1     1592      209     4673     4247     3757     6550     4387     8478
#> g2    12184    18741    35068    17810    25154    13972     4360   128737
#> g3   260007   103582   317782   378701   274371   185898   463242   275008
#> g4   571280   223488   540140   422815   333927   301526   500482   393359
#> g5   280611   146903   193362   123534   178312   108266   203832   244179
#> g6   355536   179789   309571   209285   191192   217913   257002   322348
#>    sample18 sample19 sample20
#> g1     6483     5326     2732
#> g2    28149     5985    35416
#> g3   409651   258623   255804
#> g4   424278   405219   444973
#> g5   271481   162803   253701
#> g6   174805   384250   245187
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
#>             pVals         pAdj
#> g1   1.000000e+00 1.000000e+00
#> g2   2.966035e-03 4.893958e-02
#> g3   6.494791e-01 1.000000e+00
#> g4   1.000000e+00 1.000000e+00
#> g5   8.305709e-01 1.000000e+00
#> g6   1.000000e+00 1.000000e+00
#> g7   5.342199e-06 2.644389e-04
#> g8   3.748528e-02 4.123381e-01
#> g9   7.479758e-01 1.000000e+00
#> g10  4.316162e-01 1.000000e+00
#> g11  8.305709e-01 1.000000e+00
#> g12  1.000000e+00 1.000000e+00
#> g13  1.000000e+00 1.000000e+00
#> g14  2.341152e-01 8.277645e-01
#> g15  1.377423e-01 8.277645e-01
#> g16  1.000000e+00 1.000000e+00
#> g17  2.341152e-01 8.277645e-01
#> g18  8.305709e-01 1.000000e+00
#> g19  1.000000e+00 1.000000e+00
#> g20  1.000000e+00 1.000000e+00
#> g21  1.000000e+00 1.000000e+00
#> g22  1.000000e+00 1.000000e+00
#> g23  1.338080e-02 1.655874e-01
#> g24  8.305709e-01 1.000000e+00
#> g25  1.000000e+00 1.000000e+00
#> g26  7.479758e-01 1.000000e+00
#> g27  1.000000e+00 1.000000e+00
#> g28  4.316162e-01 1.000000e+00
#> g29  1.000000e+00 1.000000e+00
#> g30  1.000000e+00 1.000000e+00
#> g31  7.479758e-01 1.000000e+00
#> g32  4.316162e-01 1.000000e+00
#> g33  1.000000e+00 1.000000e+00
#> g34  1.000000e+00 1.000000e+00
#> g35  7.479758e-01 1.000000e+00
#> g36  2.966035e-03 4.893958e-02
#> g37  1.000000e+00 1.000000e+00
#> g38  1.000000e+00 1.000000e+00
#> g39  1.000000e+00 1.000000e+00
#> g40  3.804410e-01 1.000000e+00
#> g41  1.000000e+00 1.000000e+00
#> g42  1.000000e+00 1.000000e+00
#> g43  7.479758e-01 1.000000e+00
#> g44  9.516552e-01 1.000000e+00
#> g45  1.338080e-02 1.655874e-01
#> g46  4.935994e-01 1.000000e+00
#> g47  8.305709e-01 1.000000e+00
#> g48  2.031778e-01 8.277645e-01
#> g49  1.000000e+00 1.000000e+00
#> g50  1.000000e+00 1.000000e+00
#> g51  1.000000e+00 1.000000e+00
#> g52  1.000000e+00 1.000000e+00
#> g53  8.255709e-02 6.035047e-01
#> g54  4.992210e-01 1.000000e+00
#> g55  9.975115e-04 3.291788e-02
#> g56  1.000000e+00 1.000000e+00
#> g57  2.341152e-01 8.277645e-01
#> g58  1.000000e+00 1.000000e+00
#> g59  1.000000e+00 1.000000e+00
#> g60  9.144010e-02 6.035047e-01
#> g61  7.479758e-01 1.000000e+00
#> g62  8.305709e-01 1.000000e+00
#> g63  2.341152e-01 8.277645e-01
#> g64  1.000000e+00 1.000000e+00
#> g65  2.031778e-01 8.277645e-01
#> g66  8.305709e-01 1.000000e+00
#> g67  4.241277e-02 4.198864e-01
#> g68  1.000000e+00 1.000000e+00
#> g69  2.341152e-01 8.277645e-01
#> g70  8.305709e-01 1.000000e+00
#> g71  7.479758e-01 1.000000e+00
#> g72  9.144010e-02 6.035047e-01
#> g73  2.966035e-03 4.893958e-02
#> g74  7.479758e-01 1.000000e+00
#> g75  7.479758e-01 1.000000e+00
#> g76  3.804410e-01 1.000000e+00
#> g77  7.479758e-01 1.000000e+00
#> g78  5.843986e-02 4.821289e-01
#> g79  7.479758e-01 1.000000e+00
#> g80  2.076129e-01 8.277645e-01
#> g81  2.341152e-01 8.277645e-01
#> g82  1.000000e+00 1.000000e+00
#> g83  7.479758e-01 1.000000e+00
#> g84  7.479758e-01 1.000000e+00
#> g85  7.479758e-01 1.000000e+00
#> g86  1.000000e+00 1.000000e+00
#> g87  8.305709e-01 1.000000e+00
#> g88  1.000000e+00 1.000000e+00
#> g89  2.341152e-01 8.277645e-01
#> g90  7.479758e-01 1.000000e+00
#> g91  1.000000e+00 1.000000e+00
#> g92  4.992210e-01 1.000000e+00
#> g94  1.000000e+00 1.000000e+00
#> g95  1.000000e+00 1.000000e+00
#> g96  6.339757e-01 1.000000e+00
#> g97  5.843986e-02 4.821289e-01
#> g98  5.117139e-07 5.065967e-05
#> g99  2.031778e-01 8.277645e-01
#> g100 2.031778e-01 8.277645e-01
```
