
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ROSeq

Modeling expression ranks for noise-tolerant differential expression
analysis of scRNA-Seq data

## Introduction

ROSeq - A rank based approach to modeling gene expression with filtered
and normalized read count matrix. ROSeq takes filtered and normalized
read matrix and cell-annotation/condition as input and determines the
differentially expressed genes between the contrasting groups of single
cells. One of the input parameters is the number of cores to be used.

## Installation

The developer’s version of the R package can be installed with the
following R commands:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("ROSeq")
```

The github’s version of the R package can be installed with the
following R commands:

``` r
library(devtools)
install_github('krishan57gupta/ROSeq')
```

## Vignette tutorial

This vignette uses the Tung dataset, which is already inbuilt in the
package, to demonstrate a standard pipeline.

## Example

Libraries need to be loaded before running.

``` r
library(ROSeq)
library(edgeR)
#> Loading required package: limma
library(limma)
```

### Loading tung dataset

``` r
samples<-list()
samples$count<-ROSeq::L_Tung_single$NA19098_NA19101_count
samples$group<-ROSeq::L_Tung_single$NA19098_NA19101_group
samples$count[1:5,1:5]
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

### Data Preprocessing:

#### Cells and genes filtering then voom transformation after TMM normalization

Below commands can be used for Cell/gene filtering, TMM normalization
and voom transformation. The user is free to use an alternative
preprocessing strategy while using different filtering/normalization
methods.

``` r
gene_names<-rownames(samples$count)
samples$count<-apply(samples$count,2,function(x) as.numeric(x))
rownames(samples$count)<-gene_names
samples$count<-samples$count[,colSums(samples$count> 0) > 2000]
gkeep<-apply(samples$count,1,function(x) sum(x>2)>=3)
samples$count<-samples$count[gkeep,]
samples$count<-limma::voom(ROSeq::TMMnormalization(samples$count))
```

### ROSeq analysis.

Input: gene expression matrix with genes in rows and cells in columns.
Condition/group annotation of cells also need to be supplied. User can
set numCores based the hardware specifications in her
computer.

``` r
output<-ROSeq(countData=samples$count$E, condition = samples$group, numCores=1)
```

### Showing results are in the form of pVals and pAdj

##### p\_Vals : p\_value (unadjusted)

##### p\_Adj : Adjusted p-value, based on FDR method

``` r
output[1:5,]
#>                     pVals      pAdj
#> ENSG00000237683 0.6741425 0.9321651
#> ENSG00000188976 0.7484244 0.9426495
#> ENSG00000187608 0.2282451 0.8481636
#> ENSG00000188157 0.5138812 0.9082800
#> ENSG00000131591 0.1235577 0.7438811
```
