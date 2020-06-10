
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ROSeq

Modeling expression ranks for noise-tolerant differential expression
analysis of scRNA-Seq data.

## Introduction

ROSeq - A rank based approach to modeling gene expression with filtered
and normalized read count matrix. Takes in the complete filtered and
normalized read count matrix, the location of the two sub-populations
and the number of cores to be used.

## Installation

The developer version of the R package can be installed with the
following R commands:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('ROSeq')
```

or can be installed with the following R commands:

``` r
library(devtools)
install_github('krishan57gupta/ROSeq')
```

## Vignette tutorial

This vignette uses a tung dataset already inbuilt in same package, to
demonstrate a standard pipeline. This vignette can be used as a tutorial
as well.

#### Reference of Tung data:

Tung, P.-Y.et al.Batch effects and the effective design of single-cell
geneexpression studies.Scientific reports7, 39921 (2017).

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

First convert matrix to numeric values, then cell filtering, gene
filtering. After all finally normalization and then tranformation. Note:
For filtering, normalization and tranfromation other methods can be
used, but recomended as shown in example.

``` r
gene_names<-rownames(samples$count)
samples$count<-apply(samples$count,2,function(x) as.numeric(x))
rownames(samples$count)<-gene_names
samples$count<-samples$count[,colSums(samples$count> 0) > 2000]
gkeep<-apply(samples$count,1,function(x) sum(x>2)>=3)
samples$count<-samples$count[gkeep,]
samples$count<-limma::voom(ROSeq::TMMnormalization(samples$count))
```

### ROSeq calling:

Requires a matrix with row as genes and columns and cells, and also
condition of cells, means lables for each cell. numCores can be set as
per number of core/cpu
avaialble.

``` r
output<-ROSeq(countData=samples$count$E, condition = samples$group, numCores=1)
```

### Showing results are in the form of pVals, pAdj and log2FC

##### p\_Vals : p\_value (unadjusted)

##### p\_Adj : Adjusted p-value, based on FDR method

##### log2FC : log fold-chage of the average expression between the two groups,

#### Note:

Positive values show feature is highly enriched in the first group.

``` r
output[1:5,]
#>                     pVals      pAdj      log2FC
#> ENSG00000237683 0.6741425 0.9321651 -0.02240619
#> ENSG00000188976 0.7484244 0.9426495  0.03652966
#> ENSG00000187608 0.2282451 0.8481636  0.15428280
#> ENSG00000188157 0.5138812 0.9082800 -0.06789033
#> ENSG00000131591 0.1235577 0.7438811 -0.07333149
```

## Publication (preprint available at bioRxiv):

Gupta, K., Lalit, M., Biswas, A., Maulik, U., Bandyopadhyay, S., Ahuja,
G., Ghosh, A., and Sengupta, D., 2020. ROSeq: Noise-tolerant
differential expression analysis of scRNA-Seq data. DOI:
<https://doi.org/10.1101/374025>

## All the coding files related to the analysis done in the paper:

Available at github, the link for the same: https://github.com/krishan57gupta/ROSeq/tree/master/doc/paper_analysis_code
