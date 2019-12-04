## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ROSeq)
library(compcodeR)
library(edgeR)

## ----functions, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE-------
TMMnormalization <- function(countTable){
  nf=calcNormFactors(countTable ,method= "TMM")
  nf= colSums(countTable)*nf
  scalingFactors = nf/mean(nf)
  countTableTMM <- t(t(countTable)/scalingFactors)
  return(countTableTMM)
}

## ----data, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE------------
n_samples=10
n_genes=100
n_de_genes=50
mydata <- generateSyntheticData(dataset = "mydata", n.vars = n_genes, samples.per.cond = n_samples, n.diffexp = n_de_genes)
samples=list()
samples$count=mydata@count.matrix
samples$group=mydata@sample.annotations$condition
head(samples$count)

## ----preprocesing, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE----
gkeep <- apply(samples$count,1,function(x) sum(x>0)>5)
samples$count=samples$count[gkeep,]
samples$count=TMMnormalization(samples$count)

