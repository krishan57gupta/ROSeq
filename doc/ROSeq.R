## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ROSeq)
library(compcodeR)
library(edgeR)

## ----data, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE------------
n_samples=8
n_genes=80
n_de_genes=40
mydata <- generateSyntheticData(dataset = "mydata", n.vars = n_genes, samples.per.cond = n_samples, n.diffexp = n_de_genes)
samples=list()
samples$count=mydata@count.matrix
samples$group=mydata@sample.annotations$condition
head(samples$count)

## ----preprocesing, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE----
samples$count=apply(samples$count,2,function(x) as.numeric(x))
gkeep <- apply(samples$count,1,function(x) sum(x>0)>5)
samples$count=samples$count[gkeep,]
samples$count=edgeR::cpm(samples$count)

## ----main, message=FALSE,warning = FALSE, include=TRUE, cache=FALSE-----------
output=ROSeq(countData=samples$count, condition = samples$group, numCores=1)

## ----output, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE----------
output

