## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ROSeq)
library(compcodeR)
library(edgeR)
library(limma)

## ----data, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE------------
samples<-list()
samples$count<-ROSeq::L_Tung_single$NA19098_NA19101_count
samples$group<-ROSeq::L_Tung_single$NA19098_NA19101_group
samples$count[1:5,1:5]

## ----preprocesing, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE----
samples$count<-apply(samples$count,2,function(x) as.numeric(x))
gkeep <- apply(samples$count,1,function(x) sum(x>0)>5)
samples$count<-samples$count[gkeep,]
samples$count<-limma::voom(ROSeq::TMMnormalization(samples$count))

## ----main, message=FALSE,warning = FALSE, include=TRUE, cache=FALSE-----------
output<-ROSeq(countData=samples$count, condition = samples$group, nbits=0, numCores=1)

## ----output, message=FALSE,warning = FALSE,include=TRUE, cache=FALSE----------
output[1:5,]

