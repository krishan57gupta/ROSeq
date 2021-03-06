% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROSeq.R
\name{ROSeq}
\alias{ROSeq}
\title{Modeling expression ranks for noise-tolerant differential 
expression analysis of scRNA-Seq data}
\usage{
ROSeq(countData, condition, numCores = 1)
}
\arguments{
\item{countData}{The normalised and filtered, read count matrix, with
row names as genes name/ID and column names as sample id/name}

\item{condition}{Labels for the two sub-populations}

\item{numCores}{The number of cores to be used}
}
\value{
pValues and FDR adjusted p significance values
}
\description{
Takes in the complete filtered and normalized read count
matrix, the location of the two sub-populations and the number of cores
to be used
}
\examples{
countData<-list()
countData$count<-ROSeq::L_Tung_single$NA19098_NA19101_count
countData$group<-ROSeq::L_Tung_single$NA19098_NA19101_group
head(countData$count)
gene_names<-rownames(countData$count)
countData$count<-apply(countData$count,2,function(x) as.numeric(x))
rownames(countData$count)<-gene_names
countData$count<-countData$count[,colSums(countData$count> 0) > 2000]
g_keep <- apply(countData$count,1,function(x) sum(x>2)>=3)
countData$count<-countData$count[g_keep,]
countData$count<-limma::voom(ROSeq::TMMnormalization(countData$count))
output<-ROSeq(countData=countData$count$E, condition = countData$group)
output
}
