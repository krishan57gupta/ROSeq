median_normalization <- function(data){
  cname=colnames(data)
  rname=rownames(data)
  normalized_matrix = as.matrix(data)
  total_count=colSums(data)
  med_total=median(total_count)
  for(i in 1:ncol(normalized_matrix))
  {
    normalized_matrix[,i] = data[,i] * (med_total/total_count[i])
  }
  colnames(normalized_matrix)=cname
  rownames(normalized_matrix)=rname
  return(normalized_matrix)
}
TMMnormalization <- function(countTable){
  ## TMM normalization based on edgeR package:
  require("edgeR")
  cname=colnames(countTable)
  rname=rownames(countTable)
  nf=calcNormFactors(countTable ,method= "TMM")
  nf= colSums(countTable)*nf
  scalingFactors = nf/mean(nf)
  countTableTMM <- t(t(countTable)/scalingFactors)
  colnames(countTableTMM)=cname
  rownames(countTableTMM)=rname
  return(countTableTMM)
}
gene_filter <- function(data,min.count=5,min.cell=5){
  express=apply(data,1,function(x) sum(x>min.count)>=(min.cell))
  print(paste('Filtered Genes:',sum(express),sep = " "))
  return(data[express,])
}
DESeq2 <- function(L){
  metaData=data.frame("id"=colnames(L$count),"dex"=L$condt)
  rownames(metaData)=colnames(L$count)
  dds=DESeqDataSetFromMatrix(countData=L$count, 
                             colData=metaData, 
                             design=~dex)
  dds <- DESeq(dds)
  res=results(dds)[,5:6]
  colnames(res)<-c("pval","padj")
  return(res)
}
exe_time=c()
pred_label_founder<-function(L,file)
{
  L$count=gene_filter(L$count,min.count=2,min.cell=3)
  colnames(L$count)=paste(colnames(L$count),1:dim(L$count)[2],sep="_")
  library(ROSeq)
  library(Linnorm)
  library(limma)
  library(edgeR)
  library(BPSC)
  library(MAST)
  library(DESeq2, lib.loc = "/usr/local/lib/R/site-library")
  library(DESeq2)
  source("~/ROSeq/DE_Methods/apply_edgeRQLF.R")
  source("~/ROSeq/DE_Methods/apply_limmatrend.R")
  source("~/ROSeq/DE_Methods/apply_voomlimma.R")
  source("~/ROSeq/DE_Methods/apply_MASTcpm.R")
  source("~/ROSeq/DE_Methods/apply_BPSC.R")
  source("~/ROSeq/DE_Methods/apply_Wilcoxon.R")
  source("~/ROSeq/DE_Methods/apply_SCDE.R")
  A_time <- Sys.time()
  R4=run_MASTcpm(L)
  R4$df$padj<-p.adjust(R4$df$pval, method = "fdr")
  B_time <- Sys.time()
  R5=run_BPSC(L)
  R5$df$padj<-p.adjust(R5$df$pval, method = "fdr")
  C_time <- Sys.time()
  R6=run_Wilcoxon(L)
  R6$df$padj<-p.adjust(R6$df$pval, method = "fdr")
  D_time <- Sys.time()
  R7=ROSeq(countData = Linnorm(L$count), condition = L$condt, numCores = 8)[,c(1,2)]
  E_time <- Sys.time()
  R8=ROSeq(countData = edgeR::cpm(L$count), condition = L$condt, numCores = 8)[,c(1,2)]
  F_time <- Sys.time()
  R9=ROSeq(countData = TMMnormalization(L$count), condition = L$condt, numCores = 8)[,c(1,2)]
  G_time <- Sys.time()
  R10=ROSeq(countData = limma::voom(TMMnormalization(L$count))$E, condition = L$condt, numCores = 8)[,c(1,2)]
  H_time <- Sys.time()
  R11=ROSeq(countData = median_normalization(L$count), condition = L$condt, numCores = 8)[,c(1,2)]
  I_time <- Sys.time()
  R12=DESeq2(L)
  R12$padj=p.adjust(R12$pval, method = "fdr")
  J_time <- Sys.time()
  R13=run_SCDE(L)
  R13$df$padj<-p.adjust(R13$df$pval, method = "fdr")
  K_time <- Sys.time()
  all_time=c(B_time-A_time,
             C_time-B_time,
             D_time-C_time,
             E_time-D_time,
             F_time-E_time,
             G_time-F_time,
             H_time-G_time,
             I_time-H_time,
             J_time-I_time,
             K_time-J_time)
  names(all_time)=c("MAST",
                "BPSC",
                "Wilcoxon",
                "ROSeq_0_linnorm",
                "ROSeq_0_cpm",
                "ROSeq_0_tmm",
                "ROSeq_0_tmmvoom",
                "ROSeq_0_median",
                "DEseq2",
                "SCDE")
  saveRDS(all_time,paste("~/ROSeq/paper_work/latest_data/pair/",file,"_","time",".rds",sep=""))
  R=cbind(R4$df,R5$df,R6$df,R7,R8,R9,R10,R11,R12,R13$df[,c(1,2)])
  colnames(R)=c("MASTcpm_pval","MASTcpm_padj",
                "BPSC_pval","BPSC_padj",
                "Wilcoxon_pval","Wilcoxon_padj",
                "ROSeq_0_linnorm_pval","ROSeq_0_linnorm_padj",
                "ROSeq_0_cpm_pval","ROSeq_0_cpm_padj",
                "ROSeq_0_tmm_pval","ROSeq_0_tmm_padj",
                "ROSeq_0_tmmvoom_pval","ROSeq_0_tmmvoom_padj",
                "ROSeq_0_median_pval","ROSeq_0_median_padj",
                "DEseq2_pval","DEseq2_padj",
                "SCDE_pval","SCDE_padj")
  return(R)
}



# ------------------------------------------------Trap dataset

load("~/ROSeq/paper_work/latest_data/pair/data2.RData")
bulk_pseudo_count=good_pseudo_count
name=rownames(bulk_pseudo_count)
bulk_pseudo_count=apply(bulk_pseudo_count, 2, as.numeric)
rownames(bulk_pseudo_count)<- name
dim(bulk_pseudo_count)
bulk_pseudo_count<-bulk_pseudo_count[,colSums(bulk_pseudo_count> 0) > 2000]
dim(bulk_pseudo_count)
L_Trap=list()
L_Trap$T0_count=bulk_pseudo_count[,which(substr(colnames(bulk_pseudo_count),1,2)%in% "T0")]
L_Trap$T24_count=bulk_pseudo_count[,which(substr(colnames(bulk_pseudo_count),1,2)%in% "T2")]
L_Trap$T0_group=rep(1,dim(L_Trap$T0_count)[2])
L_Trap$T24_group=rep(2,dim(L_Trap$T24_count)[2])
L_Trap$T0_T24_count=cbind(L_Trap$T0_count,L_Trap$T24_count)
L_Trap$T0_T24_group=c(L_Trap$T0_group,L_Trap$T24_group)



L=list()
L$count=L_Trap$T0_T24_count
L$condt=L_Trap$T0_T24_group
L_Trap$T0_T24_pred_label=pred_label_founder(L,"L_Trap_T0_T24")

saveRDS(L_Trap,"~/ROSeq/paper_work/latest_data/pair/L_Trap_single.rds")


# ------------------------------------------------Tung dataset
molecules_1 <- read.table("~/ROSeq/Functions/singleCellSeq-gh-pages/data/molecules.txt", sep = "\t")
molecules=molecules_1[-1,-1]
bulk_counts=apply(molecules,2,as.numeric)
rownames(bulk_counts)=molecules_1[,1][-1]
colnames(bulk_counts)=as.matrix(molecules_1[1,][-1])
dim(bulk_counts)
bulk_counts<-bulk_counts[,colSums(bulk_counts> 0) > 2000]
dim(bulk_counts)
anno <- read.table("~/ROSeq/Functions/singleCellSeq-gh-pages/data/annotation.txt", sep = "\t", header = TRUE)[,1]

L_Tung=list()
L_Tung$NA19098_count=bulk_counts[,which(anno %in% "NA19098")]
L_Tung$NA19101_count=bulk_counts[,which(anno %in% "NA19101")]
L_Tung$NA19239_count=bulk_counts[,which(anno %in% "NA19239")]
L_Tung$NA19098_group=rep(1,dim(L_Tung$NA19098_count)[2])
L_Tung$NA19101_group=rep(2,dim(L_Tung$NA19101_count)[2])
L_Tung$NA19239_group=rep(3,dim(L_Tung$NA19239_count)[2])
L_Tung$NA19098_NA19101_count=cbind(L_Tung$NA19098_count,L_Tung$NA19101_count)
L_Tung$NA19098_NA19101_group=c(L_Tung$NA19098_group,L_Tung$NA19101_group)
L_Tung$NA19098_NA19239_count=cbind(L_Tung$NA19098_count,L_Tung$NA19239_count)
L_Tung$NA19098_NA19239_group=c(L_Tung$NA19098_group,L_Tung$NA19239_group)
L_Tung$NA19239_NA19101_count=cbind(L_Tung$NA19239_count,L_Tung$NA19101_count)
L_Tung$NA19239_NA19101_group=c(L_Tung$NA19239_group,L_Tung$NA19101_group)

L=list()
L$count=L_Tung$NA19098_NA19101_count
L$condt=L_Tung$NA19098_NA19101_group
L_Tung$NA19098_NA19101_pred_label=pred_label_founder(L,"L_Tung_NA19098_NA19101")

L=list()
L$count=L_Tung$NA19098_NA19239_count
L$condt=L_Tung$NA19098_NA19239_group
L_Tung$NA19098_NA19239_pred_label=pred_label_founder(L,"L_Tung_NA19098_NA19239")

L=list()
L$count=L_Tung$NA19239_NA19101_count
L$condt=L_Tung$NA19239_NA19101_group
L_Tung$NA19239_NA19101_pred_label=pred_label_founder(L,"L_Tung_NA19239_NA19101")

saveRDS(L_Tung,"~/ROSeq/paper_work/latest_data/pair/L_Tung_single.rds")

# # ------------------------------------------------simulated data
# library(SimSeq)
# data(kidney)
# counts <- kidney$counts # Matrix of read counts from KIRC dataset
# replic <- kidney$replic # Replic vector indicating paired columns
# treatment <- kidney$treatment # Treatment vector indicating Non-Tumor or Tumor columns
# nf <- apply(counts, 2, quantile, 0.75)
# library(fdrtool)
# ## Not run:
# ### Example 1: Simulate Matrix with 1000 DE genes and 4000 EE genes
# data.sim <- SimData(counts = counts, replic = replic, treatment = treatment,
#                     sort.method = "paired", k.ind = 30, n.genes = 10000, n.diff = 1000,
#                     norm.factors = nf)
# L_SyntheticData_10=list()
# L_SyntheticData_10$count=data.sim$counts
# colnames(L_SyntheticData_10$count)<-paste("C",1:dim(L_SyntheticData_10$count)[2],sep="_")
# L_SyntheticData_10$group=data.sim$treatment
# L_SyntheticData_10$truth_label=data.sim$DE.ind
# L_SyntheticData_10$truth_label[L_SyntheticData_10$truth_label==TRUE]<-1
# L_SyntheticData_10$truth_label[L_SyntheticData_10$truth_label==FALSE]<-2
# L_SyntheticData_10$truth_label=L_SyntheticData_10$truth_label-1
# names(L_SyntheticData_10$truth_label)<-rownames(L_SyntheticData_10$count)
# 
# L=list()
# L$count=L_SyntheticData_10$count
# L$condt=L_SyntheticData_10$group
# L_SyntheticData_10$pred_label=pred_label_founder(L,"L_SyntheticData_10")
# saveRDS(L_SyntheticData_10,"~/ROSeq/paper_work/latest_data/pair/L_SyntheticData_10.rds")
# 
# 
# # ------------------------------------------------simulated data
# 
# library(compcodeR)
# n_samples=100
# n_genes=10000
# n_de_genes=1000
# mydata <- generateSyntheticData(dataset = "mydata", n.vars = n_genes, samples.per.cond = n_samples, n.diffexp = n_de_genes)
# mydata@count.matrix<-mydata@count.matrix[rowSums(mydata@count.matrix> 0) > 5,]
# L_SyntheticData_100=list()
# L_SyntheticData_100$count=mydata@count.matrix
# L_SyntheticData_100$group=mydata@sample.annotations$condition
# L_SyntheticData_100$truth_label=abs(mydata@variable.annotations$differential.expression-1)
# names(L_SyntheticData_100$truth_label)<-rownames(L_SyntheticData_100$count)
# 
# 
# L=list()
# L$count=L_SyntheticData_100$count
# L$condt=L_SyntheticData_100$group
# L_SyntheticData_100$pred_label=pred_label_founder(L,"L_SyntheticData_100")
# saveRDS(L_SyntheticData_100,"~/ROSeq/paper_work/latest_data/pair/L_SyntheticData_100.rds")


