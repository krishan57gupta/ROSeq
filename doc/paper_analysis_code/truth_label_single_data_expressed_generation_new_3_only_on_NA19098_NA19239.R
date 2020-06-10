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
pred_label_founder<-function(L,file,i,j)
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
  R10=ROSeq(countData = limma::voom(TMMnormalization(L$count))$E, condition = L$condt, numCores = 1)[,c(1,2)]
  H_time <- Sys.time()
  R12=DESeq2(L)
  R12$padj=p.adjust(R12$pval, method = "fdr")
  J_time <- Sys.time()
  R13=run_SCDE(L)
  R13$df$padj<-p.adjust(R13$df$pval, method = "fdr")
  K_time <- Sys.time()
  all_time=c(B_time-A_time,
             C_time-B_time,
             D_time-C_time,
             H_time-D_time,
             J_time-H_time,
             K_time-J_time)
  names(all_time)=c("MAST",
                    "BPSC",
                    "Wilcoxon",
                    "ROSeq_0_tmmvoom",
                    "DEseq2",
                    "SCDE")
  saveRDS(all_time,paste("~/ROSeq/paper_work/latest_data/NA19098_NA19239/time/",file,"_",j,"_",i,".rds",sep=""))
  R=cbind(R4$df,R5$df,R6$df,R10,R12,R13$df[,c(1,2)])
  colnames(R)=c("MASTcpm_pval","MASTcpm_padj",
                "BPSC_pval","BPSC_padj",
                "Wilcoxon_pval","Wilcoxon_padj",
                "ROSeq_0_tmmvoom_pval","ROSeq_0_tmmvoom_padj",
                "DEseq2_pval","DEseq2_padj",
                "SCDE_pval","SCDE_padj")
  return(R)
}


# ------------------------------------------------Tung dataset
molecules_1 <- read.table("~/ROSeq/Functions/singleCellSeq-gh-pages/data/molecules.txt", sep = "\t")
molecules=molecules_1[-1,-1]
bulk_counts=apply(molecules,2,as.numeric)
rownames(bulk_counts)=molecules_1[,1][-1]
colnames(bulk_counts)=as.matrix(molecules_1[1,][-1])
dim(bulk_counts)
bulk_counts<-bulk_counts[,colSums(bulk_counts> 0) > 20]
dim(bulk_counts)
anno <- read.table("~/ROSeq/Functions/singleCellSeq-gh-pages/data/annotation.txt", sep = "\t", header = TRUE)[,1]

L_Tung=list()
L_Tung$NA19098_count=bulk_counts[,which(anno %in% "NA19098")]
L_Tung$NA19239_count=bulk_counts[,which(anno %in% "NA19239")]
Truth_Tung <- readRDS("~/ROSeq/Functions/05_Tung_31/Truth_trap_NA19098_NA19239.rds")
Truth_Tung <- Truth_Tung[Truth_Tung=="0"]
comman_gene_id=intersect(rownames(L_Tung$NA19098_count),names(Truth_Tung))
L_Tung$NA19098_count=L_Tung$NA19098_count[which(rownames(L_Tung$NA19098_count) %in% names(Truth_Tung)),]
L_Tung$NA19239_count=L_Tung$NA19239_count[which(rownames(L_Tung$NA19239_count) %in% names(Truth_Tung)),]

L_Tung$NA19098_group=rep(1,dim(L_Tung$NA19098_count)[2])
L_Tung$NA19239_group=rep(2,dim(L_Tung$NA19239_count)[2])


L_Tung$NA19098_NA19239_pred_label=list()
for (j in c(10,25,50,75,100,125,150,175,200,225,250,275))
{
  if(file.exists(paste("~/ROSeq/paper_work/latest_data/NA19098_NA19239/L_Tung_single_expressed_NA19098_NA19239","_",j,".rds",sep="")))
  {
    temp_object<-readRDS(paste("~/ROSeq/paper_work/latest_data/NA19098_NA19239/L_Tung_single_expressed_NA19098_NA19239","_",j,".rds",sep=""))
  L_Tung$NA19098_NA19239_pred_label<-  temp_object$NA19098_NA19239_pred_label
  }
  else
  {
      L_Tung$NA19098_NA19239_pred_label=list()
      }
  if(length(L_Tung$NA19098_NA19239_pred_label)>=100)
  {
    next
  }
  for(i in (length(L_Tung$NA19098_NA19239_pred_label)+1):100)
  {
    print(paste(j,i,sep="_"))
  rand_samp=sample(1:dim(L_Tung$NA19098_count)[2],j)
  rand_samp_1=sample(1:dim(L_Tung$NA19239_count)[2],j)
  L_Tung$NA19098_NA19239_count=cbind(L_Tung$NA19098_count[,c(rand_samp)],L_Tung$NA19239_count[,c(rand_samp_1)])
  L_Tung$NA19098_NA19239_group=c(rep(1,length(L_Tung$NA19098_group[c(rand_samp)])),rep(2,length(L_Tung$NA19239_group[c(rand_samp_1)])))
  
  L=list()
  L$count=L_Tung$NA19098_NA19239_count
  L$condt=L_Tung$NA19098_NA19239_group
  L_Tung$NA19098_NA19239_pred_label[[i]]=pred_label_founder(L,"expressed_NA19098_NA19239_time",i,j)
  saveRDS(L_Tung,paste("~/ROSeq/paper_work/latest_data/NA19098_NA19239/L_Tung_single_expressed_NA19098_NA19239","_",j,".rds",sep=""))
  }
  
}
