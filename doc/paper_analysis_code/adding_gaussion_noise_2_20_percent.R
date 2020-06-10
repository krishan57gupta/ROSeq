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
pred_label_founder<-function(L,file,i,j,ii,jj)
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
  E_time <- Sys.time()
  R12=DESeq2(L)
  R12$padj=p.adjust(R12$pval, method = "fdr")
  F_time <- Sys.time()
  R13=run_SCDE(L)
  R13$df$padj<-p.adjust(R13$df$pval, method = "fdr")
  G_time <- Sys.time()
  all_time=c(B_time-A_time,
             C_time-B_time,
             D_time-C_time,
             E_time-D_time,
             F_time-E_time,
             G_time-F_time
             )
  names(all_time)=c("MAST",
                    "BPSC",
                    "Wilcoxon",
                    "ROSeq_0_tmmvoom",
                    "DEseq2",
                    "SCDE"
                    )
  saveRDS(all_time,paste("~/ROSeq/paper_work/latest_data/noise_jurkat/time/",file,"_",i,"_",j,"_",ii,"_",jj,".rds",sep=""))
  R=cbind(R4$df,
          R5$df,
          R6$df,
          R10,
          R12,
          R13$df[,c(1,2)]
          )
  colnames(R)=c("MASTcpm_pval","MASTcpm_padj",
                "BPSC_pval","BPSC_padj",
                "Wilcoxon_pval","Wilcoxon_padj",
                "ROSeq_0_tmmvoom_pval","ROSeq_0_tmmvoom_padj",
                "DEseq2_pval","DEseq2_padj",
                "SCDE_pval","SCDE_padj"
                )
  return(R)
}


# ------------------------------------------------Tung dataset
library(Matrix)
for (k in c(25,50,75,100,125,150,175,200,225,250))
{
  if(file.exists(paste("~/ROSeq/paper_work/latest_data/noise_jurkat/null_noise_2_20","_",k,".rds",sep="")))
  { 
    next
   }
  else
    {
      result=list()
      result$gene=list()
      result$sample=list()
      library(Matrix)
      filtered_data=readMM("~/ROSeq/paper_work/latest_data/jurkat_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/matrix.mtx")
      dim(filtered_data)
      r_name=read.table("~/ROSeq/paper_work/latest_data/jurkat_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/genes.tsv")
      c_name=read.table("~/ROSeq/paper_work/latest_data/jurkat_filtered_gene_bc_matrices/filtered_matrices_mex/hg19/barcodes.tsv")
      dim(r_name)
      dim(c_name)
      rownames(filtered_data)=r_name[,2]
      colnames(filtered_data)=c_name[,1]
      # filtered_data=read.table("~/ROSeq/paper_work/latest_data/jurkat_FiRE_data/jurkat_two_species_1580.txt")[,1:1540]
      # dim(filtered_data)
      # c_name=read.table("~/ROSeq/paper_work/latest_data/jurkat_FiRE_data/labels_jurkat_two_species_1580.txt")
      # dim(c_name)
      # rownames(filtered_data)=paste("Gene",1:dim(filtered_data)[1],sep="_")
      # colnames(filtered_data)=paste("Cell",1:1540,c_name[1:1540,1],sep="_")
      filtered_data=filtered_data[,sample(1:dim(filtered_data)[2],k*4)]
      filtered_data=apply(filtered_data,2,function(x) as.integer(x))
      data=apply(filtered_data,2,function(x) as.numeric(x))
      print(dim(data))
      data=data[,colSums(data> 0) > 2000]
      print(dim(data))
      F_G=sample(1:dim(data)[2],k)
      S_G=sample(1:dim(data)[2],k)
      for (ii in c(2,4))
      {
        for (jj in c(.2,.4))
        {
          # on genes
          print(paste(k,ii,jj,sep="_"))
          
          F_D=data[,F_G]                
          S_D=data[,S_G]
          F_D=data[,F_G]                
          S_D=data[,S_G]
          random_subgroup=sample(1:2,1)
          if(random_subgroup==1)
          {
            # F_sum=colSums(F_D) 
            # F_med=median(F_sum)
            # F_D=t((t(F_D)/F_sum)*F_med)
            print(sum(is.na(S_D)))
            for(i in 1:dim(F_D)[1])
            {
              zero_samples=sample(0,length(F_D[i,]),replace = TRUE)
              random_10_percent_index=sample(1:length(F_D[i,]),as.integer(length(F_D[i,])/5))
              random_10_percent_sample=rnorm(as.integer(length(F_D[i,])/5),mean=(mean(F_D[i,])+sd(F_D[i,])*ii),sd=sd(F_D[i,])*jj)
              zero_samples[random_10_percent_index]<-random_10_percent_sample
              F_D[i,]=F_D[i,]+zero_samples
            }
            # F_D=t((t(F_D)*F_sum)/F_med)
            F_D=apply(F_D,2,function(x) as.numeric(as.integer(x)))
            F_D[F_D<0]<-0
            print(sum(is.na(F_D)))
          }
          if(random_subgroup==2)
          {
            # S_sum=colSums(F_D)
            # S_med=median(S_sum)
            # S_D=t((t(S_D)/S_sum)*S_med)
            print(sum(is.na(F_D)))
            for(i in 1:dim(S_D)[1])
            {
              zero_samples=sample(0,length(S_D[i,]),replace = TRUE)
              random_10_percent_index=sample(1:length(S_D[i,]),as.integer(length(S_D[i,])/5))
              random_10_percent_sample=rnorm(as.integer(length(S_D[i,])/5),mean=(mean(S_D[i,])+sd(S_D[i,])*ii),sd=sd(S_D[i,])*jj)
              zero_samples[random_10_percent_index]<-random_10_percent_sample
              S_D[i,]=S_D[i,]+zero_samples
            }
            # S_D=t((t(S_D)*S_sum)/S_med)
            S_D=apply(S_D,2,function(x) as.numeric(as.integer(x)))
            S_D[S_D<0]<-0
            print(sum(is.na(S_D)))
          }
          L=list()
          L$count=cbind(F_D,S_D)
          L$condt=c(rep(1,dim(F_D)[2]),rep(2,dim(S_D)[2]))
          D_E_2=pred_label_founder(L,"null_noise_2_20","gene",k,ii,jj)
          result$gene=append(result$gene,list(D_E_2))
          ############################## next experiment on samples
          F_D=data[,F_G]                
          S_D=data[,S_G]
          F_D=data[,F_G]                
          S_D=data[,S_G]
          random_subgroup=sample(1:2,1)
          if(random_subgroup==1)
          {
            # F_sum=colSums(F_D)
            # F_med=median(F_sum)
            # F_D=t((t(F_D)/F_sum)*F_med)
            print(sum(is.na(F_D)))
            for(i in 1:dim(F_D)[2])
            {
              zero_samples=sample(0,length(F_D[,i]),replace = TRUE)
              random_10_percent_index=sample(1:length(F_D[,i]),as.integer(length(F_D[,i])/5))
              random_10_percent_sample=rnorm(as.integer(length(F_D[,i])/5),mean=(mean(F_D[,i])+sd(F_D[,i])*ii),sd=sd(F_D[,i])*jj)
              zero_samples[random_10_percent_index]<-random_10_percent_sample
              F_D[,i]=F_D[,i]+zero_samples
            }
            # F_D=t((t(F_D)*F_sum)/F_med)
            F_D=apply(F_D,2,function(x) as.numeric(as.integer(x)))
            F_D[F_D<0]<-0
            print(sum(is.na(F_D)))
          }
          if(random_subgroup==2)
          {
            # S_sum=colSums(F_D)
            # S_med=median(S_sum)
            # S_D=t((t(S_D)/S_sum)*S_med)
            print(sum(is.na(S_D)))
            for(i in 1:dim(S_D)[2])
            {
              zero_samples=sample(0,length(S_D[,i]),replace = TRUE)
              random_10_percent_index=sample(1:length(S_D[,i]),as.integer(length(S_D[,i])/5))
              random_10_percent_sample=rnorm(as.integer(length(S_D[,i])/5),mean=(mean(S_D[,i])+sd(S_D[,i])*ii),sd=sd(S_D[,i])*jj)
              zero_samples[random_10_percent_index]<-random_10_percent_sample
              S_D[,i]=S_D[,i]+zero_samples
            }
            # S_D=t((t(S_D)*S_sum)/S_med)
            S_D=apply(S_D,2,function(x) as.numeric(as.integer(x)))
            S_D[S_D<0]<-0
            print(sum(is.na(S_D)))
          }
          L=list()
          L$count=cbind(F_D,S_D)
          L$condt=c(rep(1,dim(F_D)[2]),rep(2,dim(S_D)[2]))
          D_E_2=pred_label_founder(L,"null_noise_2_20","sample",k,ii,jj)
          result$sample=append(result$sample,list(D_E_2))
          saveRDS(result,paste("~/ROSeq/paper_work/latest_data/noise_jurkat/null_noise_2_20","_",k,".rds",sep=""))
        }
      }
  }
}
