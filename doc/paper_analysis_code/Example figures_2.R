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
library(DESeq2, lib.loc = "/usr/local/lib/R/site-library")
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
  



# ------------------------------------------------Trap dataset




# ------------------------------------------------Tung dataset
# molecules_1 <- read.table("~/ROSeq/Functions/singleCellSeq-gh-pages/data/molecules.txt", sep = "\t")
# molecules=molecules_1[-1,-1]
# bulk_counts=apply(molecules,2,as.numeric)
# rownames(bulk_counts)=molecules_1[,1][-1]
# colnames(bulk_counts)=as.matrix(molecules_1[1,][-1])
# anno <- read.table("~/ROSeq/Functions/singleCellSeq-gh-pages/data/annotation.txt", sep = "\t", header = TRUE)[,1]
# 
# L_Tung=list()
# L_Tung$NA19098_count=bulk_counts[,which(anno %in% "NA19098")]
# L_Tung$NA19101_count=bulk_counts[,which(anno %in% "NA19101")]
# L_Tung$NA19239_count=bulk_counts[,which(anno %in% "NA19239")]
# L_Tung$NA19098_group=rep(1,dim(L_Tung$NA19098_count)[2])
# L_Tung$NA19101_group=rep(2,dim(L_Tung$NA19101_count)[2])
# L_Tung$NA19239_group=rep(3,dim(L_Tung$NA19239_count)[2])
# L_Tung$NA19098_NA19101_count=cbind(L_Tung$NA19098_count,L_Tung$NA19101_count)
# L_Tung$NA19098_NA19101_group=c(L_Tung$NA19098_group,L_Tung$NA19101_group)
# L_Tung$NA19098_NA19239_count=cbind(L_Tung$NA19098_count,L_Tung$NA19239_count)
# L_Tung$NA19098_NA19239_group=c(L_Tung$NA19098_group,L_Tung$NA19239_group)
# L_Tung$NA19239_NA19101_count=cbind(L_Tung$NA19239_count,L_Tung$NA19101_count)
# L_Tung$NA19239_NA19101_group=c(L_Tung$NA19239_group,L_Tung$NA19101_group)
# L=list()
# L$count=L_Tung$NA19098_NA19101_count
# L$condt=L_Tung$NA19098_NA19101_group

load("~/ROSeq/paper_work/latest_data/pair/data2.RData")
bulk_pseudo_count=good_pseudo_count
name=rownames(bulk_pseudo_count)
bulk_pseudo_count=apply(bulk_pseudo_count, 2, as.numeric)
rownames(bulk_pseudo_count)<- name
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
L$count=L$count[,colSums(L$count> 0) > 2000]
L$count=gene_filter(L$count,min.count=2,min.cell=3)
colnames(L$count)=paste(colnames(L$count),1:dim(L$count)[2],sep="_")

genes=dim(L$count)[1]
iter=10
L$count=L$count[sample(1:dim(L$count)[1],genes),]

source("~/ROSeq/paper_work/latest_script/ROSeq_2.R")
R10=ROSeq(countData = limma::voom(TMMnormalization(L$count))$E, condition = L$condt, numCores = 1, new_step=0.05)
saveRDS(R10,"~/ROSeq/paper_work/latest_data/Example/R10.rds")
R10 <- readRDS("~/ROSeq/paper_work/latest_data/Example/R10.rds")

DE_genes=rownames(R10$stats)[R10$stats[,2]<0.05]
NDE_genes=rownames(R10$stats)[R10$stats[,2]>0.05]
n_clusters=8
n_times=10
ht=100
pheat=pheatmap::pheatmap(R10$stats[DE_genes,c("F_a","F_b","S_a","S_b" )],kmeans_k = n_clusters)
pheat$kmeans$cluster
DE_cluster_genes=c()
for(i in 1:n_clusters)
DE_cluster_genes=c(DE_cluster_genes,names(pheat$kmeans$cluster)[pheat$kmeans$cluster==i][c(1:n_times)])
pheat=pheatmap::pheatmap(R10$stats[NDE_genes,c("F_a","F_b","S_a","S_b" )],kmeans_k = n_clusters)
pheat$kmeans$cluster
NDE_cluster_genes=c()
for(i in 1:n_clusters)
NDE_cluster_genes=c(NDE_cluster_genes,names(pheat$kmeans$cluster)[pheat$kmeans$cluster==i][c(1:n_times)])
############################################ next figures combine both group and model for combine group debsity 
######################## for DE
new_plot=list()
k=1
for (i in c(match(DE_cluster_genes,rownames(R10$stats))))
{
    df=data.frame("Exp"=c(R10$plots_1[[i]]$ds,R10$plots_2[[i]]$ds),
                  "Group"=c(rep("G1",length(R10$plots_1[[i]]$ds)),rep("G2",length(R10$plots_2[[i]]$ds))))
    new_plot[[k]]=ggplot(df, aes(x=Exp,color=Group),size=.5) +
          geom_density() +
           ggtitle(rownames(R10$stats)[[i]]) + xlab("expression") + ylab("density") +
          theme_classic() + theme(
              plot.title = element_text(color="black", size=14, face="bold.italic"),
              #axis.title.x=element_blank(),
              #axis.text.x=element_blank(),
              #axis.ticks.x=element_blank(),
              #axis.line.x = element_blank(),
              #axis.title.y=element_blank(),
              #axis.text.y=element_blank(),
              #axis.ticks.y=element_blank(),
              #axis.line.y = element_blank(),
              legend.position="none")
    k=k+1
    df=data.frame("NF"=c(R10$plots_1[[i]]$normalized_read_count_sorted,R10$plots_1[[i]]$f),
                  "NFL"=factor(c(rep("Actual",length(R10$plots_1[[i]]$normalized_read_count_sorted)),
                                 rep("Model",length(R10$plots_1[[i]]$f)))),
                  "rank"=c(R10$plots_1[[i]]$rank,R10$plots_1[[i]]$rank))
    new_plot[[k]]=ggplot(df, aes(x=rank,y=NF,color=NFL)) + geom_line(size=.5) +
      ggtitle(rownames(R10$stats)[[i]])  + xlab("Rank of bin") + ylab("norm. freq.") +
      theme_classic() + theme(
        plot.title = element_text(color="black", size=14, face="bold.italic"),
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        #axis.line.x = element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        #axis.line.y = element_blank(),
        legend.position="none")
    k=k+1
    df=data.frame("NF"=c(R10$plots_2[[i]]$normalized_read_count_sorted,R10$plots_2[[i]]$f),
                  "NFL"=factor(c(rep("Actual",length(R10$plots_2[[i]]$normalized_read_count_sorted)),
                                 rep("Model",length(R10$plots_2[[i]]$f)))),
                  "rank"=c(R10$plots_2[[i]]$rank,R10$plots_2[[i]]$rank))
    new_plot[[k]]=ggplot(df, aes(x=rank,y=NF,color=NFL)) + geom_line(size=.5) +
      ggtitle(rownames(R10$stats)[[i]])  + xlab("Rank of bin") + ylab("norm. freq.") +
      theme_classic() + theme(
        plot.title = element_text(color="black", size=14, face="bold.italic"),
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        #axis.line.x = element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        #axis.line.y = element_blank(),
        legend.position="none")
    k=k+1
}
pdf(paste("~/ROSeq/paper_work/latest_plot/Example/combine_DE_exp_model_combine_density_a_b_",iter,".pdf",sep=""),
    height = ht)
print(cowplot::plot_grid(plotlist = new_plot,nrow=n_clusters*n_times,ncol=3))
dev.off()
######################## for non NDE
new_plot=list()
k=1
for (i in c(match(NDE_cluster_genes,rownames(R10$stats))))
{
    df=data.frame("Exp"=c(R10$plots_1[[i]]$ds,R10$plots_2[[i]]$ds),
                  "Group"=c(rep("Group 1",length(R10$plots_1[[i]]$ds)),rep("Group 2",length(R10$plots_2[[i]]$ds))))
    new_plot[[k]]=ggplot(df, aes(x=Exp,color=Group),size=.5) +
      geom_density() +
      ggtitle(rownames(R10$stats)[[i]]) + xlab("expression") + ylab("density") +
      theme_classic() + theme(
        plot.title = element_text(color="black", size=14, face="bold.italic"),
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        #axis.line.x = element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        #axis.line.y = element_blank(),
        legend.position="none")
    k=k+1
    df=data.frame("NF"=c(R10$plots_1[[i]]$normalized_read_count_sorted,R10$plots_1[[i]]$f),
                  "NFL"=factor(c(rep("Actual",length(R10$plots_1[[i]]$normalized_read_count_sorted)),
                                 rep("Model",length(R10$plots_1[[i]]$f)))),
                  "rank"=c(R10$plots_1[[i]]$rank,R10$plots_1[[i]]$rank))
    new_plot[[k]]=ggplot(df, aes(x=rank,y=NF,color=NFL)) + geom_line(size=.5) +
      ggtitle(rownames(R10$stats)[[i]])  + xlab("Rank of bin") + ylab("norm. freq.") +
      theme_classic() + theme(
        plot.title = element_text(color="black", size=14, face="bold.italic"),
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        #axis.line.x = element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        #axis.line.y = element_blank(),
        legend.position="none")
    k=k+1
    df=data.frame("NF"=c(R10$plots_2[[i]]$normalized_read_count_sorted,R10$plots_2[[i]]$f),
                  "NFL"=factor(c(rep("Actual",length(R10$plots_2[[i]]$normalized_read_count_sorted)),
                                 rep("Model",length(R10$plots_2[[i]]$f)))),
                  "rank"=c(R10$plots_2[[i]]$rank,R10$plots_2[[i]]$rank))
    new_plot[[k]]=ggplot(df, aes(x=rank,y=NF,color=NFL)) + geom_line(size=.5) +
      ggtitle(rownames(R10$stats)[[i]])  + xlab("Rank of bin") + ylab("norm. freq.") +
      theme_classic() + theme(
        plot.title = element_text(color="black", size=14, face="bold.italic"),
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        #axis.line.x = element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        #axis.line.y = element_blank(),
        legend.position="none")
    k=k+1
}
pdf(paste("~/ROSeq/paper_work/latest_plot/Example/combine_NDE_exp_model_combine_density_a_b_",iter,".pdf",sep=""),
    height=ht)
print(cowplot::plot_grid(plotlist = new_plot,nrow=n_clusters*n_times,ncol=3))
dev.off()
###########################################################################################
NDE_cluster_genes=c(
  #"AGTR1","LINC00152","EIF1AX","ZNF207","HERC1","TPRG1L","HCFC1R1","POLR2K",
  "ZNF207","POLR2F","CLCN4","RSL1D1","METTL7A")
DE_cluster_genes=c(
  #"STMN3","COPS6","GPX8","ANAPC1","SLC20A1","DKK1","CHD1L","LGR4",
  "ERCC6L","HSPA13","TRMT6","STMN3","SERPINE1")
n_clusters=5
n_times=1
iter=11
ht=7
######################## for non NDE and DE combine
new_plot_1=list()
k=1
for (i in c(match(NDE_cluster_genes,rownames(R10$stats))))
{
  df=data.frame("Exp"=c(R10$plots_1[[i]]$ds,R10$plots_2[[i]]$ds),
                "Group"=c(rep("Gr1",length(R10$plots_1[[i]]$ds)),rep("Gr2",length(R10$plots_2[[i]]$ds))))
  new_plot[[k]]=ggplot(df, aes(x=Exp,color=Group),size=.5) +
    geom_density() +
    ggtitle(rownames(R10$stats)[[i]]) + xlab("expression") + ylab("density") +
    theme_classic() + theme(
      plot.title = element_text(color="black", size=14, face="bold.italic"),
      legend.position="none")
  k=k+1
  df=data.frame("NF"=c(R10$plots_1[[i]]$f,R10$plots_2[[i]]$f),
                "NFL"=factor(c(rep("Gr1",length(R10$plots_1[[i]]$f)),
                               rep("Gr2",length(R10$plots_2[[i]]$f)))),
                "rank"=c(R10$plots_1[[i]]$rank,R10$plots_2[[i]]$rank))
  new_plot[[k]]=ggplot(df, aes(x=rank,y=NF,color=NFL)) + geom_line(size=.5) +
    ggtitle(rownames(R10$stats)[[i]])  + xlab("Rank of bin") + ylab("norm. freq.") +
    theme_classic() + theme(
      plot.title = element_text(color="black", size=14, face="bold.italic")
      ,legend.position="none"
      )
  k=k+1
}
A1=cowplot::plot_grid(plotlist = new_plot,nrow=n_clusters*n_times,ncol=2)

new_plot_1=list()
k=1
for (i in c(match(DE_cluster_genes,rownames(R10$stats))))
{
  df=data.frame("Exp"=c(R10$plots_1[[i]]$ds,R10$plots_2[[i]]$ds),
                "Group"=c(rep("Gr1",length(R10$plots_1[[i]]$ds)),rep("Gr2",length(R10$plots_2[[i]]$ds))))
  new_plot[[k]]=ggplot(df, aes(x=Exp,color=Group),size=.5) +
    geom_density() +
    ggtitle(rownames(R10$stats)[[i]]) + xlab("expression") + ylab("density") +
    theme_classic() + theme(
      plot.title = element_text(color="black", size=14, face="bold.italic"),
      legend.position="none")
  k=k+1
  df=data.frame("NF"=c(R10$plots_1[[i]]$f,R10$plots_2[[i]]$f),
                "NFL"=factor(c(rep("Gr1",length(R10$plots_1[[i]]$f)),
                               rep("Gr2",length(R10$plots_2[[i]]$f)))),
                "rank"=c(R10$plots_1[[i]]$rank,R10$plots_2[[i]]$rank))
  new_plot[[k]]=ggplot(df, aes(x=rank,y=NF,color=NFL)) + geom_line(size=.5) +
    ggtitle(rownames(R10$stats)[[i]])  + xlab("Rank of bin") + ylab("norm. freq.") +
    theme_classic() + theme(
      plot.title = element_text(color="black", size=14, face="bold.italic")
      ,legend.position="none"
      )
  k=k+1
}
A2=cowplot::plot_grid(plotlist = new_plot,nrow=n_clusters*n_times,ncol=2)

pdf(paste("~/ROSeq/paper_work/latest_plot/Example/combine_NDE_DE_exp_model_combine_density_a_b_",iter,".pdf",sep=""),
    height=ht)
print(cowplot::plot_grid(plotlist = list(A1,A2),nrow=1,ncol=2
                         #,labels = c("non-differential genes","differential genes")
                         ))
dev.off()
#################################################33
library(ggplot2)
df=data.frame(R10$stats)
gg_hist<-ggplot(df, aes(x=F_R2)) + geom_histogram(color="red",fill="red") + 
  ggtitle("")  + xlab("R2") + ylab("Gene Frequency") +
  theme_classic() 
print(gg_hist)
pdf("~/ROSeq/paper_work/latest_plot/Example/G1_R2.pdf")
print(gg_hist)
dev.off()
gg_hist<-ggplot(df, aes(x=S_R2)) + geom_histogram(color="green",fill="green") +
  ggtitle("")  + xlab("R2") + ylab("Gene Frequency") +
  theme_classic() 
print(gg_hist)
pdf("~/ROSeq/paper_work/latest_plot/Example/G2_R2.pdf")
print(gg_hist)
dev.off()
##########
df=data.frame("R2"=c(R10$stats[,5],R10$stats[,10]),"Group"=c(rep("Gr1",length(R10$stats[,5])),rep("Gr2",length(R10$stats[,5]))))
gg_hist<-ggplot(df, aes(x=R2, color=Group)) + geom_density() + 
  ggtitle("")  + xlab("R2") + ylab("Gene Frequency") +
  theme_classic() 
print(gg_hist)
pdf("~/ROSeq/paper_work/latest_plot/Example/G1_G2_combine_R2.pdf")
print(gg_hist)
dev.off()
##########################################
i=match("COPS6",rownames(R10$stats))
R10$stats[i,c(5,6,9)]
df=data.frame("NF"=c(R10$plots_1[[i]]$normalized_read_count_sorted,R10$plots_1[[i]]$f),
              "NFL"=factor(c(rep("Actual",length(R10$plots_1[[i]]$normalized_read_count_sorted)),
                      rep("Model",length(R10$plots_1[[i]]$f)))),
              "rank"=c(R10$plots_1[[i]]$rank,R10$plots_1[[i]]$rank))
gg_dot<-ggplot(df, aes(x=rank,y=NF,color=NFL)) + geom_line(size=2) + 
  ggtitle(paste("Gene name: ",rownames(R10$stats)[[i]],sep=""))  + 
  xlab("Rank of bin") + ylab("Normalised frequency") +
  theme_classic()
print(gg_dot)
pdf("~/ROSeq/paper_work/latest_plot/Example/model_1.pdf")
print(gg_dot)
dev.off()
df=data.frame("NF"=c(R10$plots_2[[i]]$normalized_read_count_sorted,R10$plots_2[[i]]$f),
              "NFL"=factor(c(rep("Actual",length(R10$plots_2[[i]]$normalized_read_count_sorted)),
                             rep("Model",length(R10$plots_2[[i]]$f)))),
              "rank"=c(R10$plots_2[[i]]$rank,R10$plots_2[[i]]$rank))
gg_dot<-ggplot(df, aes(x=rank,y=NF,color=NFL)) + geom_line(size=2) + 
  ggtitle(paste("Gene name: ",rownames(R10$stats)[[i]],sep=""))  + 
  xlab("Rank of bin") + ylab("Normalised frequency") +
  theme_classic()
print(gg_dot)
pdf("~/ROSeq/paper_work/latest_plot/Example/model_2.pdf")
print(gg_dot)
dev.off()
