file="~/ROSeq/paper_work/latest_data/noise_jurkat/null_noise_uniform_2_10_"
titles_names_1=c("MAST","BPSC","Wilcoxon","ROSeq","DESeq2","SCDE")
titles_names_2=c("MAST","BPSC","Wilc","ROSe","DEse","SCDE")
library(ggplot2)
n_noise=4
n_set=10
normd=c()
for (k in c(1:n_set)*25)
{
  print(k)
  for (j in 1:n_noise)
  {
    data=readRDS(paste(file,k,".rds",sep=""))
    normd=c(normd,colSums(as.matrix(data$gene[[j]])<.05)[c(1,3,5,7,9,11)])
  }
}
normd=unlist(normd)
normd=(normd-min(normd))/(max(normd)-min(normd))
plot=list()
for (k in c(1:6))
{
  print(k)
  # data=readRDS(paste("~/ROSeq/paper_work/latest_data/noise_jurkat/result_noise_noise",k*25,".rds",sep=""))
  data$gene=normd[as.integer(c((1:(length(normd)/6))-1)*6)+k]
  # print(data$gene)
  if(substr(names(data$gene)[1],1,4)==titles_names_2[k])
    title=titles_names_1[k]
  mat=matrix(data$gene,ncol=n_set,nrow=n_noise)
  colnames(mat)=c(1:n_set)*25
  rownames(mat)=paste("1, ",c(1:n_noise)*.1,"m",sep="")
  plot_data <- expand.grid(n_cells=colnames(mat),noise=rownames(mat))
  plot_data$Z <- c(t(mat))
  #plot_data$Z <- data$gene_2
  if(k<6)
    plot[[k]]=ggplot(plot_data, aes(n_cells, noise, fill=Z)) +
    geom_tile() + geom_text(aes(label = round(Z*100, 3))) +
    ggtitle(title) + theme_classic() + 
    theme(axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())+
    # scale_fill_gradient(low = "white", high = "red")
    scale_fill_continuous(low = "white", high = "red",limits=c(0, 1), breaks=seq(0,1,by=1))
  else
    plot[[k]]=ggplot(plot_data, aes(n_cells, noise, fill=Z)) +
    geom_tile() + geom_text(aes(label = round(Z*100, 3))) +
    labs(x = "# Cells") +
    ggtitle(title) + theme_classic() +
    theme(axis.title.y = element_blank())+
    # scale_fill_gradient(low = "white", high = "red")
    scale_fill_continuous(low = "white", high = "red",limits=c(0, 1), breaks=seq(0,1,by=1))
}
library(ggpubr)
a=ggarrange(plot[[1]],plot[[2]],plot[[3]],plot[[4]],plot[[5]],plot[[6]],
            ncol = 1,nrow=6,common.legend = TRUE)
pdf("~/ROSeq/paper_work/latest_plot/noise_jurkat/uniform_noise_combined_all_gene_2_10_pval.pdf",height = 10)
a
dev.off()
a
# ######################################## another way
file="~/ROSeq/paper_work/latest_data/noise_jurkat/null_noise_uniform_2_10_"
library(ggplot2)
n_noise=4
n_set=10
normd=c()
for (k in c(1:n_set)*25)
{
  print(k)
  for (j in 1:n_noise)
  {
    data=readRDS(paste(file,k,".rds",sep=""))
    normd=c(normd,colSums(as.matrix(data$sample[[j]])<.05)[c(1,3,5,7,9,11)])
  }
}
normd=unlist(normd)
normd=(normd-min(normd))/(max(normd)-min(normd))
plot=list()
for (k in c(1:6))
{
  print(k)
  data$sample=normd[as.integer(c((1:(length(normd)/6))-1)*6)+k]
  # print(data$sample)
  if(substr(names(data$sample)[1],1,4)==titles_names_2[k])
    title=titles_names_1[k]
  mat=matrix(data$sample,ncol=n_set,nrow=n_noise)
  colnames(mat)=c(1:n_set)*25
  rownames(mat)=paste("1, ",c(1:n_noise)*.1,"m",sep="")
  plot_data <- expand.grid(n_cells=colnames(mat),noise=rownames(mat))
  plot_data$Z <- c(t(mat))
  #plot_data$Z <- data$sample_2
  if(k<6)
    plot[[k]]=ggplot(plot_data, aes(n_cells, noise, fill=Z)) +
    geom_tile() + geom_text(aes(label = round(Z*100, 3))) +
    ggtitle(title) + theme_classic() + 
    theme(axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())+
    # scale_fill_gradient(low = "white", high = "red")
    scale_fill_continuous(low = "white", high = "red",limits=c(0, 1), breaks=seq(0,1,by=1))
  else
    plot[[k]]=ggplot(plot_data, aes(n_cells, noise, fill=Z)) +
    geom_tile() + geom_text(aes(label = round(Z*100, 3))) +
    labs(x = "# Cells") +
    ggtitle(title) + theme_classic() +
    theme(axis.title.y = element_blank())+
    # scale_fill_gradient(low = "white", high = "red")
    scale_fill_continuous(low = "white", high = "red",limits=c(0, 1), breaks=seq(0,1,by=1))
}
library(ggpubr)
a=ggarrange(plot[[1]],plot[[2]],plot[[3]],plot[[4]],plot[[5]],plot[[6]],
            ncol = 1,nrow=6,common.legend = TRUE)
pdf("~/ROSeq/paper_work/latest_plot/noise_jurkat/uniform_noise_combined_all_sample_2_10_pval.pdf",height = 10)
a
dev.off()
a
