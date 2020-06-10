cm <- function(ground_labels, predicted_labels)
{
  ground_labels = factor(as.matrix(ground_labels), levels=c(0,1))
  predicted_labels = factor(as.matrix(predicted_labels),levels=c(0,1))
  tp = 0
  tn = 0
  fn = 0
  fp = 0
  for(k in 1:length(ground_labels)){
    if(ground_labels[k]==predicted_labels[k] && ground_labels[k]==0)
      tp <- tp + 1
    else if(ground_labels[k]==predicted_labels[k] && ground_labels[k]==1)
      tn <- tn + 1
    else if(ground_labels[k]!=predicted_labels[k] && ground_labels[k]==0)
      fp <- fp + 1
    else
      fn <- fn + 1
  }
  recall = tp/(tp+fn)
  precision = tp/(tp+fp)
  total = tp+tn+fp+fn
  acc = (tp+tn)/total
  f1_score = 2/((1/precision)+(1/recall))
  den = sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc = (tp*tn-fp*fn)/den
  po = acc
  pe = ((tp+fp)*(tp+fn)+(fn+tn)*(tn+fp))/(total*total)
  kappa = (po-pe)/(1-pe)
  
  stats=c(mcc,f1_score,precision,recall,kappa,acc)
  names(stats)<-c("mcc","f1_score","precision","recall","kappa","acc")
  return(stats )
}
L_Trap_single <- readRDS("~/ROSeq/paper_work/latest_data/pair/L_Trap_single.rds")
L_Tung_single <- readRDS("~/ROSeq/paper_work/latest_data/pair/L_Tung_single.rds")
data=list()
plot=list()
plot_seq=c(1,2,3,7,9,10)
plot_seq_1=c(1,2,6,8)
alpha=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
################################################################# L_Trap_single Data ###########################################
###################################### L_Trap_single$T0_T24 Data #########
library(gtools)
library(ggplot2)
library(plotROC)
library(pheatmap)
library(EnhancedVolcano)
library(biomaRt)
known_pval=c()
known_padj=c()
pred_pval=c()
pred_padj=c()
method_pval=c()
method_padj=c()
j=0
for(i in plot_seq)
{
  j=j+1
  Truth_trap <- readRDS("~/ROSeq/Functions/02_Trapnell/Truth_trap.rds")
  common=intersect(rownames(L_Trap_single$T0_T24_pred_label),names(Truth_trap))
  t_index=match(common,names(Truth_trap))
  p_index=match(common,rownames(L_Trap_single$T0_T24_pred_label))
  known_pval=c(known_pval,as.numeric(Truth_trap[t_index]))
  known_padj=c(known_padj,as.numeric(Truth_trap[t_index]))
  pred_pval=c(pred_pval,L_Trap_single$T0_T24_pred_label[p_index,i*2-1])
  pred_padj=c(pred_padj,L_Trap_single$T0_T24_pred_label[p_index,i*2])
  method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_Trap_single$T0_T24_pred_label)[i*2-1],sep = "_"),dim(L_Trap_single$T0_T24_pred_label[p_index,])[1]))
  method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_Trap_single$T0_T24_pred_label)[i*2],sep = "_"),dim(L_Trap_single$T0_T24_pred_label[p_index,])[1]))
}
data$Trap_single_T0_T24=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)

## draw plots
basicplot_pval <- ggplot(data$Trap_single_T0_T24, aes(d = truth_pval, m = prediction_pval, color = model_pval)) + 
  geom_roc(n.cuts = 0) + 
  style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") 
## calculate auc
calc_auc(basicplot_pval)

basicplot_padj <- ggplot(data$Trap_single_T0_T24, aes(d = truth_padj, m = prediction_padj, color = model_padj)) + 
  geom_roc(n.cuts = 0) + 
  style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") 
## calculate auc
calc_auc(basicplot_padj)

known_pval=c()
known_padj=c()
pred_pval=c()
pred_padj=c()
method_pval=c()
method_padj=c()
j=0
for(i in plot_seq)
{
  j=j+1
  Truth_trap <- readRDS("~/ROSeq/Functions/02_Trapnell/Truth_trap.rds")
  common=intersect(rownames(L_Trap_single$T0_T24_pred_label),names(Truth_trap))
  t_index=match(common,names(Truth_trap))
  p_index=match(common,rownames(L_Trap_single$T0_T24_pred_label))
  known_pval=c(known_pval,as.numeric(Truth_trap[t_index]))
  known_padj=c(known_padj,as.numeric(Truth_trap[t_index]))
  pred_pval=c(pred_pval,L_Trap_single$T0_T24_pred_label[p_index,i*2-1])
  pred_padj=c(pred_padj,L_Trap_single$T0_T24_pred_label[p_index,i*2])
  method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_Trap_single$T0_T24_pred_label)[i*2-1],round(calc_auc(basicplot_pval)[j,3],2),sep = "_"),dim(L_Trap_single$T0_T24_pred_label[p_index,])[1]))
  method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_Trap_single$T0_T24_pred_label)[i*2],round(calc_auc(basicplot_padj)[j,3],2),sep = "_"),dim(L_Trap_single$T0_T24_pred_label[p_index,])[1]))
}
data$Trap_single_T0_T24=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)

############################# draw plots pval
basicplot_pval <- ggplot(data$Trap_single_T0_T24, aes(d = truth_pval, m = prediction_pval, color = model_pval)) + 
  geom_roc(n.cuts = 0) + style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") +
  scale_color_manual(values=c("red", "blue", "dark green","black","orange","magenta"))
## calculate auc
calc_auc(basicplot_pval)

plot$Trap_single_T0_T24_AUC <- basicplot_pval + 
  ggtitle("ROC") + theme_classic()

count_data=L_Trap_single$T0_T24_count[rownames(L_Trap_single$T0_T24_pred_label),]
log2FC=c()
for(i in 1:dim(count_data)[1])
{
  fav=mean(count_data[i,L_Trap_single$T0_T24_group==1])
  sav=mean(count_data[i,L_Trap_single$T0_T24_group==2])
  if(fav==0)
    fav=0.000001
  if(sav==0)
    sav=0.000001
  log2FC=c(log2FC,log2(fav/sav))
}
names(log2FC)<-rownames(count_data)

DE_genes<-rownames(L_Trap_single$T0_T24_pred_label)[L_Trap_single$T0_T24_pred_label$ROSeq_0_tmmvoom_padj<.05 & log2FC>6.6 ]
DE_genes=c(DE_genes,rownames(L_Trap_single$T0_T24_pred_label)[L_Trap_single$T0_T24_pred_label$ROSeq_0_tmmvoom_padj<.05 &  log2FC<(-23.6)])
DE_data<-L_Trap_single$T0_T24_count[DE_genes,]
DE_info=data.frame("log2FC"=log2FC,"pval"=L_Trap_single$T0_T24_pred_label$ROSeq_0_tmmvoom_pval)
plot$Trap_single_T0_T24_volcano <- EnhancedVolcano(DE_info,
                                                   lab = rownames(DE_info),
                                                   x = 'log2FC',
                                                   y = 'pval',
                                                   xlim = c(-10, 10),
                                                   title = "Trap_1_2",
                                                   pCutoff = 0.05,
                                                   FCcutoff = 2,
                                                   # col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                                                   colAlpha = 1)

col_ann=data.frame("group"=L_Trap_single$T0_T24_group)
col_ann[col_ann==1]<-"T0"
col_ann[col_ann==2]<-"T24"
rownames(col_ann)<-colnames(DE_data)

plot$Trap_single_T0_T24_pheatmap <- pheatmap::pheatmap(log2(DE_data+1),
                                                       annotation_col = col_ann,show_rownames = T,show_colnames = F,cluster_cols=F,cluster_rows = F,
                                                       main = "heatmap")
library(grid)
library(pROC)
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_new_1_pval.pdf")
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[1]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[1]], 
            smooth=TRUE,auc=TRUE),col="red")
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[2]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[2]], 
            smooth=TRUE),col="blue",add = TRUE)
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[3]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[3]], 
            smooth=TRUE),col="dark green",add = TRUE)
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[4]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[4]], 
            smooth=TRUE),col="black",add = TRUE)
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[5]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[5]], 
            smooth=TRUE),col="orange",add = TRUE)
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[6]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[6]], 
            smooth=TRUE),col="magenta",
        legend(.65, .3, legend=unique(data$Trap_single_T0_T24$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
               cex=1,text.font=2),add = TRUE)
dev.off()
library(png)
png("Rlogo.png")
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[1]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[1]], 
            smooth=TRUE,auc=TRUE),col="red")
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[2]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[2]], 
            smooth=TRUE),col="blue",add = TRUE)
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[3]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[3]], 
            smooth=TRUE),col="dark green",add = TRUE)
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[4]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[4]], 
            smooth=TRUE),col="black",add = TRUE)
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[5]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[5]], 
            smooth=TRUE),col="orange",add = TRUE)
g<-plot(roc(data$Trap_single_T0_T24$truth_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[6]], 
            data$Trap_single_T0_T24$prediction_pval[data$Trap_single_T0_T24$model_pval==unique(data$Trap_single_T0_T24$model_pval)[6]], 
            smooth=TRUE),col="magenta",
        legend(.65, .3, legend=unique(data$Trap_single_T0_T24$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
               cex=1,text.font=2),add = TRUE)
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Trap_single_T0_T24_AUC<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
plot$Trap_single_T0_T24_AUC
png("Rlogo.png")
plot$Trap_single_T0_T24_volcano
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Trap_single_T0_T24_volcano<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
png("Rlogo.png")
plot$Trap_single_T0_T24_pheatmap
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Trap_single_T0_T24_pheatmap<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
################################################################# L_Tung_single Data ###########################################
###################################### L_Tung_single$NA19098_NA19101 Data #########
library(gtools)
library(ggplot2)
library(plotROC)
library(pheatmap)
known_pval=c()
known_padj=c()
pred_pval=c()
pred_padj=c()
method_pval=c()
method_padj=c()
j=0
for(i in plot_seq)
{
  j=j+1
  Truth_Tung <- readRDS("~/ROSeq/Functions/03_Tung_12/Truth_trap_NA19098_NA19101.rds")
  common=intersect(rownames(L_Tung_single$NA19098_NA19101_pred_label),names(Truth_Tung))
  t_index=match(common,names(Truth_Tung))
  p_index=match(common,rownames(L_Tung_single$NA19098_NA19101_pred_label))
  known_pval=c(known_pval,as.numeric(Truth_Tung[t_index]))
  known_padj=c(known_padj,as.numeric(Truth_Tung[t_index]))
  pred_pval=c(pred_pval,L_Tung_single$NA19098_NA19101_pred_label[p_index,i*2-1])
  pred_padj=c(pred_padj,L_Tung_single$NA19098_NA19101_pred_label[p_index,i*2])
  method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_Tung_single$NA19098_NA19101_pred_label)[i*2-1],sep = "_"),dim(L_Tung_single$NA19098_NA19101_pred_label[p_index,])[1]))
  method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_Tung_single$NA19098_NA19101_pred_label)[i*2],sep = "_"),dim(L_Tung_single$NA19098_NA19101_pred_label[p_index,])[1]))
}
data$Tung_single_NA19098_NA19101=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)

## draw plots
basicplot_pval <- ggplot(data$Tung_single_NA19098_NA19101, aes(d = truth_pval, m = prediction_pval, color = model_pval)) + 
  geom_roc(n.cuts = 0) + 
  style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") 
## calculate auc
calc_auc(basicplot_pval)

basicplot_padj <- ggplot(data$Tung_single_NA19098_NA19101, aes(d = truth_padj, m = prediction_padj, color = model_padj)) + 
  geom_roc(n.cuts = 0) + 
  style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") 
## calculate auc
calc_auc(basicplot_padj)

known_pval=c()
known_padj=c()
pred_pval=c()
pred_padj=c()
method_pval=c()
method_padj=c()
j=0
for(i in plot_seq)
{
  j=j+1
  Truth_Tung <- readRDS("~/ROSeq/Functions/03_Tung_12/Truth_trap_NA19098_NA19101.rds")
  common=intersect(rownames(L_Tung_single$NA19098_NA19101_pred_label),names(Truth_Tung))
  t_index=match(common,names(Truth_Tung))
  p_index=match(common,rownames(L_Tung_single$NA19098_NA19101_pred_label))
  known_pval=c(known_pval,as.numeric(Truth_Tung[t_index]))
  known_padj=c(known_padj,as.numeric(Truth_Tung[t_index]))
  pred_pval=c(pred_pval,L_Tung_single$NA19098_NA19101_pred_label[p_index,i*2-1])
  pred_padj=c(pred_padj,L_Tung_single$NA19098_NA19101_pred_label[p_index,i*2])
  method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_Tung_single$NA19098_NA19101_pred_label)[i*2-1],round(calc_auc(basicplot_pval)[j,3],2),sep = "_"),dim(L_Tung_single$NA19098_NA19101_pred_label[p_index,])[1]))
  method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_Tung_single$NA19098_NA19101_pred_label)[i*2],round(calc_auc(basicplot_padj)[j,3],2),sep = "_"),dim(L_Tung_single$NA19098_NA19101_pred_label[p_index,])[1]))
}
data$Tung_single_NA19098_NA19101=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)

############################# draw plots pval
basicplot_pval <- ggplot(data$Tung_single_NA19098_NA19101, aes(d = truth_pval, m = prediction_pval, color = model_pval)) + 
  geom_roc(n.cuts = 0) + style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") +
  scale_color_manual(values=c("red", "blue", "dark green","black","orange","magenta"))
## calculate auc
calc_auc(basicplot_pval)

plot$Tung_single_NA19098_NA19101_AUC <- basicplot_pval + 
  ggtitle("ROC") + theme_classic()

count_data=L_Tung_single$NA19098_NA19101_count[rownames(L_Tung_single$NA19098_NA19101_pred_label),]
log2FC=c()
for(i in 1:dim(count_data)[1])
{
  fav=mean(count_data[i,L_Tung_single$NA19098_NA19101_group==1])
  sav=mean(count_data[i,L_Tung_single$NA19098_NA19101_group==2])
  if(fav==0)
    fav=0.000001
  if(sav==0)
    sav=0.000001
  log2FC=c(log2FC,log2(fav/sav))
}
names(log2FC)<-rownames(count_data)

DE_genes<-rownames(L_Tung_single$NA19098_NA19101_pred_label)[L_Tung_single$NA19098_NA19101_pred_label$ROSeq_0_tmmvoom_padj<.05 & log2FC>2.8 ]
DE_genes=c(DE_genes,rownames(L_Tung_single$NA19098_NA19101_pred_label)[L_Tung_single$NA19098_NA19101_pred_label$ROSeq_0_tmmvoom_padj<.05 &  log2FC<(-4.5)])
DE_data<-L_Tung_single$NA19098_NA19101_count[DE_genes,]
DE_info=data.frame("log2FC"=log2FC,"pval"=L_Tung_single$NA19098_NA19101_pred_label$ROSeq_0_tmmvoom_pval)
plot$Tung_single_NA19098_NA19101_volcano <- EnhancedVolcano(DE_info,
                                                            lab = rownames(DE_info),
                                                            x = 'log2FC',
                                                            y = 'pval',
                                                            xlim = c(-10, 10),
                                                            title = "Tung_1_2",
                                                            pCutoff = 0.05,
                                                            FCcutoff = 2,
                                                            # col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                                                            colAlpha = 1)

col_ann=data.frame("group"=L_Tung_single$NA19098_NA19101_group)
col_ann[col_ann==1]<-"NA19098"
col_ann[col_ann==2]<-"NA19101"
rownames(col_ann)<-colnames(DE_data)
plot$Tung_single_NA19098_NA19101_pheatmap <- pheatmap::pheatmap(log2(DE_data+1),
                                                                annotation_col = col_ann,show_rownames = T,show_colnames = F,cluster_cols=F,cluster_rows = F, 
                                                                main = "heatmap")
library(grid)
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_new_2_pval.pdf")
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[1]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[1]], 
            smooth=TRUE,auc=TRUE),col="red")
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[2]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[2]], 
            smooth=TRUE),col="blue",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[3]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[3]], 
            smooth=TRUE),col="dark green",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[4]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[4]], 
            smooth=TRUE),col="black",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[5]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[5]], 
            smooth=TRUE),col="orange",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[6]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[6]], 
            smooth=TRUE),col="magenta",
        legend(.65, .3, legend=unique(data$Tung_single_NA19098_NA19101$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
               cex=1,text.font=2),add = TRUE)
dev.off()
library(png)
png("Rlogo.png")
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[1]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[1]], 
            smooth=TRUE,auc=TRUE),col="red")
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[2]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[2]], 
            smooth=TRUE),col="blue",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[3]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[3]], 
            smooth=TRUE),col="dark green",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[4]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[4]], 
            smooth=TRUE),col="black",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[5]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[5]], 
            smooth=TRUE),col="orange",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19101$truth_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[6]], 
            data$Tung_single_NA19098_NA19101$prediction_pval[data$Tung_single_NA19098_NA19101$model_pval==unique(data$Tung_single_NA19098_NA19101$model_pval)[6]], 
            smooth=TRUE),col="magenta",
        legend(.65, .3, legend=unique(data$Tung_single_NA19098_NA19101$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
               cex=1,text.font=2),add = TRUE)
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Tung_single_NA19098_NA19101_AUC<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
plot$Tung_single_NA19098_NA19101_AUC
png("Rlogo.png")
plot$Tung_single_NA19098_NA19101_volcano
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Tung_single_NA19098_NA19101_volcano<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
png("Rlogo.png")
plot$Tung_single_NA19098_NA19101_pheatmap
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Tung_single_NA19098_NA19101_pheatmap<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())

###################################### L_Tung_single$NA19098_NA19239 Data #########
library(gtools)
library(ggplot2)
library(plotROC)
library(pheatmap)
known_pval=c()
known_padj=c()
pred_pval=c()
pred_padj=c()
method_pval=c()
method_padj=c()
j=0
for(i in plot_seq)
{
  j=j+1
  Truth_Tung <- readRDS("~/ROSeq/Functions/05_Tung_31/Truth_trap_NA19098_NA19239.rds")
  common=intersect(rownames(L_Tung_single$NA19098_NA19239_pred_label),names(Truth_Tung))
  t_index=match(common,names(Truth_Tung))
  p_index=match(common,rownames(L_Tung_single$NA19098_NA19239_pred_label))
  known_pval=c(known_pval,as.numeric(Truth_Tung[t_index]))
  known_padj=c(known_padj,as.numeric(Truth_Tung[t_index]))
  pred_pval=c(pred_pval,L_Tung_single$NA19098_NA19239_pred_label[p_index,i*2-1])
  pred_padj=c(pred_padj,L_Tung_single$NA19098_NA19239_pred_label[p_index,i*2])
  method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_Tung_single$NA19098_NA19239_pred_label)[i*2-1],sep = "_"),dim(L_Tung_single$NA19098_NA19239_pred_label[p_index,])[1]))
  method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_Tung_single$NA19098_NA19239_pred_label)[i*2],sep = "_"),dim(L_Tung_single$NA19098_NA19239_pred_label[p_index,])[1]))
}
data$Tung_single_NA19098_NA19239=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)

## draw plots
basicplot_pval <- ggplot(data$Tung_single_NA19098_NA19239, aes(d = truth_pval, m = prediction_pval, color = model_pval)) + 
  geom_roc(n.cuts = 0) + 
  style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") 
## calculate auc
calc_auc(basicplot_pval)

known_pval=c()
known_padj=c()
pred_pval=c()
pred_padj=c()
method_pval=c()
method_padj=c()
j=0
for(i in plot_seq)
{
  j=j+1
  Truth_Tung <- readRDS("~/ROSeq/Functions/05_Tung_31/Truth_trap_NA19098_NA19239.rds")
  common=intersect(rownames(L_Tung_single$NA19098_NA19239_pred_label),names(Truth_Tung))
  t_index=match(common,names(Truth_Tung))
  p_index=match(common,rownames(L_Tung_single$NA19098_NA19239_pred_label))
  known_pval=c(known_pval,as.numeric(Truth_Tung[t_index]))
  known_padj=c(known_padj,as.numeric(Truth_Tung[t_index]))
  pred_pval=c(pred_pval,L_Tung_single$NA19098_NA19239_pred_label[p_index,i*2-1])
  pred_padj=c(pred_padj,L_Tung_single$NA19098_NA19239_pred_label[p_index,i*2])
  method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_Tung_single$NA19098_NA19239_pred_label)[i*2-1],round(calc_auc(basicplot_pval)[j,3],2),sep = "_"),dim(L_Tung_single$NA19098_NA19239_pred_label[p_index,])[1]))
  method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_Tung_single$NA19098_NA19239_pred_label)[i*2],round(calc_auc(basicplot_padj)[j,3],2),sep = "_"),dim(L_Tung_single$NA19098_NA19239_pred_label[p_index,])[1]))
}
data$Tung_single_NA19098_NA19239=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)

############################# draw plots pval
basicplot_pval <- ggplot(data$Tung_single_NA19098_NA19239, aes(d = truth_pval, m = prediction_pval, color = model_pval)) + 
  geom_roc(n.cuts = 0) + style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") +
  scale_color_manual(values=c("red", "blue", "dark green","black","orange","magenta"))
## calculate auc
calc_auc(basicplot_pval)

basicplot_padj <- ggplot(data$Tung_single_NA19098_NA19239, aes(d = truth_padj, m = prediction_padj, color = model_padj)) + 
  geom_roc(n.cuts = 0) + 
  style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") 
## calculate auc
calc_auc(basicplot_padj)

plot$Tung_single_NA19098_NA19239_AUC <- basicplot_pval + 
  ggtitle("ROC") + theme_classic()

count_data=L_Tung_single$NA19098_NA19239_count[rownames(L_Tung_single$NA19098_NA19239_pred_label),]
log2FC=c()
for(i in 1:dim(count_data)[1])
{
  fav=mean(count_data[i,L_Tung_single$NA19098_NA19239_group==1])
  sav=mean(count_data[i,L_Tung_single$NA19098_NA19239_group==3])
  if(fav==0)
    fav=0.000001
  if(sav==0)
    sav=0.000001
  log2FC=c(log2FC,log2(fav/sav))
}
names(log2FC)<-rownames(count_data)

DE_genes<-rownames(L_Tung_single$NA19098_NA19239_pred_label)[L_Tung_single$NA19098_NA19239_pred_label$ROSeq_0_tmmvoom_padj<.05 & log2FC>5.322 ]
DE_genes=c(DE_genes,rownames(L_Tung_single$NA19098_NA19239_pred_label)[L_Tung_single$NA19098_NA19239_pred_label$ROSeq_0_tmmvoom_padj<.05 &  log2FC<(-4.15)])
DE_data<-L_Tung_single$NA19098_NA19239_count[DE_genes,]
DE_info=data.frame("log2FC"=log2FC,"pval"=L_Tung_single$NA19098_NA19239_pred_label$ROSeq_0_tmmvoom_pval)
plot$Tung_single_NA19098_NA19239_volcano <- EnhancedVolcano(DE_info,
                                                            lab = rownames(DE_info),
                                                            x = 'log2FC',
                                                            y = 'pval',
                                                            xlim = c(-10, 10),
                                                            title = "Tung_1_3",
                                                            pCutoff = 0.05,
                                                            FCcutoff = 2,
                                                            # col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                                                            colAlpha = 1)

col_ann=data.frame("group"=L_Tung_single$NA19098_NA19239_group)
col_ann[col_ann==1]<-"NA19098"
col_ann[col_ann==3]<-"NA19239"
rownames(col_ann)<-colnames(DE_data)
plot$Tung_single_NA19098_NA19239_pheatmap <- pheatmap::pheatmap(log2(DE_data+1),
                                                                annotation_col = col_ann,show_rownames = T,show_colnames = F,cluster_cols=F,cluster_rows = F, 
                                                                main = "heatmap")
library(grid)
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_new_3_pval.pdf")
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[1]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[1]], 
            smooth=TRUE,auc=TRUE),col="red")
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[2]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[2]], 
            smooth=TRUE),col="blue",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[3]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[3]], 
            smooth=TRUE),col="dark green",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[4]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[4]], 
            smooth=TRUE),col="black",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[5]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[5]], 
            smooth=TRUE),col="orange",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[6]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[6]], 
            smooth=TRUE),col="magenta",
        legend(.65, .3, legend=unique(data$Tung_single_NA19098_NA19239$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
               cex=1,text.font=2),add = TRUE)
dev.off()
library(png)
png("Rlogo.png")
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[1]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[1]], 
            smooth=TRUE,auc=TRUE),col="red")
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[2]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[2]], 
            smooth=TRUE),col="blue",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[3]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[3]], 
            smooth=TRUE),col="dark green",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[4]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[4]], 
            smooth=TRUE),col="black",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[5]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[5]], 
            smooth=TRUE),col="orange",add = TRUE)
g<-plot(roc(data$Tung_single_NA19098_NA19239$truth_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[6]], 
            data$Tung_single_NA19098_NA19239$prediction_pval[data$Tung_single_NA19098_NA19239$model_pval==unique(data$Tung_single_NA19098_NA19239$model_pval)[6]], 
            smooth=TRUE),col="magenta",
        legend(.65, .3, legend=unique(data$Tung_single_NA19098_NA19239$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
               cex=1,text.font=2),add = TRUE)
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Tung_single_NA19098_NA19239_AUC<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
plot$Tung_single_NA19098_NA19239_AUC
png("Rlogo.png")
plot$Tung_single_NA19098_NA19239_volcano
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Tung_single_NA19098_NA19239_volcano<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
png("Rlogo.png")
plot$Tung_single_NA19098_NA19239_pheatmap
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Tung_single_NA19098_NA19239_pheatmap<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
###################################### L_Tung_single$NA19239_NA19101 Data #########
library(gtools)
library(ggplot2)
library(plotROC)
library(NMF)
library(pheatmap)
known_pval=c()
known_padj=c()
pred_pval=c()
pred_padj=c()
method_pval=c()
method_padj=c()
j=0
for(i in plot_seq)
{
  j=j+1
  Truth_Tung <- readRDS("~/ROSeq/Functions/04_Tung_23/Truth_trap_NA19101_NA19239.rds")
  common=intersect(rownames(L_Tung_single$NA19239_NA19101_pred_label),names(Truth_Tung))
  t_index=match(common,names(Truth_Tung))
  p_index=match(common,rownames(L_Tung_single$NA19239_NA19101_pred_label))
  known_pval=c(known_pval,as.numeric(Truth_Tung[t_index]))
  known_padj=c(known_padj,as.numeric(Truth_Tung[t_index]))
  pred_pval=c(pred_pval,L_Tung_single$NA19239_NA19101_pred_label[p_index,i*2-1])
  pred_padj=c(pred_padj,L_Tung_single$NA19239_NA19101_pred_label[p_index,i*2])
  method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_Tung_single$NA19239_NA19101_pred_label)[i*2-1],sep = "_"),dim(L_Tung_single$NA19239_NA19101_pred_label[p_index,])[1]))
  method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_Tung_single$NA19239_NA19101_pred_label)[i*2],sep = "_"),dim(L_Tung_single$NA19239_NA19101_pred_label[p_index,])[1]))
}
data$Tung_single_NA19239_NA19101=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)

## draw plots
basicplot_pval <- ggplot(data$Tung_single_NA19239_NA19101, aes(d = truth_pval, m = prediction_pval, color = model_pval)) + 
  geom_roc(n.cuts = 0) + 
  style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") 
## calculate auc
calc_auc(basicplot_pval)

basicplot_padj <- ggplot(data$Tung_single_NA19239_NA19101, aes(d = truth_padj, m = prediction_padj, color = model_padj)) + 
  geom_roc(n.cuts = 0) + 
  style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") 
## calculate auc
calc_auc(basicplot_padj)

known_pval=c()
known_padj=c()
pred_pval=c()
pred_padj=c()
method_pval=c()
method_padj=c()
j=0
for(i in plot_seq)
{
  j=j+1
  Truth_Tung <- readRDS("~/ROSeq/Functions/04_Tung_23/Truth_trap_NA19101_NA19239.rds")
  common=intersect(rownames(L_Tung_single$NA19239_NA19101_pred_label),names(Truth_Tung))
  t_index=match(common,names(Truth_Tung))
  p_index=match(common,rownames(L_Tung_single$NA19239_NA19101_pred_label))
  known_pval=c(known_pval,as.numeric(Truth_Tung[t_index]))
  known_padj=c(known_padj,as.numeric(Truth_Tung[t_index]))
  pred_pval=c(pred_pval,L_Tung_single$NA19239_NA19101_pred_label[p_index,i*2-1])
  pred_padj=c(pred_padj,L_Tung_single$NA19239_NA19101_pred_label[p_index,i*2])
  method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_Tung_single$NA19239_NA19101_pred_label)[i*2-1],round(calc_auc(basicplot_pval)[j,3],2),sep = "_"),dim(L_Tung_single$NA19239_NA19101_pred_label[p_index,])[1]))
  method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_Tung_single$NA19239_NA19101_pred_label)[i*2],round(calc_auc(basicplot_padj)[j,3],2),sep = "_"),dim(L_Tung_single$NA19239_NA19101_pred_label[p_index,])[1]))
}
data$Tung_single_NA19239_NA19101=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)

############################# draw plots pval
basicplot_pval <- ggplot(data$Tung_single_NA19239_NA19101, aes(d = truth_pval, m = prediction_pval, color = model_pval)) + 
  geom_roc(n.cuts = 0) + style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") +
  scale_color_manual(values=c("red", "blue", "dark green","black","orange","magenta"))
## calculate auc
calc_auc(basicplot_pval)

plot$Tung_single_NA19239_NA19101_AUC <- basicplot_pval + 
  ggtitle("ROC") + theme_classic()

count_data=L_Tung_single$NA19239_NA19101_count[rownames(L_Tung_single$NA19239_NA19101_pred_label),]
log2FC=c()
for(i in 1:dim(count_data)[1])
{
  fav=mean(count_data[i,L_Tung_single$NA19239_NA19101_group==3])
  sav=mean(count_data[i,L_Tung_single$NA19239_NA19101_group==2])
  if(fav==0)
    fav=0.000001
  if(sav==0)
    sav=0.000001
  log2FC=c(log2FC,log2(fav/sav))
}
names(log2FC)<-rownames(count_data)

DE_genes<-rownames(L_Tung_single$NA19239_NA19101_pred_label)[L_Tung_single$NA19239_NA19101_pred_label$ROSeq_0_tmmvoom_padj<.05 & log2FC>1.9 ]
DE_genes=c(DE_genes,rownames(L_Tung_single$NA19239_NA19101_pred_label)[L_Tung_single$NA19239_NA19101_pred_label$ROSeq_0_tmmvoom_padj<.05 &  log2FC<(-5.1)])
DE_data<-L_Tung_single$NA19239_NA19101_count[DE_genes,]
DE_info=data.frame("log2FC"=log2FC,"pval"=L_Tung_single$NA19239_NA19101_pred_label$ROSeq_0_tmmvoom_pval)
plot$Tung_single_NA19239_NA19101_volcano <- EnhancedVolcano(DE_info,
                                                            lab = rownames(DE_info),
                                                            x = 'log2FC',
                                                            y = 'pval',
                                                            xlim = c(-10, 10),
                                                            title = "Tung_3_2",
                                                            pCutoff = 0.05,
                                                            FCcutoff = 2,
                                                            # col=c('grey75', 'grey50', 'grey25', 'red','blue'),
                                                            colAlpha = 1)


col_ann=data.frame("group"=L_Tung_single$NA19239_NA19101_group)
col_ann[col_ann==3]<-"NA19239"
col_ann[col_ann==2]<-"NA19101"
rownames(col_ann)<-colnames(DE_data)
plot$Tung_single_NA19239_NA19101_pheatmap <- pheatmap::pheatmap(log2(DE_data+1),
                                                                annotation_col = col_ann,show_rownames = T,show_colnames = F,cluster_cols=F,cluster_rows = F, 
                                                                main = "heatmap")
library(grid)
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_new_4_pval.pdf")
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[1]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[1]], 
            smooth=TRUE,auc=TRUE),col="red")
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[2]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[2]], 
            smooth=TRUE),col="blue",add = TRUE)
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[3]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[3]], 
            smooth=TRUE),col="dark green",add = TRUE)
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[4]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[4]], 
            smooth=TRUE),col="black",add = TRUE)
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[5]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[5]], 
            smooth=TRUE),col="orange",add = TRUE)
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[6]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[6]], 
            smooth=TRUE),col="magenta",
        legend(.65, .3, legend=unique(data$Tung_single_NA19239_NA19101$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
               cex=1,text.font=2),add = TRUE)
dev.off()
library(png)
png("Rlogo.png")
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[1]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[1]], 
            smooth=TRUE,auc=TRUE),col="red")
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[2]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[2]], 
            smooth=TRUE),col="blue",add = TRUE)
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[3]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[3]], 
            smooth=TRUE),col="dark green",add = TRUE)
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[4]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[4]], 
            smooth=TRUE),col="black",add = TRUE)
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[5]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[5]], 
            smooth=TRUE),col="orange",add = TRUE)
g<-plot(roc(data$Tung_single_NA19239_NA19101$truth_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[6]], 
            data$Tung_single_NA19239_NA19101$prediction_pval[data$Tung_single_NA19239_NA19101$model_pval==unique(data$Tung_single_NA19239_NA19101$model_pval)[6]], 
            smooth=TRUE),col="magenta",
        legend(.65, .3, legend=unique(data$Tung_single_NA19239_NA19101$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
               cex=1,text.font=2),add = TRUE)
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Tung_single_NA19239_NA19101_AUC<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
plot$Tung_single_NA19239_NA19101_AUC
png("Rlogo.png")
plot$Tung_single_NA19239_NA19101_volcano
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Tung_single_NA19239_NA19101_volcano<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
png("Rlogo.png")
plot$Tung_single_NA19239_NA19101_pheatmap
dev.off()
g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)

plot$Tung_single_NA19239_NA19101_pheatmap<-qplot(geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() +
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_blank())
#############################################################################################################################
# ###################################### Simulated Data 10 #########
# L_SyntheticData_10 <- readRDS("~/ROSeq/paper_work/latest_data/pair/L_SyntheticData_10.rds")
# library(gtools)
# library(ggplot2)
# library(plotROC)
# library(NMF)
# library(pheatmap)
# known_pval=c()
# known_padj=c()
# pred_pval=c()
# pred_padj=c()
# method_pval=c()
# method_padj=c()
# j=0
# for(i in plot_seq)
# {
#   j=j+1
#   Truth_Tung <- readRDS("~/ROSeq/paper_work/latest_data/pair/L_SyntheticData_10.rds")$truth_label
#   common=intersect(rownames(L_SyntheticData_10$pred_label),names(Truth_Tung))
#   t_index=match(common,names(Truth_Tung))
#   p_index=match(common,rownames(L_SyntheticData_10$pred_label))
#   known_pval=c(known_pval,as.numeric(Truth_Tung[t_index]))
#   known_padj=c(known_padj,as.numeric(Truth_Tung[t_index]))
#   pred_pval=c(pred_pval,L_SyntheticData_10$pred_label[p_index,i*2-1])
#   pred_padj=c(pred_padj,L_SyntheticData_10$pred_label[p_index,i*2])
#   method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_SyntheticData_10$pred_label)[i*2-1],sep = "_"),dim(L_SyntheticData_10$pred_label[p_index,])[1]))
#   method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_SyntheticData_10$pred_label)[i*2],sep = "_"),dim(L_SyntheticData_10$pred_label[p_index,])[1]))
# }
# data$SyntheticData_10=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)
# 
# ## draw plots
# basicplot_pval <- ggplot(data$SyntheticData_10, aes(d = truth_pval, m = prediction_pval, color = model_pval)) +
#   geom_roc(n.cuts = 0) +
#   style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity")
# ## calculate auc
# calc_auc(basicplot_pval)
# 
# basicplot_padj <- ggplot(data$SyntheticData_10, aes(d = truth_padj, m = prediction_padj, color = model_padj)) +
#   geom_roc(n.cuts = 0) +
#   style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity")
# ## calculate auc
# calc_auc(basicplot_padj)
# 
# known_pval=c()
# known_padj=c()
# pred_pval=c()
# pred_padj=c()
# method_pval=c()
# method_padj=c()
# j=0
# for(i in plot_seq)
# {
#   j=j+1
#   Truth_Tung <- readRDS("~/ROSeq/paper_work/latest_data/pair/L_SyntheticData_10.rds")$truth_label
#   common=intersect(rownames(L_SyntheticData_10$pred_label),names(Truth_Tung))
#   t_index=match(common,names(Truth_Tung))
#   p_index=match(common,rownames(L_SyntheticData_10$pred_label))
#   known_pval=c(known_pval,as.numeric(Truth_Tung[t_index]))
#   known_padj=c(known_padj,as.numeric(Truth_Tung[t_index]))
#   pred_pval=c(pred_pval,L_SyntheticData_10$pred_label[p_index,i*2-1])
#   pred_padj=c(pred_padj,L_SyntheticData_10$pred_label[p_index,i*2])
#   method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_SyntheticData_10$pred_label)[i*2-1],round(calc_auc(basicplot_pval)[j,3],2),sep = "_"),dim(L_SyntheticData_10$pred_label[p_index,])[1]))
#   method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_SyntheticData_10$pred_label)[i*2],round(calc_auc(basicplot_padj)[j,3],2),sep = "_"),dim(L_SyntheticData_10$pred_label[p_index,])[1]))
# }
# data$SyntheticData_10=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)
# 
# ############################# draw plots pval
# basicplot_pval <- ggplot(data$SyntheticData_10, aes(d = truth_pval, m = prediction_pval, color = model_pval)) +
#   geom_roc(n.cuts = 0) + style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") +
#   scale_color_manual(values=c("red", "blue", "dark green","black","orange","magenta"))
# ## calculate auc
# calc_auc(basicplot_pval)
# 
# plot$SyntheticData_10_AUC <- basicplot_pval +
#   ggtitle("ROC") + theme_classic()
# 
# count_data=L_SyntheticData_10$count[rownames(L_SyntheticData_10$pred_label),]
# log2FC=c()
# for(i in 1:dim(count_data)[1])
# {
#   fav=mean(count_data[i,L_SyntheticData_10$group==0])
#   sav=mean(count_data[i,L_SyntheticData_10$group==1])
#   if(fav==0)
#     fav=0.000001
#   if(sav==0)
#     sav=0.000001
#   log2FC=c(log2FC,log2(fav/sav))
# }
# names(log2FC)<-rownames(count_data)
# 
# DE_genes<-rownames(L_SyntheticData_10$pred_label)[L_SyntheticData_10$pred_label$ROSeq_0_tmmvoom_pval<.05 & log2FC>2 ]
# DE_genes=c(DE_genes,rownames(L_SyntheticData_10$pred_label)[L_SyntheticData_10$pred_label$ROSeq_0_tmmvoom_pval<.05 &  log2FC<(-5)])
# DE_data<-L_SyntheticData_10$count[DE_genes,]
# DE_info=data.frame("log2FC"=log2FC,"pval"=L_SyntheticData_10$pred_label$ROSeq_0_tmmvoom_pval)
# plot$SyntheticData_10_volcano <- EnhancedVolcano(DE_info,
#                                                             lab = rownames(DE_info),
#                                                             x = 'log2FC',
#                                                             y = 'pval',
#                                                             xlim = c(-10, 10),
#                                                             title = "Tung_3_2",
#                                                             pCutoff = 0.05,
#                                                             FCcutoff = 2,
#                                                             # col=c('grey75', 'grey50', 'grey25', 'red','blue'),
#                                                             colAlpha = 1)
# 
# 
# col_ann=data.frame("group"=L_SyntheticData_10$group)
# rownames(col_ann)<-colnames(DE_data)
# plot$SyntheticData_10_pheatmap <- pheatmap::pheatmap(log2(DE_data+1),
#                                                                 annotation_col = col_ann,show_rownames = T,show_colnames = F,cluster_cols=F,cluster_rows = F,
#                                                                 main = "heatmap")
# library(grid)
# pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_new_4_pval.pdf")
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[1]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[1]],
#             smooth=TRUE,auc=TRUE),col="red")
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[2]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[2]],
#             smooth=TRUE),col="blue",add = TRUE)
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[3]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[3]],
#             smooth=TRUE),col="dark green",add = TRUE)
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[4]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[4]],
#             smooth=TRUE),col="black",add = TRUE)
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[5]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[5]],
#             smooth=TRUE),col="orange",add = TRUE)
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[6]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[6]],
#             smooth=TRUE),col="magenta",
#         legend(.65, .3, legend=unique(data$SyntheticData_10$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
#                cex=1,text.font=2),add = TRUE)
# dev.off()
# library(png)
# png("Rlogo.png")
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[1]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[1]],
#             smooth=TRUE,auc=TRUE),col="red")
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[2]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[2]],
#             smooth=TRUE),col="blue",add = TRUE)
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[3]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[3]],
#             smooth=TRUE),col="dark green",add = TRUE)
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[4]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[4]],
#             smooth=TRUE),col="black",add = TRUE)
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[5]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[5]],
#             smooth=TRUE),col="orange",add = TRUE)
# g<-plot(roc(data$SyntheticData_10$truth_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[6]],
#             data$SyntheticData_10$prediction_pval[data$SyntheticData_10$model_pval==unique(data$SyntheticData_10$model_pval)[6]],
#             smooth=TRUE),col="magenta",
#         legend(.65, .3, legend=unique(data$SyntheticData_10$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
#                cex=1,text.font=2),add = TRUE)
# dev.off()
# g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)
# 
# plot$SyntheticData_10_AUC<-qplot(geom="blank") +
#   annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
#   geom_point() +
#   theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.line = element_blank())
# plot$SyntheticData_10_AUC
# png("Rlogo.png")
# plot$SyntheticData_10_volcano
# dev.off()
# g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)
# 
# plot$SyntheticData_10_volcano<-qplot(geom="blank") +
#   annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
#   geom_point() +
#   theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.line = element_blank())
# png("Rlogo.png")
# plot$SyntheticData_10_pheatmap
# dev.off()
# g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)
# 
# plot$SyntheticData_10_pheatmap<-qplot(geom="blank") +
#   annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
#   geom_point() +
#   theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.line = element_blank())
# 
# ###################################### Simulated Data 100 #########
# L_SyntheticData_100 <- readRDS("~/ROSeq/paper_work/latest_data/pair/L_SyntheticData_100.rds")
# library(gtools)
# library(ggplot2)
# library(plotROC)
# library(NMF)
# library(pheatmap)
# known_pval=c()
# known_padj=c()
# pred_pval=c()
# pred_padj=c()
# method_pval=c()
# method_padj=c()
# j=0
# for(i in plot_seq)
# {
#   j=j+1
#   Truth_Tung <- readRDS("~/ROSeq/paper_work/latest_data/pair/L_SyntheticData_100.rds")$truth_label
#   common=intersect(rownames(L_SyntheticData_100$pred_label),names(Truth_Tung))
#   t_index=match(common,names(Truth_Tung))
#   p_index=match(common,rownames(L_SyntheticData_100$pred_label))
#   known_pval=c(known_pval,as.numeric(Truth_Tung[t_index]))
#   known_padj=c(known_padj,as.numeric(Truth_Tung[t_index]))
#   pred_pval=c(pred_pval,L_SyntheticData_100$pred_label[p_index,i*2-1])
#   pred_padj=c(pred_padj,L_SyntheticData_100$pred_label[p_index,i*2])
#   method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_SyntheticData_100$pred_label)[i*2-1],sep = "_"),dim(L_SyntheticData_100$pred_label[p_index,])[1]))
#   method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_SyntheticData_100$pred_label)[i*2],sep = "_"),dim(L_SyntheticData_100$pred_label[p_index,])[1]))
# }
# data$SyntheticData_100=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)
# 
# ## draw plots
# basicplot_pval <- ggplot(data$SyntheticData_100, aes(d = truth_pval, m = prediction_pval, color = model_pval)) +
#   geom_roc(n.cuts = 0) +
#   style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity")
# ## calculate auc
# calc_auc(basicplot_pval)
# 
# basicplot_padj <- ggplot(data$SyntheticData_100, aes(d = truth_padj, m = prediction_padj, color = model_padj)) +
#   geom_roc(n.cuts = 0) +
#   style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity")
# ## calculate auc
# calc_auc(basicplot_padj)
# 
# known_pval=c()
# known_padj=c()
# pred_pval=c()
# pred_padj=c()
# method_pval=c()
# method_padj=c()
# j=0
# for(i in plot_seq)
# {
#   j=j+1
#   Truth_Tung <- readRDS("~/ROSeq/paper_work/latest_data/pair/L_SyntheticData_100.rds")$truth_label
#   common=intersect(rownames(L_SyntheticData_100$pred_label),names(Truth_Tung))
#   t_index=match(common,names(Truth_Tung))
#   p_index=match(common,rownames(L_SyntheticData_100$pred_label))
#   known_pval=c(known_pval,as.numeric(Truth_Tung[t_index]))
#   known_padj=c(known_padj,as.numeric(Truth_Tung[t_index]))
#   pred_pval=c(pred_pval,L_SyntheticData_100$pred_label[p_index,i*2-1])
#   pred_padj=c(pred_padj,L_SyntheticData_100$pred_label[p_index,i*2])
#   method_pval=c(method_pval,rep(paste(alpha[j],colnames(L_SyntheticData_100$pred_label)[i*2-1],round(calc_auc(basicplot_pval)[j,3],2),sep = "_"),dim(L_SyntheticData_100$pred_label[p_index,])[1]))
#   method_padj=c(method_padj,rep(paste(alpha[j],colnames(L_SyntheticData_100$pred_label)[i*2],round(calc_auc(basicplot_padj)[j,3],2),sep = "_"),dim(L_SyntheticData_100$pred_label[p_index,])[1]))
# }
# data$SyntheticData_100=data.frame("truth_pval"=known_pval,"truth_padj"=known_padj, "prediction_pval"=pred_pval, "prediction_padj"=pred_padj, "model_pval"=method_pval, "model_padj"=method_padj)
# 
# ############################# draw plots pval
# basicplot_pval <- ggplot(data$SyntheticData_100, aes(d = truth_pval, m = prediction_pval, color = model_pval)) +
#   geom_roc(n.cuts = 0) + style_roc(theme = theme_bw, xlab = "1-Specificity", ylab = "Sensitivity") +
#   scale_color_manual(values=c("red", "blue", "dark green","black","orange","magenta"))
# ## calculate auc
# calc_auc(basicplot_pval)
# 
# plot$SyntheticData_100_AUC <- basicplot_pval +
#   ggtitle("ROC") + theme_classic()
# 
# count_data=L_SyntheticData_100$count[rownames(L_SyntheticData_100$pred_label),]
# log2FC=c()
# for(i in 1:dim(count_data)[1])
# {
#   fav=mean(count_data[i,L_SyntheticData_100$group==1])
#   sav=mean(count_data[i,L_SyntheticData_100$group==2])
#   if(fav==0)
#     fav=0.000001
#   if(sav==0)
#     sav=0.000001
#   log2FC=c(log2FC,log2(fav/sav))
# }
# names(log2FC)<-rownames(count_data)
# 
# DE_genes<-rownames(L_SyntheticData_100$pred_label)[L_SyntheticData_100$pred_label$ROSeq_0_tmmvoom_pval<.05 & log2FC>1 ]
# DE_genes=c(DE_genes,rownames(L_SyntheticData_100$pred_label)[L_SyntheticData_100$pred_label$ROSeq_0_tmmvoom_pval<.05 &  log2FC<(-1)])
# DE_data<-L_SyntheticData_100$count[DE_genes,]
# DE_info=data.frame("log2FC"=log2FC,"pval"=L_SyntheticData_100$pred_label$ROSeq_0_tmmvoom_pval)
# plot$SyntheticData_100_volcano <- EnhancedVolcano(DE_info,
#                                                  lab = rownames(DE_info),
#                                                  x = 'log2FC',
#                                                  y = 'pval',
#                                                  xlim = c(-10, 10),
#                                                  title = "Tung_3_2",
#                                                  pCutoff = 0.05,
#                                                  FCcutoff = 2,
#                                                  # col=c('grey75', 'grey50', 'grey25', 'red','blue'),
#                                                  colAlpha = 1)
# 
# 
# col_ann=data.frame("group"=L_SyntheticData_100$group)
# rownames(col_ann)<-colnames(DE_data)
# plot$SyntheticData_100_pheatmap <- pheatmap::pheatmap(log2(DE_data+1),
#                                                      annotation_col = col_ann,show_rownames = T,show_colnames = F,cluster_cols=F,cluster_rows = F,
#                                                      main = "heatmap")
# library(grid)
# pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_new_4_pval.pdf")
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[1]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[1]],
#             smooth=TRUE,auc=TRUE),col="red")
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[2]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[2]],
#             smooth=TRUE),col="blue",add = TRUE)
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[3]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[3]],
#             smooth=TRUE),col="dark green",add = TRUE)
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[4]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[4]],
#             smooth=TRUE),col="black",add = TRUE)
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[5]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[5]],
#             smooth=TRUE),col="orange",add = TRUE)
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[6]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[6]],
#             smooth=TRUE),col="magenta",
#         legend(.65, .3, legend=unique(data$SyntheticData_100$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
#                cex=1,text.font=2),add = TRUE)
# dev.off()
# library(png)
# png("Rlogo.png")
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[1]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[1]],
#             smooth=TRUE,auc=TRUE),col="red")
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[2]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[2]],
#             smooth=TRUE),col="blue",add = TRUE)
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[3]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[3]],
#             smooth=TRUE),col="dark green",add = TRUE)
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[4]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[4]],
#             smooth=TRUE),col="black",add = TRUE)
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[5]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[5]],
#             smooth=TRUE),col="orange",add = TRUE)
# g<-plot(roc(data$SyntheticData_100$truth_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[6]],
#             data$SyntheticData_100$prediction_pval[data$SyntheticData_100$model_pval==unique(data$SyntheticData_100$model_pval)[6]],
#             smooth=TRUE),col="magenta",
#         legend(.65, .3, legend=unique(data$SyntheticData_100$model_pval),fill=c("red", "blue","dark green","black","orange","magenta"),
#                cex=1,text.font=2),add = TRUE)
# dev.off()
# g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)
# 
# plot$SyntheticData_100_AUC<-qplot(geom="blank") +
#   annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
#   geom_point() +
#   theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.line = element_blank())
# plot$SyntheticData_100_AUC
# png("Rlogo.png")
# plot$SyntheticData_100_volcano
# dev.off()
# g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)
# 
# plot$SyntheticData_100_volcano<-qplot(geom="blank") +
#   annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
#   geom_point() +
#   theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.line = element_blank())
# png("Rlogo.png")
# plot$SyntheticData_100_pheatmap
# dev.off()
# g <- rasterGrob(readPNG("Rlogo.png"), interpolate=TRUE)
# 
# plot$SyntheticData_100_pheatmap<-qplot(geom="blank") +
#   annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
#   geom_point() +
#   theme_classic()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.line = element_blank())
# 
# #############################################################################################################################
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_1_pval.pdf")
plot$Trap_single_T0_T24_AUC
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_volcano_1_pval.pdf")
plot$Trap_single_T0_T24_volcano
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_2_pval.pdf")
plot$Tung_single_NA19098_NA19101_AUC
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_volcano_2_pval.pdf")
plot$Tung_single_NA19098_NA19101_volcano
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_3_pval.pdf")
plot$Tung_single_NA19098_NA19239_AUC
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_volcano_3_pval.pdf")
plot$Tung_single_NA19098_NA19239_volcano
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_4_pval.pdf")
plot$Tung_single_NA19239_NA19101_AUC
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_volcano_4_pval.pdf")
plot$Tung_single_NA19239_NA19101_volcano
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_5_pval.pdf")
plot$SyntheticData_10_AUC
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_volcano_5_pval.pdf")
plot$SyntheticData_10_volcano
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_ROC_6_pval.pdf")
plot$SyntheticData_100_AUC
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_volcano_6_pval.pdf")
plot$SyntheticData_100_volcano
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_heatmap_1_padj.pdf")
plot$Trap_single_T0_T24_pheatmap
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_heatmap_2_padj.pdf")
plot$Tung_single_NA19098_NA19101_pheatmap
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_heatmap_3_padj.pdf")
plot$Tung_single_NA19098_NA19239_pheatmap
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_heatmap_4_padj.pdf")
plot$Tung_single_NA19239_NA19101_pheatmap
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_heatmap_5_pval.pdf")
plot$SyntheticData_10_pheatmap
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_heatmap_6_pval.pdf")
plot$SyntheticData_100_pheatmap
dev.off()

library(ggpubr)
A<-ggarrange(plot$Trap_single_T0_T24_AUC,
             plot$Tung_single_NA19098_NA19101_AUC,
             plot$Tung_single_NA19098_NA19239_AUC,
             plot$Tung_single_NA19239_NA19101_AUC,
             labels = c("A1","B1","C1","D1"),ncol = 1)
B<-ggarrange(plot$Trap_single_T0_T24_volcano,
             plot$Tung_single_NA19098_NA19101_volcano,
             plot$Tung_single_NA19098_NA19239_volcano,
             plot$Tung_single_NA19239_NA19101_volcano,
             labels = c("A2","B2","C2","D2"),ncol = 1)
C<-ggarrange(plot$Trap_single_T0_T24_pheatmap,
             plot$Tung_single_NA19098_NA19101_pheatmap,
             plot$Tung_single_NA19098_NA19239_pheatmap,
             plot$Tung_single_NA19239_NA19101_pheatmap,
             labels = c("A3","B3","C3","D3"),ncol = 1)

pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_combined_pval.pdf")
ggarrange(A,B,C,ncol = 3)
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_combined_ROC_pval.pdf")
ggarrange(plot$Trap_single_T0_T24_AUC,
          plot$Tung_single_NA19098_NA19101_AUC,
          plot$Tung_single_NA19098_NA19239_AUC,
          plot$Tung_single_NA19239_NA19101_AUC,
          labels = c("A","B","C","D"),ncol =2,nrow=2)
dev.off()
pdf("~/ROSeq/paper_work/latest_plot/pair/Final_fig_combined_volvano_heatmap_pval.pdf")
ggarrange(B,C,ncol = 2)
dev.off()
##########################################################33 calculating kappa's score
stat=matrix(0,nrow = 24,ncol = 6)
rn=c()
for (j in 0:(length(names(data))-1))
{
  uv=unique(data[[j+1]]$model_pval)
  for(i in 1:length(uv))
  {
    print(paste(j,i,sep="_"))
    rn=c(rn,paste(names(data)[j+1],uv[i],sep="_"))
    gt=data[[j+1]]$truth_pval[data[[j+1]]$model_pval==uv[i]]
    pv=data[[j+1]]$prediction_pval[data[[j+1]]$model_pval==uv[i]]
    pv[pv<=0.05]=0
    pv[pv>0.05]=1
    stat[(j*6+i),]<-cm(gt,pv)
  }
}
rownames(stat)<-rn
colnames(stat)<-c("mcc","f1_score","precision","recall","kappa","acc")
stat
write.csv(stat,"~/ROSeq/paper_work/latest_plot/pair/stats.csv")
