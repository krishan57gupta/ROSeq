data=list()
plot=list()
plot_seq_2=c(1,2,3,7,9,10)
p_type_="pval"
process_type_="count"
alpha=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
################################################################# L_Trap_single Data ####################################
###################################### L_Trap_single$T0_T24 Data #########
library(gtools)
library(ggplot2)
library(plotROC)
library(pheatmap)
library(EnhancedVolcano)
library(biomaRt)
################################################################ null data set ########################################
##########
data_median=c()
data_mad=c()
data_med_mad_method=c()
data_samples=c()
################################################### plots according to different number of samples in NA
L_Trap_single <- readRDS("~/ROSeq/paper_work/latest_data/NA19098/L_Tung_single_null_NA19098_10.rds")
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                       cells_=10,
                       plot_seq_=plot_seq_2,
                       p_type_=p_type_,
                       process_type_=process_type_)
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_10 <-framed_data$plots_
###########################
L_Tung_single <- readRDS("~/ROSeq/paper_work/latest_data/NA19098/L_Tung_single_null_NA19098_25.rds")
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=25,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_)
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_25 <-framed_data$plots_
###########################
L_Trap_single <- readRDS("~/ROSeq/paper_work/latest_data/NA19098/L_Tung_single_null_NA19098_50.rds")
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=50,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_)
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_50 <-framed_data$plots_
###########################
L_Trap_single <- readRDS("~/ROSeq/paper_work/latest_data/NA19098/L_Tung_single_null_NA19098_75.rds")
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=75,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_)
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_75 <-framed_data$plots_
###########################
L_Trap_single <- readRDS("~/ROSeq/paper_work/latest_data/NA19098/L_Tung_single_null_NA19098_100.rds")
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=100,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_)
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_100 <-framed_data$plots_
###########################
L_Trap_single <- readRDS("~/ROSeq/paper_work/latest_data/NA19098/L_Tung_single_null_NA19098_125.rds")
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=125,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_)
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_125 <-framed_data$plots_
###########################
L_Trap_single <- readRDS("~/ROSeq/paper_work/latest_data/NA19098/L_Tung_single_null_NA19098_150.rds")
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=150,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_)
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_150 <-framed_data$plots_
#############################################################################################################################
data_med_mad=data.frame("median"=data_median,"MAD"=data_mad,"samples"=data_samples,"methods"=data_med_mad_method)
plot$L_Trap_single_NA19098_NA19098_median <-ggplot(data_med_mad, aes(x=samples, y=median, group=methods)) + 
  geom_line(aes(color=methods)) +
  geom_point(aes(color=methods))+ ggtitle("Tung_single$NA19098_NA19098_Median") +
  scale_color_manual(values=c("red","green","blue","black","cyan","orange"))+
  theme_classic()+
  theme(axis.ticks.x=element_blank())
plot$L_Trap_single_NA19098_NA19098_MAD <-ggplot(data_med_mad, aes(x=samples, y=MAD, group=methods)) + 
  geom_line(aes(color=methods)) +
  geom_point(aes(color=methods))+ ggtitle("Tung_single$NA19098_NA19098_MAD") +
  scale_color_manual(values=c("red","green","blue","black","cyan","orange"))+
  theme_classic()+
  theme(axis.ticks.x=element_blank())
###########---- Time Plot ---- ####################################################
################################################# box plot with time #########################################################
process_type_="time"
data_median=c()
data_mad=c()
data_med_mad_method=c()
data_samples=c()
######################################################################### 10 samples
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=10,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_,
                          time_file_="~/ROSeq/paper_work/latest_data/NA19098/time/null_NA19098_time_10_")
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_10_time <-framed_data$plots_
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=25,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_,
                          time_file_="~/ROSeq/paper_work/latest_data/NA19098/time/null_NA19098_time_25_")
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_25_time <-framed_data$plots_
######################################################################### 50 samples
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=50,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_,
                          time_file_="~/ROSeq/paper_work/latest_data/NA19098/time/null_NA19098_time_50_")
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_50_time <-framed_data$plots_
######################################################################### 75 samples
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=75,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_,
                          time_file_="~/ROSeq/paper_work/latest_data/NA19098/time/null_NA19098_time_75_")
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_75_time <-framed_data$plots_
######################################################################### 100 samples
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=100,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_,
                          time_file_="~/ROSeq/paper_work/latest_data/NA19098/time/null_NA19098_time_100_")
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_100_time <-framed_data$plots_
######################################################################### 125 samples
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=125,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_,
                          time_file_="~/ROSeq/paper_work/latest_data/NA19098/time/null_NA19098_time_125_")
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_125_time <-framed_data$plots_
######################################################################### 150 samples
source("~/ROSeq/paper_work/latest_script/Figures_functions.R")
framed_data<-data_framing(L_=L_Trap_single$NA19098_NA19098_pred_label,
                          cells_=150,
                          plot_seq_=plot_seq_2,
                          p_type_=p_type_,
                          process_type_=process_type_,
                          time_file_="~/ROSeq/paper_work/latest_data/NA19098/time/null_NA19098_time_150_")
data_median=c(data_median,framed_data$data_median_)
data_mad=c(data_mad,framed_data$data_mad_)
data_med_mad_method=c(data_med_mad_method,framed_data$data_med_mad_method_)
data_samples=c(data_samples,framed_data$data_samples_)
plot$L_Trap_single_NA19098_NA19098_boxplot_150_time <-framed_data$plots_
##################################################################################### median mad time
data_med_mad=data.frame("median"=data_median,"MAD"=data_mad,"samples"=factor(data_samples),"methods"=data_med_mad_method)
plot$L_Trap_single_NA19098_NA19098_time_median <-ggplot(data_med_mad, aes(x=samples, y=median, group=methods)) + 
  geom_line(aes(color=methods)) +
  geom_point(aes(color=methods))+ ggtitle("Tung_single$NA19098_NA19098_time_Median") +
  scale_color_manual(values=c("red","green","blue","black","cyan","orange"))+
  theme_classic()+
  theme(axis.ticks.x=element_blank())
plot$L_Trap_single_NA19098_NA19098_time_MAD <-ggplot(data_med_mad, aes(x=samples, y=MAD, group=methods)) + 
  geom_line(aes(color=methods)) +
  geom_point(aes(color=methods))+ ggtitle("Tung_single$NA19098_NA19098_time_MAD") +
  scale_color_manual(values=c("red","green","blue","black","cyan","orange"))+
  theme_classic()+
  theme(axis.ticks.x=element_blank())
#############################################################################################################################
pdf(paste("~/ROSeq/paper_work/latest_plot/NA19098/Final_fig_nulldataset_NA19098_median_",p_type_,".pdf",sep=""),
    width = 7, onefile=FALSE)
plot$L_Trap_single_NA19098_NA19098_median
dev.off()
pdf(paste("~/ROSeq/paper_work/latest_plot/NA19098/Final_fig_nulldataset_NA19098_MAD_",p_type_,".pdf",sep=""),
    width = 7, onefile=FALSE)
plot$L_Trap_single_NA19098_NA19098_MAD
dev.off()
pdf(paste("~/ROSeq/paper_work/latest_plot/NA19098/Final_fig_nulldataset_NA19098_time_median_",p_type_,".pdf",sep=""),
    width = 7, onefile=FALSE)
plot$L_Trap_single_NA19098_NA19098_time_median
dev.off()
pdf(paste("~/ROSeq/paper_work/latest_plot/NA19098/Final_fig_nulldataset_NA19098_time_MAD_",p_type_,".pdf",sep=""),
    width = 7, onefile=FALSE)
plot$L_Trap_single_NA19098_NA19098_time_MAD
dev.off()

library(ggpubr)
pdf(paste("~/ROSeq/paper_work/latest_plot/NA19098/Final_fig_nulldataset_combine_",p_type_,".pdf",sep=""),
    width = 12, onefile=FALSE)
ggarrange(plot$L_Trap_single_NA19098_NA19098_boxplot_25,
          plot$L_Trap_single_NA19098_NA19098_boxplot_50,
          plot$L_Trap_single_NA19098_NA19098_boxplot_75,
          plot$L_Trap_single_NA19098_NA19098_boxplot_100,
          plot$L_Trap_single_NA19098_NA19098_boxplot_125,
          plot$L_Trap_single_NA19098_NA19098_boxplot_150,
          labels = c("A","B","C","D","E","F"),
          common.legend = TRUE)
dev.off()
pdf(paste("~/ROSeq/paper_work/latest_plot/NA19098/Final_fig_nulldataset_combine_time_",p_type_,".pdf",sep=""),
    width = 12, onefile=FALSE)
ggarrange(plot$L_Trap_single_NA19098_NA19098_boxplot_25_time,
          plot$L_Trap_single_NA19098_NA19098_boxplot_50_time,
          plot$L_Trap_single_NA19098_NA19098_boxplot_75_time,
          plot$L_Trap_single_NA19098_NA19098_boxplot_100_time,
          plot$L_Trap_single_NA19098_NA19098_boxplot_125_time,
          plot$L_Trap_single_NA19098_NA19098_boxplot_150_time,
          labels = c("A","B","C","D","E","F"),
          common.legend = TRUE)
dev.off()


