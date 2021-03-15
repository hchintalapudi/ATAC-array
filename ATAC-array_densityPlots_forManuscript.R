##title: "ATAC-array analysis visualization"
## This script processes the final file produced by the ATAC-array pipeline and plots the probe distributions in density plots for all samples
#author: "Himanshu Chintalapudi"
#date: "3/1/2020"
library("magicfor")
library("ggplot2")
library("reshape")
library("tidyverse")
library("Rmpfr")
library("matrixStats")
library("readxl")
library("writexl")
library("survival")
library("survminer")
### Read the main LR file produced by the ATACarray pipeline
ATAC_LRs<- read.table("/Users/hchintalapudi/Desktop/ATAC-array/MSK training set/ATACregionLRs_master_edited.txt", header = T, sep = " ")
ATAC_LRs$type <- factor(ATAC_LRs$type, levels = c("CGH", "ctrl", "blue", "red"))
## Read the supplementary Sample clsiification file, also produced by the pipeline 
SampleClassification<- read.table("/Users/hchintalapudi/Desktop/ATAC-array/MSK training set/SampleClassification_MSK_edited.txt", header = T, sep = ",")

SampleClassification$pvalue<- 10^(-SampleClassification$mlogp)
SampleClassification= SampleClassification%>%
  mutate(Sig = (ifelse(pvalue>0.001 & pvalue< 0.01, "t-test p<0.01", ifelse(pvalue>0.01 & pvalue < 0.05, "t-test p<0.05", ifelse(pvalue<0.001, "t-test p<0.001", "t-test: ns")))))
sig<- SampleClassification$Sig

colNames <- names(ATAC_LRs)[6:ncol(ATAC_LRs)]
#map colnames to mlogp values:
names(sig) = colNames

## Generate and save plots for all samples through iteration:
for(i in colNames){
  plt <- ggplot(ATAC_LRs, aes_string(x = i,"type"), mapping = aes_string(x = i, color = "type")) +
    labs(x = "Normalized intensities", y = "Density", title = paste0(gsub(".L2R","",i))) +
    scale_color_manual(labels = c( "CGH (Negative control)", "Ctrl (Positive control)","Blue","Red"), values = c("black","green", "blue", "red")) +
    theme_bw() + 
    geom_line(stat = "density",size = 1.8, linejoin = 'bevel', position = "identity", lineend = "round") + theme(axis.text=element_text(size=25),axis.title=element_text(size=25),legend.text=element_text(size=25),plot.title = element_text(size = 25), legend.title = element_blank())
  plt<- plt + geom_vline(aes_(xintercept= paste0("median(",i, ", na.rm = TRUE)")), linetype = "solid", size =1)
  #print(plt)
  ggsave(plt, file = paste0("plot_", i,".pdf"), width = 14, height = 8.5, dpi = 300)
  Sys.sleep(2)
}


## For Individual Sample plots like in Fig.4b (i and ii):
g<-ggplot(ATAC_LRs, aes_string(x = ATAC_LRs$PT60.L2R,"type"), mapping = aes_string(x = ATAC_LRs$PT60.L2R, color = "type")) +
    labs(x = "Normalized intensities", y = "Density", title = "PT60") +
    scale_color_manual(labels = c( "CGH (Negative control)", "Ctrl (Positive control)","Blue","Red"), values = c("black","green", "blue", "red")) +
    theme_bw() + 
    geom_line(stat = "density",size = 1.8, linejoin = 'bevel', position = "identity", lineend = "round") + theme(axis.text=element_text(size=25),axis.title=element_text(size=25),legend.text=element_text(size=25),plot.title = element_text(size = 25), legend.title = element_blank())
  #geom_vline(data = med, aes(xintercept= med), linetype = "dashed", size =1) +
g +geom_segment(x= 0, y=0, xend = 0, yend=0.52, linetype = "dotted")+
   geom_segment(x=0.53, y=0, xend= 0.53, yend=0.51, linetype = "dotted", color = "blue") +
   geom_segment(x=1.8, y=0, xend= 1.8, yend=0.515, linetype = "dotted", color = "red") + 
   geom_segment(x=2, y=0, xend= 2, yend=0.74, linetype = "dotted", color = "green") +
  annotate("text", x = -0.13, y= 0.3, label = "CGH", angle = 90, size = 7)+
  annotate("text", x = 0.41, y= 0.3, label = "BLUE", angle = 90, size = 7, color = "blue")+
  annotate("text", x = 1.68, y= 0.3, label = "RED", angle = 90, size = 7, color = "red")+
  annotate("text", x = 1.88, y= 0.39, label = "CTRL", angle = 90, size = 7, color = "green")
