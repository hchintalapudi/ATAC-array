##title: "ATAC-array analysis visualization"
## This script processes the final file produced by the ATAC-array pipeline and plots the probe distributions in density plots for all samples
#author: "Himanshu Chintalapudi"
#date: "3/1/2020"
library("magicfor")
library("reshape")
library("tidyverse")
library("Rmpfr")
library("matrixStats")
library("readxl")
library("writexl")
library("survival")
library("survminer")

ATAC_LRs<- read.table("/Users/hchintalapudi/Desktop/ATAC-array/MSK training set/ATACregionLRs_master_edited.txt", header = T, sep = " ")
ATAC_LRs$type <- factor(ATAC_LRs$type, levels = c("CGH", "ctrl", "blue", "red"))
SampleClassification<- read.table("/Users/hchintalapudi/Desktop/ATAC-array/MSK training set/SampleClassification_MSK_edited.txt", header = T, sep = ",")

SampleClassification$pvalue<- 10^(-SampleClassification$mlogp)
SampleClassification= SampleClassification%>%
  mutate(Sig = (ifelse(pvalue>0.001 & pvalue< 0.01, "t-test p<0.01", ifelse(pvalue>0.01 & pvalue < 0.05, "t-test p<0.05", ifelse(pvalue<0.001, "t-test p<0.001", "t-test: ns")))))
sig<- SampleClassification$Sig

colNames <- names(ATAC_LRs)[6:ncol(ATAC_LRs)]
#map colnames to mlogp values:
names(sig) = colNames

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