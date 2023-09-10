#============Heat map of differentially expressed genes===========
#=======================GSE30718======================
#==================Data preprocessing=================
library(reshape2) 
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggunchained)
library(cowplot)
library(pheatmap)

data <- read.csv("D:\\AKI文章\\原始数据\\Figure 1\\GSE30718_Gene_exp_data.csv",header = TRUE)
gene_data <- read.table("D:\\AKI文章\\原始数据\\Figure 1\\GSE30718_DEGs_Gene_list.txt",header = TRUE)
exp_data <- data[data$X %in% gene_data$Gene,]
rownames(exp_data) <- exp_data$X
exp_data <- exp_data[,-1]

result1 <- as.data.frame(t(exp_data))
result1 <- result1[-c(29:36),]
result1$Group <- c(rep('AKI',28),rep('control',11))

dois <- result1
#=================== heat map====================
slh <- data.frame(Group=dois$Group) 
slh$Group <- factor(slh$Group,levels = c('AKI','control'))
dois <- select(dois,-Group)
dois <- na.omit(dois)
dois <- log2(dois+1)
dois_new<- data.frame(t(dois))
rownames(slh) <- colnames(dois_new)
ann_colors = list(
  Group = c(AKI ="#c82423", control="#2878B5")
)

p <- pheatmap(dois_new, 
         show_colnames = F, 
         show_rownames=F,  
         fontsize=10, 
         color = colorRampPalette(c("#2878B5", "#ffffff","#c82423"))(50), 
         annotation_col = slh,
         annotation_legend=T,
         annotation_colors=ann_colors,
         #legend_breaks=c(-4,-2,0,2,4),
         #legend_labels=c("-4","-2","0","2","4"),
         border_color=NA,  
         scale="row",  
         cluster_rows = T, 
         cluster_cols = F)

pdf("GSE30718_heatmap.pdf",width=10,height=12)
p
dev.off()

#===============================GSE139061==============================
data0 <- read.csv("D:\\AKI文章\\原始数据\\Figure S1\\GSE139061_Gene_exp_data.csv",header = TRUE,row.names = 1)
data0 <- log2(data0 + 1)
gene_data <- read.table("D:\\AKI文章\\原始数据\\Figure S1\\GSE139061_DEGs_Gene_list.txt",header = TRUE)
gene_list <- gene_data$Gene

data1 <- data0[gene_list , ]
data2 <- as.data.frame(t(data1))

data2$Group <- c(rep('control',9),rep('AKI',39))

dois <- data2

#===========heatmap============
slh <- data.frame(Group=dois$Group) 
slh$Group <- factor(slh$Group,levels = c('AKI','control'))
dois <- select(dois,-Group)
dois <- na.omit(dois)
dois <- log2(dois+1)
dois_new<- data.frame(t(dois))
rownames(slh) <- colnames(dois_new)
ann_colors = list(
  Group = c(AKI ="#c82423", control="#2878B5")
)

p <- pheatmap(dois_new, 
         show_colnames = F, 
         show_rownames=F,  
         fontsize=10, 
         color = colorRampPalette(c("#2878B5", "#ffffff","#c82423"))(50), 
         annotation_col = slh,
         annotation_legend=T,
         annotation_colors=ann_colors,
         #legend_breaks=c(-4,-2,0,2,4),
         #legend_labels=c("-4","-2","0","2","4"),
         border_color=NA,  
         scale="row",  
         cluster_rows = T, 
         cluster_cols = F)

pdf("GSE319061_heatmap.pdf",width=10,height=12)
p
dev.off()



