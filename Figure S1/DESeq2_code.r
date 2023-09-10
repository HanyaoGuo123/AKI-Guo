library("DESeq2")
library('org.Mm.eg.db')
library('clusterProfiler')

setwd("D:\\肾损伤数据\\人的数据\\GSE139061")
coldata <- data.frame(condition = factor(
    c(rep('control', 9),rep('treat',39)), levels = c('control', 'treat')))
#coldata <- data.frame(condition = factor(rep(c('control', 'treat'), each = 3), levels = c('control', 'treat')))


dat <- read.csv('GSE139061_Eadon_processed_QN_101419.csv',header = TRUE,row.names = 1)
dat <- round(dat)


#DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design= ~condition)

#For p
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('condition', 'treat', 'control'))
res
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.csv(res1, 'control_treat.DESeq2.csv')




