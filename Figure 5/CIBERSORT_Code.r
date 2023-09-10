library('e1071')
library('parallel')
library('preprocessCore')
setwd("D:\\肾损伤数据\\CIBERSORT免疫浸润")
source('CIBERSORT.r')
result1 <- CIBERSORT('LM22.txt','GSE139061.txt', perm = 1000, QN = T) 

#===================Visualization of results====================
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(pheatmap)

res <- data.frame(result2[,1:21])%>%
  mutate(group = c(rep('control',9),rep('AKI',39)))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:22],
               names_to = "cell.type",
               values_to = 'value')
head(res,n=6)


p <- ggplot(res,aes(cell.type,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black",trim = FALSE) + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = group),label = "p.signif",size=3,method = "t.test")  #wilcox.test


pdf("Cibersort_result.pdf",width=12,height=8)
p
dev.off()