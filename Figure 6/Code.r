#Grouping samples based on the expression levels of three genes
library(reshape2) 
library(ggplot2)
library(ggpubr)
library(ggsci)

setwd('D:\\AKI文章\\免疫浸润\\根据基因表达量对样本进行分组')
data <- read.table('表达量与免疫细胞丰度.txt' , sep = '\t',header = TRUE , row.names = 1)
data <- data[-c(1:9),]  #Remove the control group

i <- 'COL1A1'
aver <- mean(data[,i])
high_group <- data[data[,i] > aver,]
low_group <- data[data[,i] < aver, ]

high_data <- high_group[,c(4:25)]
high_data$group <- rep('high',dim(high_data)[1])
low_data <- low_group[,c(4:25)]
low_data$group <- rep('low',dim(low_data)[1])

high_low <- rbind(high_data , low_data)

result<-melt(high_low, id.vars = c("group"), 
             measure.vars = colnames(high_low)[1:22],
             variable.name = c('cell'),
             value.name = 'value')

#ggplot2
p <- ggplot(result,aes(x = cell,y = value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = group),label = "p.signif",size=3,method = "t.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")))  #wilcox.test

names <- paste('高表达组与低表达组差异',i,'.pdf',sep = '')
pdf(names,width=12,height=8)
p
dev.off()
