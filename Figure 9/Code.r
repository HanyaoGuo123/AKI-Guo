#表达趋势验证
library(reshape2)  #wide data => long data by function 'melt(dataset , id = 'Group')'
library(dplyr)   #rownames => coloumns 
library(ggplot2)  #fig 
library(ggpubr)   # '***' P Values
library(plyr)   #计算表格中的平均值和标准差
setwd('D:\\AKI文章\\表达量验证')

#==============================GSE139的数据整理============================
gene_list <- c('CBFB','EGF','COL1A1')
GSE139 <- read.csv("D:\\肾损伤数据\\人的数据\\GSE139061\\exp_data_true.csv", header = TRUE,row.names = 1)
exp_data0 <- GSE139[gene_list,]
exp_data0 <- as.data.frame(t(exp_data0))
exp_data0$group <- c(rep('Ref' , 9), rep('AKI' , 39))
#宽数据整理成长数据格式
result_GSE139 <- melt(exp_data0, id.vars = c("group"), #需保留的不参与聚合的变量列名
             measure.vars = colnames(exp_data0)[1:3],#需要聚合的变量s1-s10
             variable.name = c('gene'),#聚合变量的新列名
             value.name = 'value')#聚合值的新列名
#调整为因子格式，调整顺序
result_GSE139$group <- factor(result_GSE139$group , levels = c('Ref','AKI'))
result_GSE139$gene <- factor(result_GSE139$gene , levels = c('CBFB','EGF','COL1A1'))


#==============================GSE307数据整理==========================
GSE307 <- read.table("D:\\肾损伤数据\\人的数据\\GSE30718\\data_for_ROC.txt",header = TRUE , sep = '\t')
GSE307 <- GSE307[,c(3:5)]
GSE307$group <- c(rep('AKI',28),rep('Ref',11))
#宽数据整理成长数据格式
result_GSE307 <- melt(GSE307, id.vars = c("group"), #需保留的不参与聚合的变量列名
             measure.vars = colnames(GSE307)[1:3],#需要聚合的变量s1-s10
             variable.name = c('gene'),#聚合变量的新列名
             value.name = 'value')#聚合值的新列名
#调整为因子格式，并指定顺序
result_GSE307$group <- factor(result_GSE307$group , levels = c('Ref','AKI'))
result_GSE307$gene <- factor(result_GSE307$gene , levels = c('CBFB','EGF','COL1A1'))


#=========================小鼠叶酸数据集GSE156686============================
GSE156 <- read.csv("D:\\肾损伤数据\\陈飞师哥数据\\GSE156686\\3个基因的log2后的表达量数据.csv",header = TRUE,row.names = 1)
GSE156 <- as.data.frame(t(GSE156))
GSE156$group <- c(rep('Ref',3),rep('AKI',3))
#宽数据整理成长数据格式
result_GSE156 <- melt(GSE156, id.vars = c("group"), #需保留的不参与聚合的变量列名
             measure.vars = colnames(GSE156)[1:3],#需要聚合的变量s1-s10
             variable.name = c('gene'),#聚合变量的新列名
             value.name = 'value')#聚合值的新列名
result_GSE156$group <- factor(result_GSE156$group , levels = c('Ref','AKI'))
result_GSE156$gene <- factor(result_GSE156$gene , levels = c('CBFB','EGF','COL1A1'))
#分组计算标准差和平均数，用来绘制直方图
result_GSE156 <- ddply(result_GSE156,c("group","gene"),summarise,Mean=mean(value,na.rm=TRUE),sd=sd(value,na.rm = TRUE))


#========================小鼠顺铂数据集GSE153625==============================
GSE153 <- read.table('D:\\肾损伤数据\\GSE153625\\用于验证的数据.txt',header = TRUE , row.names = 1)
GSE153 <- log2(GSE153 + 1)
GSE153 <- as.data.frame(t(GSE153))
GSE153$group <- c(rep('Ref',8),rep('AKI',4))
#宽数据整理成长数据格式
result_GSE153 <- melt(GSE153, id.vars = c("group"), #需保留的不参与聚合的变量列名
             measure.vars = colnames(GSE153)[1:2],#需要聚合的变量s1-s10
             variable.name = c('gene'),#聚合变量的新列名
             value.name = 'value')#聚合值的新列名
result_GSE153$group <- factor(result_GSE153$group , levels = c('Ref','AKI'))
result_GSE153$gene <- factor(result_GSE153$gene , levels = c('Egf','Col1a1'))
#分组计算标准差和平均数，用来绘制直方图
result_GSE153 <- ddply(result_GSE153,c("group","gene"),summarise,Mean=mean(value,na.rm=TRUE),sd=sd(value,na.rm = TRUE))

#=====================小鼠顺铂数据GSE165100=================================
GSE165 <- read.table('D:\\肾损伤数据\\GSE165100顺铂\\用于表达量验证.txt',header = TRUE , row.names = 1)
GSE165 <- as.data.frame(t(GSE165))
GSE165$group <- c(rep('AKI',4),rep('Ref',4))
#宽数据整理成长数据格式
result_GSE165 <- melt(GSE165, id.vars = c("group"), #需保留的不参与聚合的变量列名
             measure.vars = colnames(GSE165)[1:2],#需要聚合的变量s1-s10
             variable.name = c('gene'),#聚合变量的新列名
             value.name = 'value')#聚合值的新列名
result_GSE165$group <- factor(result_GSE165$group , levels = c('Ref','AKI'))
result_GSE165$gene <- factor(result_GSE165$gene , levels = c('Cbfb','Egf'))
#分组计算标准差和平均数，用来绘制直方图
result_GSE165 <- ddply(result_GSE165,c("group","gene"),summarise,Mean=mean(value,na.rm=TRUE),sd=sd(value,na.rm = TRUE))


#===========================箱线图的画图==========================
p <- ggplot(result_GSE307 , aes(x = gene , y = value , color = group)) + 
            geom_boxplot(position=position_dodge(0.9),outlier.colour = NA,width=0.5)+
            geom_point(position = position_jitterdodge(dodge.width=0.9),alpha=0.9,size = 2) +
            scale_color_manual(values=c("#6475e3fe", '#e36464fe'))+
            theme_bw()+
            theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),panel.border = element_rect(fill=NA,color="black", linetype="solid"))+
            stat_compare_means(aes(group=group),label="p.signif", label.y = 15,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),method = 't.test') 
            #ylim(7,13.5)

#其他的是15

pdf("GSE307_3个基因表达量验证.pdf",width=5,height=5)
p
dev.off()


#===========================直方图的绘图==========================
p <- ggplot(result_GSE165 , aes(x = gene , y = Mean , fill = group)) + 
            geom_bar(stat="identity", color="black", width = 0.6,position=position_dodge())+
            scale_fill_manual(values=c("#7b7b82", '#147c94'))+
            geom_errorbar(aes(ymin=Mean-sd,ymax=Mean+sd), width=0.2, position=position_dodge(0.6))+
            theme_bw()+
            theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),panel.border = element_rect(fill=NA,color="black", linetype="solid"))
            #stat_compare_means(aes(group=group),label="p.signif", label.y = 15,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),method = 't.test') 
            #ylim(7,13.5)

#其他的是15

pdf("GSE165_2个基因表达量验证.pdf",width=5,height=5)
p
dev.off()