library("ggplot2")
library("ggExtra")
library("ggpmisc")

setwd("D:\\AKI文章\\与其他肾病的相关性\\补充材料\\GFR_Berthier Lupus TubInt(不全)")

CBFB_data <- read.table('CBFB.txt',header = TRUE)
COL1A1_data <- read.table('COL1A1.txt',header = TRUE)
EGF_data <- read.table('EGF.txt',header = TRUE)


#====================绘制基本散点图=================
color = '#6475e3fe'   #CBFB
#color = '#e36464fe'   #COL1A1
#color = '#64e386fe' 

p <- ggplot(CBFB_data,aes(x = GFR ,y = gene ,color = color))+
            geom_point(size = 2 )+
            scale_color_manual(values=c(color))+
            #scale_color_manual(values=c("#00c000"))+
            theme_bw()+
            theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position = c(0.9, 0),legend.justification = c(0.9,0))+
            #scale_x_log10(limits=c(0.001, 0.1))+scale_y_log10(limits=c(0.001, 0.1))+
            #stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')), formula = y ~ x, parse = T,family = "SH")+
            geom_smooth(method = "lm")


#+xlim(0,0.04)+ylim(0,0.04)

fig<-ggMarginal(p,type="histogram", fill = color , color = 'black' )
#type = "density" 是密度图
#type = "histogram" 是直方图
pdf("CBFB.pdf",width=5,height=5)
fig
dev.off()



#fig<-ggMarginal(p,type="density", xparams = list(bins=40), yparams = list(bins=40),groupColour = TRUE, groupFill = TRUE)

