#===============GSE30718==================
dataset <- read.table("D:\\AKI文章\\原始数据\\Figure 1\\data_for_volcanoplot.txt",header = TRUE,sep = '\t')

cut_off_pvalue = 0.05 
cut_off_logFC = 0.58

dataset$change = ifelse(dataset$P.Value < cut_off_pvalue & abs(dataset$logFC) >= cut_off_logFC, 
                          ifelse(dataset$logFC> cut_off_logFC ,'Up','Down'),
                          'Stable')

p <- ggplot(
  dataset, aes(x = logFC, y = -log10(P.Value), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#2878B5", "#b9b9b9","#c82423"))+
  geom_vline(xintercept=c(-0.58,0.58),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

pdf("GSE30718_volcano plot.pdf",width=7,height=7)
p
dev.off()

#===============GSE139061==================
dataset <- read.csv("D:\\AKI文章\\原始数据\\Figure S1\\control_treat.DESeq2.csv",header = TRUE)

cut_off_pvalue = 0.05 
cut_off_logFC = 0.58 

dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                          ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                          'Stable')

p <- ggplot(
  dataset, aes(x = log2FoldChange, y = -log10(pvalue), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#2878B5", "#b9b9b9","#c82423"))+
  geom_vline(xintercept=c(-0.58,0.58),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

pdf("GSE139061_volcano plot.pdf",width=7,height=7)
p
dev.off()
