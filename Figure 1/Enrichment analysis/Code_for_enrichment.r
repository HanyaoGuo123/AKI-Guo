library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

setwd("D:\\肾损伤数据\\富集分析")
genes <- read.table('D:\\肾损伤数据\\富集分析\\GSE30718_DEGs.txt', header = TRUE )
genes$Gene <- toupper(genes$Gene)

##============================GO=======================

enrich.go <- enrichGO(gene = genes$Gene,  #List of genes to be enriched
    OrgDb = 'org.Hs.eg.db',  #The gene database of a specified species, with an example species being sheep (sheep)
    keyType = 'SYMBOL',  #Specify the given gene name type, for example, using entrze id as an example
    ont = 'ALL',  #GO Ontology: BP、MF、CC、or ALL
    pAdjustMethod = 'fdr',  #Specify p-value correction method
    pvalueCutoff = 1,  #Specify p-value threshold (1 can be specified to output all)
    qvalueCutoff = 1,  #Specify q value threshold (1 can be specified to output all)
    readable = FALSE)


write.table(enrich.go, 'GSE30718_GO.txt', sep = '\t', row.names = FALSE, quote = FALSE)  #output txt

#barplot(enrich.go)  
#dotplot(enrich.go,showCategory=15) 
#cnetplot(enrich.go,circular=FALSE,colorEdge = TRUE,node_label = NULL,showCategory = 10) 
#emapplot(enrich.go) 
#heatplot(enrich.go) 

#output pdf
pdf("GSE30718_GO.pdf",width=7,height=10)
cnetplot(enrich.go,circular=FALSE,colorEdge = TRUE,node_label = NULL,showCategory = 10)
dev.off()


##============================KEGG=======================
gene.df <- bitr(genes$Gene, fromType = "SYMBOL", #FromType refers to which type of data ID your data belongs to
                toType = c("ENTREZID"), #ToType refers to the type of ID you want to convert to, which can be written in multiple or only one form
                OrgDb = org.Hs.eg.db)  #Orgdb refers to which annotation package corresponds to

kegg <- enrichKEGG(
    gene = gene.df[,2], 
    keyType = 'kegg',  #KEGG
    organism = 'hsa',  #For example, oas represents sheep, and other species can change this line
    pAdjustMethod = 'fdr',  #Specify p-value correction method
    pvalueCutoff = 1,  ##Specify p-value threshold (1 can be specified to output all)
    qvalueCutoff = 1)  #Specify q value threshold (1 can be specified to output all)

write.table(kegg, 'GSE30718_KEGG.txt', sep = '\t', quote = FALSE, row.names = FALSE)

#cnetplot(kegg,circular=FALSE,colorEdge = TRUE,node_label = NULL,showCategory = 10) 
#dotplot(kegg,showCategory=15)
#cnetplot(kegg,circular=FALSE,colorEdge = TRUE,node_label = NULL,showCategory = 10)
#barplot(kegg) 

pdf("GSE30718_KEGG_result.pdf",width=7,height=8)
cnetplot(kegg,circular=FALSE,colorEdge = TRUE,node_label = NULL,showCategory = 10)
dev.off()

##============================GESA=======================
#ID
gene.df <- bitr(genes$Gene, fromType = "SYMBOL", 
                toType = c("ENTREZID"), 
                OrgDb = org.Hs.eg.db)  

colnames(genes) <- c('SYMBOL',"log2FC")
for_GESA <- merge(genes , gene.df , by = "SYMBOL",all = F)
for_GESA<- unique(for_GESA)

###################################KEGG
df_all <- for_GESA
df_all_sort <- df_all[order(df_all$log2FC, decreasing = T),] 
gene_fc = df_all_sort$log2FC 
#head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID 
#head(gene_fc)

KEGG <- gseKEGG(gene_fc, organism = "hsa",pvalueCutoff = 1,minGSSize = 5,verbose = FALSE)  #run kegg

#result
sortKEGG<-KEGG[order(KEGG$enrichmentScore, decreasing = T),]
head(sortKEGG)
dim(sortKEGG)
write.csv(sortKEGG,'GSE30718_GSEA_KEGG.csv',row.names = FALSE) 

pdf("GSE30718_GSEA_KEGG.pdf",width=8,height=10)
ridgeplot(KEGG, fill = 'pvalue' , 15 ) +scale_fill_viridis_c(name = "Exp", option = "C")
#scale_fill_continuous('viridis')
dev.off()

######################GO
GO <- gseGO(
  gene_fc,
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
)

sortGO<-GO[order(GO$enrichmentScore, decreasing = T),]
write.csv(sortGO,'GSE30718_GSEA_GO.csv',row.names = FALSE)

pdf("GSE30718_GSEA_GO.pdf",width=8,height=10)
ridgeplot(GO, fill = 'pvalue' , 15 ) +scale_fill_viridis_c(name = "Exp", option = "C")
#scale_fill_continuous('viridis')
dev.off()

#gseaplot2(KEGG, "hsa05200", color = "firebrick", rel_heights=c(1, .2, .6))




