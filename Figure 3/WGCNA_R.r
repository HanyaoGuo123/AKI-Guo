library('WGCNA')
setwd("D:\\肾损伤数据\\人的数据\\GSE30718")
#test <- read.csv("GSE139061_Eadon_processed_QN_101419.csv",header = TRUE ,row.names = 1)
#test2 <- t(test)
#test3 <- log2(test2 + 1)
#write.csv(test3 , 'exp_data.csv')



datExpr0 <- read.csv('exp_data_for_WGCNA.csv',header = TRUE,row.names = 1)
datTraits <- read.table('AKI_data.txt',header = TRUE ,row.names = 1)

test_data <- t(datExpr0)
m.vars=apply(test_data,1,var)
data.upper <- test_data[which(m.vars > quantile(m.vars, probs = seq(0,1,0.1))[4]),]
out <- t(data.upper)

datExpr0 <- as.data.frame(out) 


#=====================================WGCNA================================
{
    gsg = goodSamplesGenes(datExpr0 , verbose =3)
    gsg$allOK

    #如果gsg$allOK输出的结构为TRUE则说明没有缺失值，若为FALSE，需要用一下函数去重
    if (!gsg$allOK)
    {
        # Optionally, print the gene and sample names that were removed:
        if (sum(!gsg$goodGenes)>0)
            printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
        if (sum(!gsg$goodSamples)>0)
            printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
        # Remove the offending genes and samples from the data:
        datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
    }
}


{
    sampleTree = hclust(dist(datExpr0), method = "average")
    sizeGrWindow(12,9) 
    par(cex = 0.6);
    par(mar = c(0,4,2,0))
    pdf("WGCNA_fig1.pdf",width=12,height=6)
    plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab =1.5,cex = 1.3,
        cex.axis = 1.5, cex.main = 2)
}


{
    abline(h = 18.5, col = "red") 
    dev.off()
    clust = cutreeStatic(sampleTree, cutHeight = 18.5, minSize = 10)
    table(clust)   
    keepSamples = (clust==1)
    datExpr = datExpr0[keepSamples, ]
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
}

{
    femaleSamples = rownames(datExpr);
    traitRows = intersect(femaleSamples, rownames(datTraits));
    datTraits = as.data.frame(datTraits[traitRows,]);
    rownames(datTraits) = traitRows
    colnames(datTraits) = c('OS')
    collectGarbage()
    datTraits$OS <- as.numeric(datTraits$OS)
}

{
    sampleTree2 = hclust(dist(datExpr), method = "average")
    traitColors = numbers2colors(datTraits, signed = FALSE) 
    plotDendroAndColors(sampleTree2, traitColors,
                        #groupLabels = names(datTraits),
                        groupLabels = names(datTraits),
                        dendroLabels = FALSE,
                        main = "Sample dendrogram and trait heatmap")
    #dev.off()
    save(datExpr, datTraits, file = "exp_and_clinical.RData")
}

{
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
        sizeGrWindow(9, 5)
        par(mfrow = c(1,2));
        cex1 = 0.9;

        pdf("WGCNA_fig3.pdf",width=6,height=6)
        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
            xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
            main = paste("Scale independence"));
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
            labels=powers,cex=cex1,col="red");
            abline(h=0.90,col="red")  
            dev.off()
}

{
    pdf("WGCNA_fig4.pdf",width=6,height=6)
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
        xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
        main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()
    sft$powerEstimate
}

{
    
    net = blockwiseModules(datExpr, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3,
                        deepSplit = 2)
}

table(net$colors)

{
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
pdf("WGCNA_fig5.pdf",width=12,height=6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
}

{
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "FemaleLiver-02-networkConstruction-auto.RData")
}

{
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

pdf("WGCNA_fig6.pdf",width=5,height=12)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
}

{
weight = as.data.frame(datTraits$OS);
names(weight) = "weight";
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
}

{
module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));

pdf("WGCNA_turquoise.pdf",width=6,height=6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
}


{

names(datExpr)
names(datExpr)[moduleColors=="green"]

probes = names(datExpr) 
probes2annot = names(datExpr)
sum(is.na(probes2annot))

geneInfo0 = data.frame(substanceBXH = probes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue);

 modOrder = order(-abs(cor(MEs, weight, use = "p")));

 for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")
}

#TOM
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
plotTOM = dissTOM^7;
diag(plotTOM) = NA;
sizeGrWindow(9,9)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
pdf("WGCNA_Fig7.pdf",width=6,height=6)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes", col=myheatcol)
dev.off()

nSelect = 400
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes",col = myheatcol)

library(gplots) 
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')

pdf("WGCNA_Fig7.pdf",width=6,height=6)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=myheatcol)
dev.off()