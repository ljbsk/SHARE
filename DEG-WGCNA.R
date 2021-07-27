library("WGCNA")
library(reshape2)
library(stringr)
a=getwd()
setwd(a)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 10)

All_RPKM_anther=read.csv("rpkm-3sets.csv")
All_anther_phe=read.csv("phenotype_samples.csv")
rowname=All_RPKM_anther[,1]
rownames(All_RPKM_anther)=rowname
All_RPKM_anther=All_RPKM_anther[,-1]

traitData=All_anther_phe
traitData=as.data.frame(traitData)
rownamestrait=as.vector(traitData[,1])
rownames(traitData)=rownamestrait
traitData=traitData[,-1]

#m.mad <- apply(All_RPKM_LEAF,1,mad)
#dataExprVar <- All_RPKM_LEAF[which(m.mad > 
#                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(All_RPKM_anther))

gsg = goodSamplesGenes(dataExpr, verbose = 3);
gsg$allOK
#~~~~~~~~~~~~~~~~~~
# 如果上一步返回TRUE则跳过此步，如果返回FALSE则执行如下if语句去掉存在较多缺失值的基因所在行
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
# 再次检测
dim(dataExpr)
gsg = goodSamplesGenes(dataExpr, verbose = 3);
gsg$allOK

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)
pdf("sample_tree.pdf")
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample_tree", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
dev.off()
type="unsigned"

powers = c(c(1:10), seq(from = 12, to=40, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                         verbose=5)

pdf("softbeta.pdf")
par(mfrow = c(1,2));
cex1 = 0.5;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="green")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# 设定软阈值范围
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# 获得各个阈值下的 R方 和平均连接度
sft = pickSoftThreshold(multiExpr[[1]]$data, powerVector = powers, verbose = 5)
# 作图：
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


power = sft$powerEstimate
power

if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}

net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = "pearson", 
                       loadTOMs=TRUE,
                       saveTOMFileBase = "FPKM-TOM",
                       verbose = 3)
table(net$colors)
write.table(table(net$colors),"cluster.txt")

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
pdf("cluster.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)



# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
pdf("Eigengene adjacency.pdf")
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()


sampleName = rownames(dataExpr)
traitData = traitData[match(sampleName, rownames(traitData)), ]

modTraitCor = cor(MEs_col, traitData, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)


textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("module_trait2.pdf",width = 20,height = 20)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 2, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 2, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

save.image("WT_CSA_anther_Frist.Rdata")


library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
# GO 分析：
ego <- enrichGO(gene          = trait_hubGenes_spe,
                # universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

GO_BP <- as.data.frame(ego)
GO_BP$point_shape<-"0"
GO_BP$point_size<-"15"
# write.xlsx(GO_BP,"./results/392_genes_GO_BP.xlsx")

ggplot(data=GO_BP)+
  geom_bar(aes(x=reorder(Description,Count),y=Count, fill=-log10(qvalue)), stat='identity') + 
  coord_flip() +
  scale_fill_gradient(expression(-log["10"]("q value")),low="red", high = "blue") +
  xlab("") +
  ylab("Gene count") +
  scale_y_continuous(expand=c(0, 0))+
  theme_bw()+
  theme(
    axis.text.x=element_text(color="black",size=rel(1.5)),
    axis.text.y=element_text(color="black", size=rel(1.6)),
    axis.title.x = element_text(color="black", size=rel(1.6)),
    legend.text=element_text(color="black",size=rel(1.0)),
    legend.title = element_text(color="black",size=rel(1.1))
    
#heatmap
    
    install.packages("ggpolt2")
    install.packages("pheatmap")
    
    library(pheatmap)
    
    temp = read.table("global1.txt",header = T,row.names = 1,sep="\t")
    head(temp)
    
    bk <- c(seq(-2.5,-0.1,by=0.01),seq(0,3.5,by=0.01))
    
    pheatmap(temp,
             scale = "row",cluster_cols=FALSE,show_rownames=T,cellwidth=24,cellheigh=20,
             color = colorRampPalette(c("blue","white","red"))(600),
             legend_breaks=seq(-2.5,3.5,2),
             breaks=bk)
    
    
    data = read.table("DEGs-candidate.txt",header = T,row.names = 1,sep="\t")
    head(data)
    
    bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
    
    pheatmap(data,
             scale = "row",cluster_cols=FALSE,show_rownames=T,cellwidth=28,cellheigh=21,
             color = colorRampPalette(c("blue","white","red"))(400),
             legend_breaks=seq(-2,2,1),
             breaks=bk)
    
    
    
    data = read.table("DEGs-candidate.txt",header = T,row.names = 1,sep="\t")
    head(data)
    
    
    pheatmap(data,
             scale = "row",cluster_cols=FALSE,show_rownames=T,cellwidth=40,cellheigh=21,
             color = colorRampPalette(c("blue","white","red"))(400),
             legend_breaks = c(-2:2), legend_labels = c("-2","1","0","1","2"))
    
    bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
    pheatmap(data,
             scale = "row",cluster_cols=FALSE,show_rownames=T,cellwidth=40,cellheigh=20,
             color = colorRampPalette(c("blue","white","red"))(400),
             legend_breaks=seq(-2,2,1),
             breaks=bk)
    
# module for WGCNA    
    
    install.packages("WGCNA")
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("impute")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("preprocessCore")
    
    library("WGCNA")
    library(reshape2)
    library(stringr)
    a=getwd()
    setwd(a)
    options(stringsAsFactors = FALSE)
    enableWGCNAThreads()
    
    load("LD_SD_Leaf_Frist")
    
    ### 计算模块与基因的相关性矩阵
    geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
    MMPvalue = as.data.frame(corPvalueStudent(
      as.matrix(geneModuleMembership), nSamples))
    
    # 计算性状与基因的相关性矩阵
    geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
    geneTraitP = as.data.frame(corPvalueStudent(
      as.matrix(geneTraitCor), nSamples))
    
    # 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
    module = "blue"
    pheno = "day"
    modNames = substring(colnames(MEs_col), 3)
    # 获取关注的列
    module_column = match(module, modNames)
    pheno_column = match(pheno,colnames(traitData))
    # 获取模块内的基因
    moduleGenes = moduleColors == module
    
    sizeGrWindow(7, 7)
    par(mfrow = c(1,1))
    # 与性状高度相关的基因，也是与性状相关的模型的关键基因
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                       abs(geneTraitCor[moduleGenes, pheno_column]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", pheno),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    
    
    
    
    
    
    
    
    
    #导出部分网络用于Cytoscape
    load(net$TOMFiles[1], verbose=T)
    TOM <- as.matrix(TOM)
    
    # Select module probes
    probes = colnames(dataExpr) 
    ## 我们例子里面的probe就是基因名
    module = "green"
    inModule = (moduleColors==module);
    modProbes = probes[inModule]; 
    ## 也是提取指定模块的基因名
    # Select the corresponding Topological Overlap
    modTOM = TOM[inModule, inModule];
    dimnames(modTOM) = list(modProbes, modProbes)
    
    cyt = exportNetworkToCytoscape(
      modTOM,
      edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
      nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
      weighted = TRUE,
      threshold = 0.02,
      nodeNames = modProbes, 
      nodeAttr = moduleColors[inModule])
    
    #####################
    probes = colnames(dataExpr)
    dimnames(TOM) <- list(probes, probes)
    
    # Export the network into edge and node list files Cytoscape can read
    # threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
    # cytoscape中再调整
    cyt = exportNetworkToCytoscape(TOM,
                                   edgeFile = paste(exprMat, ".edges.txt", sep=""),
                                   nodeFile = paste(exprMat, ".nodes.txt", sep=""),
                                   weighted = TRUE, threshold = 0,
                                   nodeNames = probes, nodeAttr = moduleColors)
    
    
    ###########
    
    
    #hub geneaaa
    
    # (1) Intramodular connectivity
    #需要linux
    connet=abs(cor(dataExpr,use="p"))^6
    Alldegrees1=intramodularConnectivity(connet, moduleColors)
    head(Alldegrees1)
    # (2) Relationship between gene significance and intramodular connectivity
    which.module="purple"
    Sterile= as.data.frame(traitData[,2]); # change specific 
    names(Sterile) = "Sterile"
    GS1 = as.numeric(cor(Sterile,dataExpr, use="p"))
    GeneSignificance=abs(GS1)
    colorlevels=unique(moduleColors)
    pdf(" gene significance and intramodular connectivity.pdf")
    par(mfrow=c(5,9))
    par(mar=c(5,4,4,2))
    par(cex=0.2,cex.axis=0.2,cex.lab=0.2,cex.main=0.2)
    for (i in c(1:length(colorlevels)))
    {
      whichmodule=colorlevels[[i]];
      restrict1 = (moduleColors==whichmodule);
      verboseScatterplot(Alldegrees1$kWithin[restrict1],
                         GeneSignificance[restrict1], col=moduleColors[restrict1],
                         main=whichmodule,
                         xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
    }
    dev.off()
    
    
    #(3) Generalizing intramodular connectivity for all genes on the array
    datKME=signedKME(dataExpr, MEs_col, outputColumnName="MM.")
    # Display the first few rows of the data frame
    head(datKME)
    ##Finding genes with high gene significance and high intramodular connectivity in
    # interesting modules
    # abs(GS1)> .9 可以根据实际情况调整参数
    # abs(datKME$MM.brown)>.8 至少大于 >0.8
    FilterGenes= abs(GS1)> .8 & abs(datKME$MM.purple)>.8
    table(FilterGenes)
    #提取hub_gene in brown
    hubgene_midnightblue9=rownames(datKME)[FilterGenes]
    write.csv(hubgene_midnightblue9,"hubgene_purple.csv",row.names = F)
    
    #画图@####################################33333
    which.module="turquoise";
    plotMat(t(scale(dataExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=T,
            clabels=T,rcols=which.module,
            title=which.module )
    
    
    which.module="turquoise"
    ME=MEs_col[, paste("ME",which.module, sep="")]
    par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
    plotMat(t(scale(dataExpr[,moduleColors==which.module ]) ),
            nrgcols=30,rlabels=F,rcols=which.module,clabels=rownamestrait
    )
    
    sizeGrWindow(8,4);
    which.module="pink"
    ME=MEs_col[, paste("ME",which.module, sep="")]
    par(mar=c(5, 5, 0.8, 0.7))
    barplot(ME, col=which.module, main="", cex.main=2,cex.lab = 2,
            cex.axis= 1.5, ylab="eigengene expression",xlab="")
    
    
    
    # Select module 基因list
    module = "pink";
    # Select module probes
    probes = colnames(dataExpr) ## 我们例子里面的probe就是基因名
    inModule = (moduleColors==module);
    modProbes = probes[inModule]
    write.csv(modProbes,"pink.csv")
    
    #画一个树与热点
    traitColors = numbers2colors(traitData, signed = TRUE,centered=TRUE);
    plotDendroAndColors(sampleTree, 
                        traitColors, 
                        groupLabels = names(traitData), 
                        rowTextAlignment = "right-justified",
                        addTextGuide = TRUE ,
                        hang = 0.03,
                        dendroLabels = NULL, # 是否显示树labels
                        addGuide = FALSE,  # 显示虚线
                        guideHang = 0.05,
                        main = "Sample dendrogram and trait heatmap")
    
    #(3) Generalizing intramodular connectivity for all genes on the array
    datKME=signedKME(dataExpr, MEs_col, outputColumnName="MM.")
    # Display the first few rows of the data frame
    head(datKME)
    ##Finding genes with high gene significance and high intramodular connectivity in
    # interesting modules
    # abs(GS1)> .9 可以根据实际情况调整参数
    # abs(datKME$MM.brown)>.8 至少大于 >0.8
    
    
    FilterGenes= abs(GS1)> .9 & abs(datKME$MM.)>.8
    table(FilterGenes)
    #提取hub_gene in brown
    hubgene_blue=rownames(datKME)[FilterGenes]
    write.csv(hubgene_blue,"hubgene_yellow2.csv",row.names = F)
    
