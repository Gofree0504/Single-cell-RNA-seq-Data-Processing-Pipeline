rm(list = ls()) 

library(scRNAstat) 
library(ggplot2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(biomaRt)
library(vroom)
library(DT)
library(farver)
library(clusterProfiler)
library(monocle)
library(proxy)
library(farver)
library(reticulate)
library(Seurat)
library(ggplot2)
library(cowplot)
library(devtools)
library(plotly)
library(dplyr)
library(pheatmap)
library(Matrix)
library(data.table)
library(tidyr)
library(farver)

HSC <- subset(sce,subset = celltype %in% "HSC")
DimPlot(HSC)

pbmc <- HSC

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunHarmony(pbmc, group.by.vars = c("sample","study"))
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.2,algorithm = 4) 

#library("reticulate")
#use_python("c:\\Users\\A\\anacon~1\\python.exe", required = T)

DimPlot(pbmc, reduction = "umap",label = TRUE, pt.size = 0.5)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
deg <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pbmc <- subset(pbmc,subset = RNA_snn_res.0.2 %in% c(1:3,6))


new.cluster.ids <- c("Quiescent HSC","Inflammatory HSC","ECM HSC","Mature aHSC","Proliferative HSC")

names(new.cluster.ids) <- levels(pbmc) 
pbmc <- RenameIdents(pbmc, new.cluster.ids)       

pbmc$HSCtype <- Idents(pbmc)

DimPlot(pbmc,label = T)+NoLegend()


pdf(file = "HSC.pdf",height = 10,width = 10)
DimPlot(pbmc,label = T)+NoLegend()
dev.off()

library(Nebulosa)
pdf(file = "HSC1.pdf",height = 10,width = 30)
plot_density(pbmc,features = c("Top2a","Stmn1"),joint = T,reduction = "umap")
dev.off()

pdf(file = "HSC2.pdf",height = 10,width = 30)
plot_density(pbmc,features = c("Angptl6","Vipr1"),joint = T,reduction = "umap")
dev.off()

pdf(file = "HSC3.pdf",height = 10,width = 30)
plot_density(pbmc,features = c("Cxcl14","Ccl7"),joint = T,reduction = "umap")
dev.off()

pdf(file = "HSC4.pdf",height = 10,width = 30)
plot_density(pbmc,features = c("Col1a1","Mmp2"),joint = T,reduction = "umap")
dev.off()



pbmc$group <- factor(pbmc$group,levels = c("CONTROL","BDL","CCL4"))

pbmc$HSCtype <- factor(pbmc$HSCtype,levels = c("Quiescent HSC","Proliferative HSC","Inflammatory HSC","ECM HSC","Mature aHSC"))

VlnPlot(pbmc,features = gene,group.by = "HSCtype",stack = T,pt.size = 0,flip = T)

pdf(file = "S-HSC5.pdf",height = 10,width = 8)
VlnPlot(pbmc,features = gene,group.by = "HSCtype",stack = T,pt.size = 0,flip = T)
dev.off()

pdf(file = "S-HSC6.pdf",height = 10,width = 10)
VlnPlot(pbmc,features = gene,split.by = "HSCtype",
        cols = c("#33539E","#7FACD6","#BF88DA","#E887D4","#A5678E"),
        stack = T,pt.size = 0,flip = T,group.by = "group")
dev.off()

pbmc <- subset(sce,subset = group %in% c("BDL","CCL4","CONTROL"))
table(pbmc$HSCtype)
prop.table(table(Idents(pbmc)))
table(Idents(pbmc), pbmc$group)
Cellratio <- prop.table(table(Idents(pbmc), pbmc$group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2 <- factor(Cellratio$Var2,levels = c("CONTROL","BDL","CCL4"))

allcolour=c("#0000FF","#9370DB","#20B2AA","#FFA500","#DC143C","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
allcolour=c("#EA7369","#C5910E","#79A62C","#2FAF64","#1BB3B7","#359BD7","#9F7DB5","#D369A3")
library(ggplot2)
p <- ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

pdf(file = "HSC6.pdf",height = 8,width = 5)
p
dev.off()


pbmc1 <- subset(pbmc,subset = HSCtype %in% "Inflammatory HSC")
pbmc1 <- subset(pbmc,subset = HSCtype %in% "ECM HSC")
Idents(pbmc1) <- pbmc1$group
table(pbmc1$group)
mydeg <- FindMarkers(pbmc1,ident.1 = c("CCL4"),ident.2 = c("BDL"), 
                     logfc.threshold = 0,verbose = T,
                     test.use = 'wilcox',min.pct = 0.1)

mydeg.df <- bitr(rownames(mydeg),fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                 OrgDb = org.Mm.eg.db)
mydeg$SYMBOL <-rownames(mydeg)
genetable <- merge(mydeg.df,mydeg,by='SYMBOL')  
geneList <- genetable$avg_log2FC
names(geneList) <- genetable$ENTREZID
geneList <- sort(geneList,decreasing = T)
egseGO <- gseGO(geneList,ont = "BP", org.Mm.eg.db, exponent = 1, nPerm = 1000, 
                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.2, pAdjustMethod = "BH", 
                verbose = TRUE, seed = FALSE, by = "fgsea")
egseGO <- setReadable(egseGO,'org.Mm.eg.db','ENTREZID')
geasres <- egseGO@result

egsekegg <- gseKEGG(geneList=geneList,
                    organism = "mmu",  
                    nPerm=1000,
                    minGSSize = 10,
                    pvalueCutoff = 0.2,
                    verbose = FALSE)
egsekegg <- setReadable(egsekegg,'org.Mm.eg.db','ENTREZID')
kegggseares <-egsekegg@result




library(GseaVis)
# BDL GO:0007249 GO:0006954
# CCL4 GO:0030199 GO:0009612
geneSetID = c("GO:0009612")


p2 = gseaNb(object = egseGO,
            geneSetID = geneSetID,
            subPlot = 2,
            termWidth = 35,
            legend.position = c(0.7,0.8),
            ght.relHight = 0.1)
p2


pdf(file = "BDL.pdf",height = 10,width = 10)
p2
dev.off()





rm(list = ls())
library(monocle3)
library(Seurat)

load("HSC.Rdata")

data <- GetAssayData(pbmc,assay = "RNA",slot = "counts")
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds,num_dim = 20)
plot_pc_variance_explained(cds)


cds <- reduce_dimension(cds,preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))

p1 <-plot_cells(cds,reduction_method = "UMAP",color_cells_by = "seurat_clusters") +ggtitle(("cds.umap"))
p1

cds.embd <- cds@int_colData$reducedDims$UMAP
int.embd <- Embeddings(pbmc,reduction = "umap")
int.embd <-int.embd[rownames(cds.embd),]
cds@int_colData$reducedDims$UMAP <- int.embd

p2 <- plot_cells(cds,reduction_method = "UMAP",color_cells_by = "seurat_clusters") +ggtitle(("cds.umap"))
p2

cds <- cluster_cells(cds)
plot_cells(cds,color_cells_by = "partition")
cds <- learn_graph(cds)

p <- plot_cells(cds,color_cells_by = "seurat_clusters",label_groups_by_cluster = F,
                label_leaves = F,label_branch_points = F)
p
cds <- order_cells(cds)

p <- plot_cells(cds,color_cells_by = "pseudotime",label_groups_by_cluster = F,
                label_leaves = F,label_branch_points = F)
p

ggsave("Trajectory_HSC.pdf",plot = p,width = 8,height = 7)

cell_groups <- colData(cds)$group
group1_cds <- cds[, cell_groups == "BDL"]
group2_cds <- cds[, cell_groups == "CCL4"]

p1 <- plot_cells(group1_cds, color_cells_by = "pseudotime", 
                 label_groups_by_cluster = TRUE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE,
                 show_trajectory_graph = FALSE)
p1
ggsave("Trajectory_HSC-BDL.pdf",plot = p1,width = 8,height = 7)


p2 <- plot_cells(group2_cds, color_cells_by = "pseudotime", 
                 label_groups_by_cluster = TRUE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE,
                 show_trajectory_graph = FALSE)
ggsave("Trajectory_HSC-CCL4.pdf",plot = p2,width = 8,height = 7)










