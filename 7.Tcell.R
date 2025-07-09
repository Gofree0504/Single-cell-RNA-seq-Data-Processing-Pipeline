######################################################################
#免疫细胞亚群分析
#T细胞
load("~/LiverFib/1merge/luo/livermerge.Rdata")
setwd("~/LiverFib/1merge/luo/Tcell")
#提取T细胞
Tsubsce = subset(x = sce, celltype == c('T_cell'))
#dimplot
pdf(file="1.T_UMAP.pdf",width=8,height=6)
DimPlot(Tsubsce, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'subgroup')   
dev.off()
#对Tcell进行二次分群

Tsubsce <- SCTransform(Tsubsce, vars.to.regress = "percent.mt", verbose = FALSE)
Tsubsce <- RunPCA(Tsubsce, features = VariableFeatures(object = Tsubsce))
Tsubsce <- RunHarmony(Tsubsce, c("group"))
# #harmony去除批次效应
# RunHarmony(Tsubsce, group.by.vars="orig.ident")
#画肘图
pdf(file="2.ElbowPlot.pdf",width=10,height=6)
ElbowPlot(Tsubsce)
dev.off()
#tSNE
Tsubsce<-RunTSNE(Tsubsce,
                   #使用多少个PC，如果使用的是SCTransform建议多一些
                   dims=1:20)
#UMAP
Tsubsce<-RunUMAP(Tsubsce, reduction = "pca",dims=1:20)
Tsubsce <- FindNeighbors(Tsubsce,dims = 1:20)
# library(clustree)
Tsubsce <- FindClusters(object = Tsubsce, reduction.type="umap", algorithm= 1, resolution = 0.05,verbose = T)#,graph.name = "SCT_snn_res.")
Tsubsce <- FindClusters(object = Tsubsce, reduction.type="umap", algorithm= 1, resolution = 0.1,verbose = T)#,graph.name = "SCT_snn_res.")
Tsubsce <- FindClusters(object = Tsubsce, reduction.type="umap", algorithm= 1, resolution = 0.2,verbose = T)#,graph.name = "SCT_snn_res.")
Tsubsce <- FindClusters(object = Tsubsce, reduction.type="umap", algorithm= 1, resolution = 0.3,verbose = T)#,graph.name = "SCT_snn_res.")
Tsubsce <- FindClusters(object = Tsubsce, reduction.type="umap", algorithm= 1, resolution = 0.3,verbose = T)#, graph.name = "RNA_snn")
Tsubsce <- FindClusters(object = Tsubsce, reduction.type="umap", algorithm= 1, resolution = 0.4,verbose = T)#, graph.name = "RNA_snn")
Tsubsce <- FindClusters(object = Tsubsce, reduction.type="umap", algorithm= 1, resolution = 0.5,verbose = T)#, graph.name = "RNA_snn")
Tsubsce <- FindClusters(object = Tsubsce, reduction.type="umap", algorithm= 1, resolution = 0.6,verbose = T)#, graph.name = "RNA_snn")
Tsubsce <- FindClusters(object = Tsubsce, reduction.type="umap", algorithm= 1, resolution = 0.7,verbose = T)#, graph.name = "RNA_snn")
Tsubsce <- FindClusters(object = Tsubsce, reduction.type="umap", algorithm= 1, resolution = 0.8,verbose = T)#, graph.name = "RNA_snn")
pdf(file="3.clustree.pdf",width=10,height=6)
clustree(Tsubsce@meta.data, prefix = "SCT_snn_res.")
dev.off()
Tsubsce <- FindClusters(object = Tsubsce, 
                          reduction.type="umap", algorithm= 1, resolution = 0.1,
                          verbose = T)

# Look at cluster IDs of the first 5 cells
head(Idents(Tsubsce), 5)
table(Tsubsce$seurat_clusters) 
#Tsubsce <- RunUMAP(Tsubsce, dims = 1:10)
pdf(file="4.SubUMAP.pdf",width=8,height=6)
DimPlot(Tsubsce, reduction = 'umap',group.by = 'seurat_clusters',
        label = TRUE, pt.size = 0.5) #+ NoLegend()
dev.off()
#分别表示MODEL和CONTROL
#提取MODEL中的T细胞
Tsubsce_M = subset(x = Tsubsce, subgroup == c('MODEL'))
#提取CONTROL中的T细胞
Tsubsce_C = subset(x = Tsubsce, subgroup == c('CONTROL'))

#看T细胞在各种造模方式中的差异
pdf(file="5.sub-T-MODEL_UMAP.pdf",width=8,height=6)
DimPlot(Tsubsce, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'group')   
dev.off()

#umap图
pdf(file="5.1.sub-celltype_UMAP.pdf",width=8,height=6)
DimPlot(Tsubsce, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'celltype')   
dev.off()
######################################################
#Doheatmap
#画cluster中高表达基因
# 计算每个聚类的平均表达量
T.markers <- FindAllMarkers(Tsubsce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)
# 找到每个聚类中表达量最高的前10个基因
top10 <- T.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

Tsubsce$group <- factor(x =Tsubsce$group, levels = c('CONTROL','BDL','CCL4','NASH','TAA'))

# pdf(file="27.subDoheatmap1.pdf",width=10,height=10)
# DoHeatmap(subset(Tsubsce,downsample = 100), features = top10$gene,label = F,group.bar.height = 0.1,group.by ="group")
# dev.off()

pdf(file="6.subDoheatmap.pdf",width=10,height=10)
DoHeatmap(Tsubsce, features = top10$gene, group.by = "group")+NoLegend()
dev.off()

pdf(file="7.cluster-subDoheatmap.pdf",width=10,height=10)
DoHeatmap(Tsubsce, features = top10$gene, group.by = "seurat_clusters")+NoLegend()
dev.off()
###########################################################
#GSVA
library(msigdbr)
library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(Seurat)
library(ggplot2)
library(progress)
DefaultAssay(Tsubsce)<-"RNA"
Tsubsce <- NormalizeData(Tsubsce)
#基因集准备（Hallmark基因集）
#genesets = msigdbr(species = "Mus musculus", category = "H")
#genesets = msigdbr(species = "Mus musculus", category = "C2",subcategory = "CP")
genesets = msigdbr(species = "Mus musculus", category = "C5",subcategory = "GO:BP")
genesets <- subset(genesets,select = c("gs_name","gene_symbol"))%>% as.data.frame()
genesets <- split(genesets$gene_symbol,genesets$gs_name)
#提取分组平均表达矩阵
Idents(Tsubsce)<-Tsubsce@meta.data$seurat_clusters
#expr <- AverageExpression(Tsubsce, assays = "RNA",slot = "data")[[1]]
expr <- AverageExpression(Tsubsce, assays = "SCT",slot = "data")[[1]]
#findmarker函数
kupDEG <- FindMarkers(Tsubsce,
                      ident.1 = "MODEL",ident.2 = "CONTROL",#min.pct = 0.5,#设置min.pct = 0.5参数过滤掉那些在50%以下细胞中检测到的基因 
                      logfc.threshold = log(2),# 设置logfc.threshold = log(2)参数过滤掉那些在两个不同组之间平均表达的差异倍数低于2的基因 
                      #min.diff.pct = 0.25,# 设置min.diff.pct = 0.25参数过滤掉那些在两个不同组之间能检测到的细胞比例低于0.25的基因 
                      group.by = 'subgroup',recorrect_umi = FALSE,
                      assay = "SCT",slot="scale.data")
#设置条件筛选logFC和pvalue
kupDEGsig <- kupDEG[with(kupDEG, (abs(avg_diff)>2 & p_val_adj < 0.05 )), ]
write.table(kupDEGsig,file="T-DEG.xls",sep="\t",row.names=T,col.name=NA,quote=F)
rt=read.table("T-DEG.xls",sep="\t",header=T,check.names=F)
colnames(rt)[1] = "gene"   #第一格添加名称
genes=rt$gene 
#a=as.matrix(Tsubsce@assays$RNA@counts)
#a=as.matrix(Kupsubsce[["RNA"]]@data)
expr=expr[genes,]
#GSVA富集分析
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)
#######################################################################
#画代表性基因的小提琴图
Idents(Tsubsce)<-Tsubsce@meta.data$group
pdf(file="9.CD8_Vln.pdf",width=6,height=6)
VlnPlot(Tsubsce,features = "Cd8a",
        split.by = "group",pt.size = 0,ncol = 1)
dev.off()
T_marker <- c("Cd8a","Cd4","Cxcr3","Ccr7","Il7r","Spata2")
pdf(file="24.1.T_featureplot.pdf",width=12,height=6)
FeaturePlot(Tsubsce,features = T_marker,reduction = "umap",cols = c("lightgrey" ,"#DE1F1F"),ncol=3,raster=FALSE)
dev.off()

#########################################################
#选取每个细胞的特异性高表达基因
marker<-c("Ccr7","Cd27",         #Naive_T
          "Cd4","Cd44","Cd38",    #CD4
          "Cd8a","Cd8b1","Spata2",  #CD8
          "Foxp3","Il2ra",   #Treg
          "Cd27","Il7r","Il2rb","Ccr5","Cxcr3",   #memoryT
          "Cd274","Nfat5", "Tox","Ctla4","Cd40",  #ExhaustedT
          "Cd3e","Fcgr3","Fcgr4","Il2rb","Klrk1")  #NK
#点图
DotPlot(Tsubsce, features = unique(marker),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("25.celltype_marker_dot_SCT_scale_data_T.pdf",width = 12,height = 6)
dev.off()
#根据功能命名cluster
new.cluster.ids <- c("0"="CD8_T",
                     "1"="Exhausted_T",
                     "2"="CD4_T",
                     "3"="Treg",#
                     "4"="memory_T",
                     "5"="Naive_T")
Tsubsce@meta.data$celltype<- Tsubsce@meta.data$seurat_clusters
levels(Tsubsce@meta.data$celltype) <- new.cluster.ids#将celltype确定

#点图
DotPlot(Tsubsce, features = unique(marker),group.by = "celltype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("26.celltype_key_marker_dot.pdf",width = 12,height = 6)
dev.off()
###########################################
###############################################
##画比例图
#计算各组细胞百分比
library(ggstatsplot)
#查看各组细胞数
table(Tsubsce$group)
prop.table(table(Tsubsce@meta.data$group))
table(Tsubsce@meta.data$celltype, Tsubsce$group)
#绘制堆叠柱状图
Cellratio <- prop.table(table((Tsubsce@meta.data$celltype), Tsubsce$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
pdf(file="27.Tsubsce-stackplot.pdf",width=10,height=6)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()
#绘制饼图
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(plotrix)
library(ggsci)
library(celldex)
library(singleseqgset)
library(devtools)
head(Tsubsce@meta.data)
table(Tsubsce$celltype)
mynames <-   table(Tsubsce$celltype) %>% names()
myratio <-  table(Tsubsce$celltype) %>% as.numeric()
pielabel <- paste0(mynames," (", round(myratio/sum(myratio)*100,2), "%)")

cols <-c('#E64A35','#4DBBD4' ,'#01A187','#6BD66B','#3C5588'  ,'#F29F80'  ,
         '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
         '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)
pdf(file="28.Tsubsce-pieplot.pdf",width=10,height=6)
pie(myratio, labels=pielabel,
    radius = 1.0,clockwise=T,
    main = "celltype",col = cols)
dev.off()

#绘制3D饼图
pdf(file="29.Tsubsce-3Dpieplot.pdf",width=10,height=6)
pie3D(myratio,labels = pielabel,explode = 0.1, 
      main = "Cell Proption",
      height = 0.3,
      labelcex  = 1)
dev.off()
####################################################################
#GSEA
library(tidyverse)
library(presto)
library(fgsea)
Tsubsce.genes <- wilcoxauc(Tsubsce, 'celltype')
head(Tsubsce.genes)
dplyr::count(Tsubsce.genes, group)# 查看每个cluster中有多少基因（与矩阵的基因数一致）
# 仅选择fgsea的feature和auc列
cluster0.genes<- Tsubsce.genes %>% dplyr::filter(group == "CD4_T") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
ranks<- deframe(cluster0.genes)
head(ranks)
##查看物种的数据 msigdbr_show_species()；#我们使用C7免疫基因集
m_df<- msigdbr(species = "Mus musculus", category = "C7")
head(m_df)
##将m_df的基因与通路取出并改成一个通路对应相应基因的格式
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#以gs_name为factor对gene_symbol进行分类，统计落在每个gs_name中的gene_symbol的个数，并生成list。
summary(fgsea_sets)
#富集分析
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
#nperm设置的是permutation次数
#整理数据
fgseaResTidy <- fgseaRes %>%  as_tibble() %>% arrange(desc(NES))
fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()
#绘图
# 显示top20信号通路
pdf(file="30.Tsubsce-GSEAplot.pdf",width=10,height=6)
ggplot(fgseaResTidy %>% filter(padj < 0.005) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####以7.5进行绘图填色
dev.off()
#GSEA图
pdf(file="31.Tsubsce-GSEA-singleplot.pdf",width=8,height=6)
plotEnrichment(fgsea_sets[["GSE22886_NAIVE_CD8_TCELL_VS_DC_DN"]],
               ranks) + labs(title="GSE22886_NAIVE_CD8_TCELL_VS_DC_DN")
dev.off()

pdf(file="32.Tsubsce-GSEA-singleplot.pdf",width=8,height=6)
plotEnrichment(fgsea_sets[["GSE11057_CD4_EFF_MEM_VS_PBMC_DN"]],
               ranks) + labs(title="GSE11057_CD4_EFF_MEM_VS_PBMC_DN")
dev.off()

pdf(file="33.Tsubsce-GSEA-singleplot.pdf",width=8,height=6)
plotEnrichment(fgsea_sets[["GSE22886_NAIVE_TCELL_VS_MONOCYTE_DN"]],
               ranks) + labs(title="GSE22886_NAIVE_TCELL_VS_MONOCYTE_DN")
dev.off()
####################################################
#在T细胞做拟时序分析
#加载包
suppressMessages({
  library(reshape2)
  library(Seurat)
  library(corrplot)
  library(grid)
  library(cowplot)
  library(tidyverse)
  library(monocle)
  library(ggsci)
  library(ggpubr)
  library(dplyr)
})
#提取T细胞亚群
Tsubsce= subset(x = sce, celltype == c('T_cell'))
#将Seurat对象转为monocle对象
sce_matrix <- as(as.matrix(Tsubsce@assays$RNA@counts), 'sparseMatrix')
pd <-  new('AnnotatedDataFrame', data = Tsubsce@meta.data) 
fd  <-  new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(Tsubsce),row.names = row.names(Tsubsce)))
cds <-  newCellDataSet(sce_matrix,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
#数据均一化，这一步挺慢的
cds  <-   estimateSizeFactors(cds)
cds  <-   estimateDispersions(cds)
#过滤基因
cds <- detectGenes(cds, min_expr = 0.1)#在fDate（cds）中添加一列num_cells_expressed
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))#过滤掉小于10个细胞中表达的基因
#如果Seurat过滤过，这一步可以不用过滤
##使用monocle选择的高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 )
cds <- setOrderingFilter(cds, disp.genes$gene_id)
plot_ordering_genes(cds)
#降维分析，细胞排序
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
#install.packages("monocle_2.26.0.tar.gz", repos = NULL, type = "source")
# 我师兄的解决方法：这个问题是这个包的遗留问题，包的创作者还未解决
# trace('project2MST', edit = T, where = asNamespace("monocle")) 运行此句代码
# 在弹出窗口中把if (elass(projection) != "matrix") 注释，并保存，
# 重新运行你那一句HSMM <- orderCells()即可
cds <- orderCells(cds)
#差异基因
ordering_genes=names(tail(sort(apply(cds@assayData$exprs,1,mad)),2000)) #用作后边的差异基因分析，也可以用前边得到的disp.genes，seurat得到的差异基因，自己指定基因啥的
#原始代码报错
#diff_test_res <- differentialGeneTest(cds[ordering_genes,],fullModelFormulaStr = "~cluster") #这里可以选择用各种基因集，比如说seurat中发现的高变基因，monocle中的基因集等，这里用到了上一步中基因，很多教程中用到了expreesed gene ，但是这个基因集基因有点多，这一步跑的巨慢
#删除部分代码后
diff_test_res <- differentialGeneTest(cds[ordering_genes,],fullModelFormulaStr = "~seurat_clusters") #这里可以选择用各种基因集，比如说seurat中发现的高变基因，monocle中的基因集等，这里用到了上一步中基因，很多教程中用到了expreesed gene ，但是这个基因集基因有点多，这一步跑的巨慢
#diff_test_res <- differentialGeneTest(cds[expressed_genes,]) #这里可以选择用各种基因集，比如说seurat中发现的高变基因，monocle中的基因集等，这里用到了上一步中基因，很多教程中用到了expreesed gene ，但是这个基因集基因有点多，这一步跑的巨慢
deg <- subset(diff_test_res, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),]
##差异基因的结果文件保存
write.table(deg,file="monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
ordergene <- rownames(deg)
Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime)")
#Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #调整列的顺序，可以不用跑
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=2, show_rownames=T, return_heatmap=T)
ggsave("41.Time_heatmapAll.pdf", p, width = 5, height = 10)
#绘制拟时序分析图
p1<- plot_cell_trajectory(cds, color_by = "celltype")
ggsave("42.time-celltype.pdf", p1, width = 8, height = 6)

p2<- plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
ggsave("43.Pseudotime.pdf", p2, width = 8, height = 6)

p3<- plot_cell_trajectory(cds, color_by = "group") + facet_wrap("~State", nrow = 1)#可以看每个cluster的分布情况
ggsave("44.trajectory_group.pdf", plot = p3, width = 16, height = 8)

p4<- plot_cell_trajectory(cds, color_by = "group")
ggsave("45.time-group.pdf", p4, width = 8, height = 6)


