#看kupffer细胞与其他细胞的差异
#添加新的metadata分组信息
library(stringr)
phe = sce@meta.data
table(phe$group)
#View(phe)
#函数 str_split 用于拆分字符串：
phe$subgroup = str_split(phe$subgroup,'[_]',simplify = T)[,2] 
#添加变量信息
phe$subgroup = phe$group
phe$subgroup = gsub("BDL", "MODEL", phe$subgroup)  #\\d+指代的是数字
phe$subgroup = gsub("CCL4", "MODEL", phe$subgroup) 
phe$subgroup = gsub("NASH", "MODEL", phe$subgroup) 
phe$subgroup = gsub("TAA", "MODEL", phe$subgroup) 
phe$subgroup = gsub("CONTROL", "CONTROL", phe$subgroup) 
table(phe$subgroup)
##随后addmetadata即可
new_type = phe
subgroup<-new_type$subgroup
sce <- AddMetaData(object = sce, metadata = subgroup, col.name = 'subgroup')
#保存数据
save(sce,file="livermerge.Rdata")
########################################
#featureplot
load("~/LiverFib/1merge/luo/livermerge.Rdata")
#提取kupffer细胞
Kupsubsce = subset(x = sce, celltype == c('Kupffer_cell'))
#dimplot
pdf(file="15.kup_UMAP.pdf",width=8,height=6)
DimPlot(Kupsubsce, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'subgroup')   
dev.off()
#对Kupffer进行二次分群
#SCT可以代替之前NormalizeData, ScaleData, 和 FindVariableFeatures
Kupsubsce <- SCTransform(Kupsubsce, vars.to.regress = "percent.mt", verbose = FALSE)
Kupsubsce <- RunPCA(Kupsubsce, features = VariableFeatures(object = Kupsubsce))
Kupsubsce <- RunHarmony(Kupsubsce, c("group"))
# #harmony去除批次效应
# RunHarmony(Tsubsce, group.by.vars="orig.ident")
#画肘图
pdf(file="16.ElbowPlot.pdf",width=10,height=6)
ElbowPlot(Kupsubsce)
dev.off()
#tSNE
Kupsubsce<-RunTSNE(Kupsubsce,
                 #使用多少个PC，如果使用的是SCTransform建议多一些
                 dims=1:20)
#UMAP
Kupsubsce<-RunUMAP(Kupsubsce, reduction = "pca",dims=1:20)
Kupsubsce <- FindNeighbors(Kupsubsce,dims = 1:20)
# library(clustree)
Kupsubsce <- FindClusters(object = Kupsubsce, reduction.type="umap", algorithm= 1, resolution = 0.05,verbose = T)#,graph.name = "SCT_snn_res.")
Kupsubsce <- FindClusters(object = Kupsubsce, reduction.type="umap", algorithm= 1, resolution = 0.1,verbose = T)#,graph.name = "SCT_snn_res.")
Kupsubsce <- FindClusters(object = Kupsubsce, reduction.type="umap", algorithm= 1, resolution = 0.2,verbose = T)#,graph.name = "SCT_snn_res.")
Kupsubsce <- FindClusters(object = Kupsubsce, reduction.type="umap", algorithm= 1, resolution = 0.3,verbose = T)#,graph.name = "SCT_snn_res.")
Kupsubsce <- FindClusters(object = Kupsubsce, reduction.type="umap", algorithm= 1, resolution = 0.3,verbose = T)#, graph.name = "RNA_snn")
Kupsubsce <- FindClusters(object = Kupsubsce, reduction.type="umap", algorithm= 1, resolution = 0.4,verbose = T)#, graph.name = "RNA_snn")
Kupsubsce <- FindClusters(object = Kupsubsce, reduction.type="umap", algorithm= 1, resolution = 0.5,verbose = T)#, graph.name = "RNA_snn")
Kupsubsce <- FindClusters(object = Kupsubsce, reduction.type="umap", algorithm= 1, resolution = 0.6,verbose = T)#, graph.name = "RNA_snn")
Kupsubsce <- FindClusters(object = Kupsubsce, reduction.type="umap", algorithm= 1, resolution = 0.7,verbose = T)#, graph.name = "RNA_snn")
pdf(file="17.clustree.pdf",width=10,height=6)
clustree(Kupsubsce@meta.data, prefix = "SCT_snn_res.")
dev.off()
Kupsubsce <- FindClusters(object = Kupsubsce, 
                        reduction.type="umap", algorithm= 1, resolution = 0.05,
                        verbose = T)

# Look at cluster IDs of the first 5 cells
head(Idents(Kupsubsce), 5)
table(Kupsubsce$seurat_clusters) 
#Kupsubsce <- RunUMAP(Kupsubsce, dims = 1:10)
pdf(file="18.SubUMAP.pdf",width=10,height=6)
DimPlot(Kupsubsce, reduction = 'umap',group.by = 'subgroup',
        label = TRUE, pt.size = 0.5) #+ NoLegend()
dev.off()
#分别表示MODEL和CONTROL
#提取MODEL中的kupffer细胞
Kupsubsce_M = subset(x = Kupsubsce, subgroup == c('MODEL'))
#提取CONTROL中的kupffer细胞
Kupsubsce_C = subset(x = Kupsubsce, subgroup == c('CONTROL'))
#dimplot
pdf(file="19.Modelkup_UMAP.pdf",width=8,height=6)
DimPlot(Kupsubsce_M, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'subgroup')   
dev.off()

pdf(file="20.Control_UMAP.pdf",width=8,height=6)
DimPlot(Kupsubsce_C, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'subgroup')   
dev.off()

pdf(file="21.subclustertype_UMAP.pdf",width=8,height=6)
DimPlot(Kupsubsce, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'seurat_clusters')   
dev.off()
#看Kupffer细胞在各种造模方式中的差异
pdf(file="22.subMODEL_UMAP.pdf",width=8,height=6)
DimPlot(Kupsubsce, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'group')   
dev.off()
##########################################
#对kupffer细胞进行二次注释
#选取每个细胞的特异性高表达基因
marker<-c("Ly6C1","Nos2","Tnf",    #炎症型
          #"Ym1","Ym2","Mrc1","Cd163","Il10",  #修复型
          "Msr1","Scara1","Cd36",  #脂质代谢型
          "Cd80","Cd86","H2","Cd40","Cd274",  #免疫性Kupffer细胞
          "Waf1","Tnfaip1", #衰老型"
          "Tim4","Mmp14","Ly75","Cd14")  #经典型
#点图
DotPlot(Kupsubsce, features = unique(marker),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("23.celltype_marker_dot_SCT_scale_data_T.pdf",width = 12,height = 6)
dev.off()

new.cluster.ids <- c("0"="Lipid_metabolic",
                     "1"="Inflammatory",
                     "2"="Immune",
                     "3"="Conventional",
                     "4"="Lipid_metabolic",
                     "5"="Senescent")
Kupsubsce@meta.data$celltype<- Kupsubsce@meta.data$seurat_clusters
levels(Kupsubsce@meta.data$celltype) <- new.cluster.ids#将celltype确定

#点图
DotPlot(Kupsubsce, features = unique(marker),group.by = "celltype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("24.celltype_key_marker_dot.pdf",width = 12,height = 6)
dev.off()

pdf(file="25.clustertype_UMAP.pdf",width=8,height=6)
DimPlot(Kupsubsce, reduction = "umap", label = TRUE, repel = TRUE, 
        group.by = 'seurat_clusters',label.size=5,axis.text=5)   
dev.off()

pdf(file="26.celltype_UMAP.pdf",width=8,height=6)
DimPlot(Kupsubsce, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'celltype')   
dev.off()
##################################################
#Doheatmap
#画cluster中高表达基因
# 计算每个聚类的平均表达量
kup.markers <- FindAllMarkers(Kupsubsce, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)
# 找到每个聚类中表达量最高的前10个基因
top10 <- kup.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

Kupsubsce$group <- factor(x =Kupsubsce$group, levels = c('CONTROL','BDL','CCL4','NASH','TAA'))

# pdf(file="27.subDoheatmap1.pdf",width=10,height=10)
# DoHeatmap(subset(Kupsubsce,downsample = 100), features = top10$gene,label = F,group.bar.height = 0.1,group.by ="group")
# dev.off()

pdf(file="27.1.subDoheatmap.pdf",width=10,height=10)
DoHeatmap(Kupsubsce, features = top10$gene, group.by = "group")#+NoLegend()
dev.off()

###################################################
#不同model中kupffer亚群的比例
#画比例图
#计算各组细胞百分比
library(ggstatsplot)
#查看各组细胞数
table(Kupsubsce$group)
prop.table(table(Kupsubsce@meta.data$group))
table(Kupsubsce@meta.data$celltype, Kupsubsce$group)
#绘制堆叠柱状图
Cellratio <- prop.table(table((Kupsubsce@meta.data$celltype), Kupsubsce$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
pdf(file="28.Kupsubsce-stackplot.pdf",width=10,height=6)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()
#画代表性基因的小提琴图
Idents(Kupsubsce)<-Kupsubsce@meta.data$group
pdf(file="29.CD68_Vln.pdf",width=6,height=6)
VlnPlot(Kupsubsce,features = "Cd68",
        split.by = "group",pt.size = 0,ncol = 1)
dev.off()

pdf(file="30.CD14_Vln.pdf",width=6,height=6)
VlnPlot(Kupsubsce,features = "Cd14",
        split.by = "group",pt.size = 0,ncol = 1)
dev.off()
###################################
#GSEA
library(tidyverse)
library(presto)
library(fgsea)
library(msigdbr)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
Kupsubsce.genes <- wilcoxauc(Kupsubsce, 'group')
head(Kupsubsce.genes)
dplyr::count(Kupsubsce.genes, group)# 查看每个cluster中有多少基因（与矩阵的基因数一致）
# 仅选择fgsea的feature和auc列
cluster0.genes<- Kupsubsce.genes %>% dplyr::filter(group == "CONTROL") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
ranks<- deframe(cluster0.genes)
head(ranks)

m_df<- msigdbr(species = "Mus musculus", category = "C5",subcategory = "GO:BP")
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
pdf(file="30.Kupsubsce-GSEAplot.pdf",width=10,height=6)
ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####以7.5进行绘图填色
dev.off()
#GSEA图
pdf(file="31.Kupsubsce-GSEA-singleplot.pdf",width=8,height=6)
plotEnrichment(fgsea_sets[["GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN"]],
               ranks) + labs(title="GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN")
dev.off()
###########################################
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
DefaultAssay(Kupsubsce)<-"RNA"
Kupsubsce <- NormalizeData(Kupsubsce)
#基因集准备（Hallmark基因集）
genesets = msigdbr(species = "Mus musculus", category = "C5")      #subcategory = "HPO")
genesets <- subset(genesets,select = c("gs_name","gene_symbol"))%>% as.data.frame()
genesets <- split(genesets$gene_symbol,genesets$gs_name)
#提取分组平均表达矩阵
Idents(Kupsubsce)<-Kupsubsce@meta.data$group
#expr <- AverageExpression(Kupsubsce, assays = "RNA",slot = "data")[[1]]
expr <- AverageExpression(Kupsubsce, assays = "SCT",slot = "data")[[1]]
#findmarker函数
kupDEG <- FindMarkers(Kupsubsce,
                         ident.1 = "MODEL",ident.2 = "CONTROL",#min.pct = 0.5,#设置min.pct = 0.5参数过滤掉那些在50%以下细胞中检测到的基因 
                         logfc.threshold = log(2),# 设置logfc.threshold = log(2)参数过滤掉那些在两个不同组之间平均表达的差异倍数低于2的基因 
                         #min.diff.pct = 0.25,# 设置min.diff.pct = 0.25参数过滤掉那些在两个不同组之间能检测到的细胞比例低于0.25的基因 
                         group.by = 'subgroup',recorrect_umi = FALSE,
                         assay = "SCT",slot="scale.data")
#设置条件筛选logFC和pvalue
kupDEGsig <- kupDEG[with(kupDEG, (abs(avg_diff)>2 & p_val_adj < 0.05 )), ]
write.table(kupDEGsig,file="kupDEG.xls",sep="\t",row.names=T,col.name=NA,quote=F)
rt=read.table("kupDEG.xls",sep="\t",header=T,check.names=F)
colnames(rt)[1] = "gene"   #第一格添加名称
genes=rt$gene 
#a=as.matrix(Kupsubsce@assays$RNA@counts)
#a=as.matrix(Kupsubsce[["RNA"]]@data)
expr=expr[genes,]
#GSVA富集分析
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="gsva")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)
##########################################
#在kupffer做拟时序分析
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
#提取MODEL中的kupffer细胞
Kupsubsce_M = subset(x = Kupsubsce, subgroup == c('MODEL'))
#提取CONTROL中的kupffer细胞
Kupsubsce_C = subset(x = Kupsubsce, subgroup == c('CONTROL'))
#将Seurat对象转为monocle对象
sce_matrix <- as(as.matrix(Kupsubsce@assays$RNA@counts), 'sparseMatrix')
pd <-  new('AnnotatedDataFrame', data = Kupsubsce@meta.data) 
fd  <-  new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(Kupsubsce),row.names = row.names(Kupsubsce)))
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
p1<- plot_cell_trajectory(cds, color_by = "group")
ggsave("42.time-group.pdf", p1, width = 10, height = 6)

p2<- plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
ggsave("43.Pseudotime.pdf", p2, width = 10, height = 6)

p3<- plot_cell_trajectory(cds, color_by = "group") + facet_wrap("~State", nrow = 1)#可以看每个cluster的分布情况
ggsave("44.trajectory_group.pdf", plot = p3, width = 16, height = 8)


