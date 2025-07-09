#合并多个seurat对象
rm(list=ls())
options(stringsAsFactors = F)
setwd("~/LiverFib/1merge")
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
library(harmony)
GSE171904 <- readRDS("GSE171904.rds")
GSE199638 <- readRDS("GSE199638.rds")
GSE221481 <- readRDS("GSE221481.rds")
#合并
sce.all <- merge(GSE171904, y = c(GSE199638, GSE221481), add.cell.ids = c("GSE171904", "GSE199638", "GSE221481"), project = "merge")
#重新整合
#提取造模亚组
Idents(sce.all)<-sce.all@meta.data$orig.ident
sce.all = sce.all[, Idents(sce.all) %in% c('GSM6873874_gan1154','GSM5979704_199','GSM5979705_201',
                                           'GSM5237106_ccl4','GSM5237107_bdl')]
#添加新的metadata分组信息
library(stringr)
phe = sce.all@meta.data
table(phe$orig.ident)
#View(phe)
#函数 str_split 用于拆分字符串：
phe$group = str_split(phe$orig.ident,'[_]',simplify = T)[,2] 
#添加变量信息
phe$GSEgroup = phe$orig.ident
phe$GSEgroup = gsub("GSM\\d+_bdl", "GSE171904", phe$GSEgroup)  #\\d+指代的是数字
phe$GSEgroup = gsub("GSM\\d+_ccl4", "GSE171904", phe$GSEgroup) 
phe$GSEgroup = gsub("GSM\\d+_199", "GSE199638", phe$GSEgroup) 
phe$GSEgroup = gsub("GSM\\d+_201", "GSE199638", phe$GSEgroup) 
phe$GSEgroup = gsub("GSM\\d+_gan1154", "GSE221481", phe$GSEgroup) 
table(phe$GSEgroup)
##随后addmetadata即可
new_type = phe
GSEgroup<-new_type$GSEgroup
sce.all <- AddMetaData(object = sce.all, metadata = GSEgroup, col.name = 'GSEgroup')
#保存数据
save(sce.all, file = "sce.all.rdata")
#数据标准化
sce.all <- NormalizeData(sce.all, 
                         normalization.method = "LogNormalize",
                         scale.factor = 1e4) 
#筛选高变基因
sce.all <- FindVariableFeatures(sce.all)
# p4 <- VariableFeaturePlot(sce.all) 
# p4
#数据归一化
sce.all <- ScaleData(sce.all)
#PCA线性降维
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
sce.all <- RunHarmony(sce.all, c("orig.ident","GSEgroup"))
names(sce.all@reductions)
#然后使用UMAP/TSNE可视化Harmony去批次效果
sce.all <- RunUMAP(sce.all,  dims = 1:20, 
                   reduction = "harmony")
#DimPlot(sce.all,reduction = "umap",label=F ) 
sce.all <- RunTSNE(sce.all, dims = 1:20, 
                   reduction = "harmony")
#DimPlot(sce.all,reduction = "tsne",label=F )
names(sce.all@reductions)
#保存数据
save(sce.all, file = "sce.all.rdata")

#细胞分群注释
rm(list=ls())
sce.all = read('sce.all.rdata') 
#细胞分群
sce.all <- FindNeighbors(sce.all, reduction = "harmony",
                         dims = 1:20)

for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all=FindClusters(sce.all, #graph.name = "CCA_snn", 
                       resolution = res, algorithm = 1)
}
#画肘图
pdf(file="5.Elbow.pdf",width=10,height=6)
ElbowPlot(sce.all,ndims=20)#根据肘图决定保留多少个PC进行后续分析
dev.off()
#clustree确定分群
library(clustree)
pdf(file="6.clustree.pdf",width=10,height=6)
clustree(sce.all@meta.data, prefix = "RNA_snn_res.")
dev.off()
#确定resolution=0.1
sce.all <- FindClusters(object = sce.all,  resolution = 0.1)
table(sce.all@active.ident) 
colnames(sce.all@meta.data)