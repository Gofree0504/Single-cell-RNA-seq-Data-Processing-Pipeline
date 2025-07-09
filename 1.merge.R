#提取每个数据集的样本进行整合
#GSE171904
#1、读取文件
rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(data.table)
library(dplyr)
#处理文件到3个文件夹中
setwd("~/LiverFib/1merge")
dir='GSE171904_RAW/' 
fs=list.files('GSE171904_RAW/','^GSM')
fs
library(tidyverse)
samples=str_split(fs,'_',simplify = T)[,1]
##处理数据，将原始文件分别整理为barcodes.tsv.gz，features.tsv.gz和matrix.mtx.gz到各自的文件夹
#批量将文件名改为 Read10X()函数能够识别的名字
lapply(unique(samples),function(x){
  # x = unique(samples)[1]
  y=fs[grepl(x,fs)]
  folder=paste0("GSE171904_RAW/", paste(str_split(y[1],'_',simplify = T)[,1:2], collapse = "_"))
  dir.create(folder,recursive = T)
  #为每个样本创建子文件夹
  file.rename(paste0("GSE171904_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  #重命名文件，并移动到相应的子文件夹里
  file.rename(paste0("GSE171904_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE171904_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})
#数据读取与合并
dir='GSE171904_RAW/'
samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  tmp = Read10X(file.path(dir,pro )) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)#返回创建的Seurat对象，将其存储在sceList中。
}) 
View(sceList)
#merge所有对象
do.call(rbind,lapply(sceList, dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples) 
#names(sce.all@assays$RNA@layers)
#合并layers
sce.all[["RNA"]]$counts 
# Alternate accessor function with the same result
LayerData(sce.all, assay = "RNA", layer = "counts")
#看看合并前后的sce变化
sce.all
sce.all <- JoinLayers(sce.all)
sce.all
#查看数据
dim(sce.all)
#as.data.frame(sce.all@assays$RNA$counts[1:10,1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 
length(sce.all$orig.ident)
#添加metadata分组信息
library(stringr)
phe = sce.all@meta.data
table(phe$orig.ident)
View(phe)
#函数 str_split 用于拆分字符串：
phe$group = str_split(phe$orig.ident,'[_]',simplify = T)[,2] 
#添加变量信息
phe$treated = phe$orig.ident
phe$treated = gsub("GSM\\d+_oil", "oil", phe$treated)  #\\d+指代的是数字
phe$treated = gsub("GSM\\d+_ccl4", "ccl4", phe$treated) 
phe$treated = gsub("GSM\\d+_bdl", "bdl", phe$treated) 
table(phe$treated)
##随后addmetadata即可
new_type = phe
treated<-new_type$treated
sce.all <- AddMetaData(object = sce.all, metadata = treated, col.name = 'treated')
#保存数据
saveRDS(sce.all, file = "GSE171904.rds")

#2、Seurat V5标准流程
rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
sce.all <- readRDS("GSE171904.rds")
#计算线粒体基因比例
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all),ignore.case = T)] 
print(mito_genes) #可能是13个线粒体基因，小鼠数据基因名为小写"^mt-"
#sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
sce.all=PercentageFeatureSet(sce.all, features = mito_genes, col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)
#计算核糖体基因比例
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
print(ribo_genes)
sce.all=PercentageFeatureSet(sce.all,  features = ribo_genes, col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)
#计算红血细胞基因比例
Hb_genes=rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
print(Hb_genes)
sce.all=PercentageFeatureSet(sce.all,  features = Hb_genes,col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
head(sce.all@meta.data)
#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito",
           "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1 
w=length(unique(sce.all$orig.ident))/3+5;w
ggsave(filename="1.Vlnplot1.pdf",plot=p1,width = w,height = 5)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2 
w=length(unique(sce.all$orig.ident))/2+5;w
ggsave(filename="2.Vlnplot2.pdf",plot=p2,width = w,height = 5)
p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3
ggsave(filename="3.Scatterplot.pdf",plot=p3)
#根据上述指标，过滤低质量细胞/基因
sce.all<-subset(sce.all,
                subset=nFeature_RNA>200&
                  nFeature_RNA<5000&
                  percent_mito<20&
                  percent_ribo > 3&
                  percent_hb < 1)
#数据标准化
sce.all <- NormalizeData(sce.all, 
                         normalization.method = "LogNormalize",
                         scale.factor = 1e4) 
#筛选高变基因
sce.all <- FindVariableFeatures(sce.all)
p4 <- VariableFeaturePlot(sce.all) 
p4
#数据归一化
sce.all <- ScaleData(sce.all)
#PCA线性降维
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
##可视化PCA结果
VizDimLoadings(sce.all, dims = 1:2, reduction = "pca")
DimPlot(sce.all, reduction = "pca") + NoLegend()
DimHeatmap(sce.all, dims = 1:12, cells = 500, balanced = TRUE)
########################################
######################################
#GSE199638
#1、读取文件
rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(data.table)
library(dplyr)
#处理文件到3个文件夹中
setwd("~/LiverFib/1merge")
dir='GSE199638_RAW/' 
fs=list.files('GSE199638_RAW/','^GSM')
fs
library(tidyverse)
samples=str_split(fs,'_',simplify = T)[,1]
##处理数据，将原始文件分别整理为barcodes.tsv.gz，features.tsv.gz和matrix.mtx.gz到各自的文件夹
#批量将文件名改为 Read10X()函数能够识别的名字
lapply(unique(samples),function(x){
  # x = unique(samples)[1]
  y=fs[grepl(x,fs)]
  folder=paste0("GSE199638_RAW/", paste(str_split(y[1],'_',simplify = T)[,1:2], collapse = "_"))
  dir.create(folder,recursive = T)
  #为每个样本创建子文件夹
  file.rename(paste0("GSE199638_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  #重命名文件，并移动到相应的子文件夹里
  file.rename(paste0("GSE199638_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE199638_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})
#数据读取与合并
dir='GSE199638_RAW/'
samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  tmp = Read10X(file.path(dir,pro )) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)#返回创建的Seurat对象，将其存储在sceList中。
}) 
View(sceList)
#merge所有对象
do.call(rbind,lapply(sceList, dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples) 
#names(sce.all@assays$RNA@layers)
#合并layers
#sce.all[["RNA"]]$counts 
# Alternate accessor function with the same result
#LayerData(sce.all, assay = "RNA", layer = "counts")
#看看合并前后的sce变化
sce.all
#sce.all <- JoinLayers(sce.all)
#sce.all
#查看数据
dim(sce.all)
#as.data.frame(sce.all@assays$RNA$counts[1:10,1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 
length(sce.all$orig.ident)
#添加metadata分组信息
library(stringr)
phe = sce.all@meta.data
table(phe$orig.ident)
View(phe)
#函数 str_split 用于拆分字符串：
phe$group = str_split(phe$orig.ident,'[_]',simplify = T)[,2] 
#添加变量信息
phe$treated = phe$orig.ident
phe$treated = gsub("GSM\\d+_157", "CTRL", phe$treated)  #\\d+指代的是数字
phe$treated = gsub("GSM\\d+_158", "CTRL", phe$treated) 
phe$treated = gsub("GSM\\d+_199", "KO_MMT", phe$treated) 
phe$treated = gsub("GSM\\d+_201", "KO_MMT", phe$treated) 
table(phe$treated)
##随后addmetadata即可
new_type = phe
treated<-new_type$treated
sce.all <- AddMetaData(object = sce.all, metadata = treated, col.name = 'treated')
#保存数据
saveRDS(sce.all, file = "GSE199638.rds")

#2、Seurat V5标准流程
rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
sce.all <- readRDS("GSE199638.rds")
#计算线粒体基因比例
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all),ignore.case = T)] 
print(mito_genes) #可能是13个线粒体基因，小鼠数据基因名为小写"^mt-"
#sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
sce.all=PercentageFeatureSet(sce.all, features = mito_genes, col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)
#计算核糖体基因比例
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
print(ribo_genes)
sce.all=PercentageFeatureSet(sce.all,  features = ribo_genes, col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)
#计算红血细胞基因比例
Hb_genes=rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
print(Hb_genes)
sce.all=PercentageFeatureSet(sce.all,  features = Hb_genes,col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
head(sce.all@meta.data)
#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito",
           "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1 
w=length(unique(sce.all$orig.ident))/3+5;w
ggsave(filename="1.Vlnplot1.pdf",plot=p1,width = w,height = 5)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2 
w=length(unique(sce.all$orig.ident))/2+5;w
ggsave(filename="2.Vlnplot2.pdf",plot=p2,width = w,height = 5)
p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3
ggsave(filename="3.Scatterplot.pdf",plot=p3)
#根据上述指标，过滤低质量细胞/基因
#过滤，根据箱线图调整
sce.all<-subset(sce.all,
                subset=nFeature_RNA>200&
                  nFeature_RNA<5000&
                  percent_mito<20&
                  percent_ribo > 3&
                  percent_hb < 1)
#数据标准化
sce.all <- NormalizeData(sce.all, 
                         normalization.method = "LogNormalize",
                         scale.factor = 1e4) 
#筛选高变基因
sce.all <- FindVariableFeatures(sce.all)
p4 <- VariableFeaturePlot(sce.all) 
p4
#数据归一化
sce.all <- ScaleData(sce.all)
#PCA线性降维
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
##可视化PCA结果
VizDimLoadings(sce.all, dims = 1:2, reduction = "pca")
DimPlot(sce.all, reduction = "pca") + NoLegend()
DimHeatmap(sce.all, dims = 1:12, cells = 500, balanced = TRUE)
########################################
######################################
#GSE221481
#1、读取文件
rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(data.table)
library(dplyr)
#处理文件到3个文件夹中
setwd("~/LiverFib/1merge")
dir='GSE221481_RAW/' 
fs=list.files('GSE221481_RAW/','^GSM')
fs
library(tidyverse)
samples=str_split(fs,'_',simplify = T)[,1]
##处理数据，将原始文件分别整理为barcodes.tsv.gz，features.tsv.gz和matrix.mtx.gz到各自的文件夹
#批量将文件名改为 Read10X()函数能够识别的名字
lapply(unique(samples),function(x){
  # x = unique(samples)[1]
  y=fs[grepl(x,fs)]
  folder=paste0("GSE221481_RAW/", paste(str_split(y[1],'_',simplify = T)[,1:2], collapse = "_"))
  dir.create(folder,recursive = T)
  #为每个样本创建子文件夹
  file.rename(paste0("GSE221481_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  #重命名文件，并移动到相应的子文件夹里
  file.rename(paste0("GSE221481_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE221481_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})
#数据读取与合并
dir='GSE221481_RAW/'
samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  tmp = Read10X(file.path(dir,pro )) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)#返回创建的Seurat对象，将其存储在sceList中。
}) 
View(sceList)
#merge所有对象
do.call(rbind,lapply(sceList, dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples) 
#names(sce.all@assays$RNA@layers)
#合并layers
#sce.all[["RNA"]]$counts 
# Alternate accessor function with the same result
#LayerData(sce.all, assay = "RNA", layer = "counts")
#看看合并前后的sce变化
sce.all
#sce.all <- JoinLayers(sce.all)
#sce.all
#查看数据
dim(sce.all)
#as.data.frame(sce.all@assays$RNA$counts[1:10,1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 
length(sce.all$orig.ident)
#添加metadata分组信息
library(stringr)
phe = sce.all@meta.data
table(phe$orig.ident)
View(phe)
#函数 str_split 用于拆分字符串：
phe$group = str_split(phe$orig.ident,'[_]',simplify = T)[,2] 
#添加变量信息
phe$treated = phe$orig.ident
phe$treated = gsub("GSM\\d+_gan1134", "normal_saline", phe$treated)  #\\d+指代的是数字
phe$treated = gsub("GSM\\d+_gan1154", "Thioacetamide", phe$treated) 
phe$treated = gsub("GSM\\d+_gan1160", "Thioacetamide+Celecoxib_7.5mg/kg", phe$treated) 
phe$treated = gsub("GSM\\d+_gan1162", "Thioacetamide+Celecoxib_30mg/kg", phe$treated) 
table(phe$treated)
##随后addmetadata即可
new_type = phe
treated<-new_type$treated
sce.all <- AddMetaData(object = sce.all, metadata = treated, col.name = 'treated')
#保存数据
saveRDS(sce.all, file = "GSE221481.rds")

#2、Seurat V5标准流程
rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
sce.all <- readRDS("GSE221481.rds")
#计算线粒体基因比例
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all),ignore.case = T)] 
print(mito_genes) #可能是13个线粒体基因，小鼠数据基因名为小写"^mt-"
#sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
sce.all=PercentageFeatureSet(sce.all, features = mito_genes, col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)
#计算核糖体基因比例
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
print(ribo_genes)
sce.all=PercentageFeatureSet(sce.all,  features = ribo_genes, col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)
#计算红血细胞基因比例
Hb_genes=rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
print(Hb_genes)
sce.all=PercentageFeatureSet(sce.all,  features = Hb_genes,col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
head(sce.all@meta.data)
#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito",
           "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1 
w=length(unique(sce.all$orig.ident))/3+5;w
ggsave(filename="1.Vlnplot1.pdf",plot=p1,width = w,height = 5)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2 
w=length(unique(sce.all$orig.ident))/2+5;w
ggsave(filename="2.Vlnplot2.pdf",plot=p2,width = w,height = 5)
p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3
ggsave(filename="3.Scatterplot.pdf",plot=p3)
#根据上述指标，过滤低质量细胞/基因
#过滤，根据箱线图调整
sce.all<-subset(sce.all,
                subset=nFeature_RNA>200&
                  nFeature_RNA<5000&
                  percent_mito<20&
                  percent_ribo > 3&
                  percent_hb < 1)
#数据标准化
sce.all <- NormalizeData(sce.all, 
                         normalization.method = "LogNormalize",
                         scale.factor = 1e4) 
#筛选高变基因
sce.all <- FindVariableFeatures(sce.all)
p4 <- VariableFeaturePlot(sce.all) 
p4
#数据归一化
sce.all <- ScaleData(sce.all)
#PCA线性降维
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))
##可视化PCA结果
VizDimLoadings(sce.all, dims = 1:2, reduction = "pca")
DimPlot(sce.all, reduction = "pca") + NoLegend()
DimHeatmap(sce.all, dims = 1:12, cells = 500, balanced = TRUE)