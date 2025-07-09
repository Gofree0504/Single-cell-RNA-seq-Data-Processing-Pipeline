#手动注释dotplot
#rm(list=ls()) 
library(tidyverse)
library(patchwork)
library(Seurat)
#library(SeuratData)
library(SummarizedExperiment)
DefaultAssay(sce.all) <- 'RNA'
marker<-c(#"CD4","CTLA4","IL2RA",#CD4T_cell
  "Cd8a","Cd8b",#CD8T_cell
  "Kit",'Slc18a2',#Mast_cell
  "Krt19","Krt18","Epcam","Krt8",#Epithelial_cell
  'Jchain','Derl3',"Spc21","Fkbp19","Aoe372",#B_cell 'Cd79a',
  #"CD2","CD3D","CD3E","CD3G",Tcell
  "Cd14","Cd163","Cd68","Cd45r","Scart1",'Fcgr3a',#Macrophage
  "Pecam1","Eng","Vwf",#Endothelial_cell "Flt1",
  "Col1a2","Dcn","Col3a1","Col6a1",#Fibroblast
  #"NKG7","GZMA","KLRF1","KLRD1", #NK_cell
  "Cd83",'Cd16','Cd14a',#'CCR5','CX3CR1',#monocyte 'Elane',
  'Ms4a1','Pax5','Cd79a','Vpreb3','Fcer2')#DC
#点图
DotPlot(sce.all, features = unique(marker),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("cluster_total_marker_dot.pdf",width = 10,height = 6)
new.cluster.ids <- c("0"="Fibroblast",
                     "1"="CD8T_cell",
                     "2"="Mast_cell",
                     "3"="Macrophage",
                     "4"="DC_cell",
                     "5"="Fibroblast",
                     "6"="Macrophage",
                     "7"="Epithelial_cell",
                     "8"="Epithelial_cell",
                     "9"="Epithelial_cell",
                     "10"="Fibroblast",
                     "11"="Macrophage",
                     "12"="Fibroblast",
                     "13"="Fibroblast",
                     "14"="B_cell",
                     "15"="Macrophage",
                     "16"="Endothelial_cell",
                     "17"="Epithelial_cell",
                     "18"="Fibroblast",
                     "19"="Macrophage")
sce.all@meta.data$celltype<- sce.all@meta.data$seurat_clusters
levels(sce.all@meta.data$celltype) <- new.cluster.ids#将celltype确定
#order <- c("CD4T_cell","CD8T_cell","Macrophage","Dendritic_cell","Epithelial_cell","B_cell",
#        "Fibroblast","Endothelial_cell","Mast_cell","NK_cell")
#arrange(order)
#选取每个细胞的特异性高表达基因
#点图
DotPlot(sce.all, features = unique(marker),group.by = "celltype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("11.celltype_key_marker_dot.pdf",width = 12,height = 6)
dev.off()

pdf(file="12.clustertype_UMAP.pdf",width=8,height=6)
DimPlot(sce.all, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'seurat_clusters')   
dev.off()

pdf(file="13.celltype_UMAP.pdf",width=8,height=6)
DimPlot(sce.all, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'celltype')   
dev.off()
###########################################################
#肝脏细胞注释
marker<-c("Cd8a","Cd8b",#CD8T_cell
          "Kit",'Slc18a2',#Mast_cell
          'Ms4a1','Pax5','Cd79a','Vpreb3','Fcer2',#DC
          "Colec11", "Ecm1", "Lrat", "Gucy1b1", "Gucy1a1",#HSC
          "Cdh5", "Bmp2", "Lyve1", "Pecam1", "Cldn5",#Endothelial_cell
          "Csf1r", "Cd74", "Trem2", "Lyz2", "C1qb", "Lgals3",#Kupffer_cell
          "Spp1", "Sorbs2", "Krt18", "Epcam", "Sox9",#Chol_cell
          "Saa3", "lgfbp5", "lgfbp6", "Msln", "Muc16", "Gas1",#Chol_cell
          "Mgp", "Eln", "Mfap4", "Dpt", "Gsn", "Gas6", "Thy1",#Fibroblast
          "CcI5", "Ly6c2", "Ly6d", "Rac2", "Cd3g", "Ptprc",#Leukocyte
          "Alb", "Ttr", "Mup20", "Angptl3", "Gpx1")#Hepatocyte
#点图
DotPlot(sce.all, features = unique(marker),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("10.cluster_total_marker_dot.pdf",width = 10,height = 6)

new.cluster.ids <- c("0"="HSC",
                     "1"="Fibroblast",
                     "2"="Kupffer_cell",
                     "3"="Kupffer_cell",
                     "4"="Leukocyte",
                     "5"="Endothelial_cell",
                     "6"="Chol_cell",
                     "7"="Hepatocyte",
                     "8"="Chol_cell",
                     "9"="Leukocyte",
                     "10"="Endothelial_cell",
                     "11"="Fibroblast",
                     "12"="Leukocyte",
                     "13"="Kupffer_cell")
sce.all@meta.data$celltype<- sce.all@meta.data$seurat_clusters
levels(sce.all@meta.data$celltype) <- new.cluster.ids#将celltype确定
#order <- c("CD4T_cell","CD8T_cell","Macrophage","Dendritic_cell","Epithelial_cell","B_cell",
#        "Fibroblast","Endothelial_cell","Mast_cell","NK_cell")
#arrange(order)
#选取每个细胞的特异性高表达基因
#点图
DotPlot(sce.all, features = unique(marker),group.by = "celltype")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
ggsave("11.celltype_key_marker_dot.pdf",width = 12,height = 6)
dev.off()

pdf(file="12.clustertype_UMAP.pdf",width=8,height=6)
DimPlot(sce.all, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'seurat_clusters')   
dev.off()

pdf(file="13.celltype_UMAP.pdf",width=8,height=6)
DimPlot(sce.all, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'celltype')   
dev.off()

#保存数据
save(sce.all,file="sce.all.rdata")