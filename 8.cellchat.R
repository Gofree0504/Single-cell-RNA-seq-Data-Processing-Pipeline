##############################################
#细胞通讯
#加载R包
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(tidyverse)
rm(list = ls())
options(stringsAsFactors = FALSE)

sce
sce@commands$FindClusters
#在Seurat对象中提出CellChat需要的数据
data.input  <- sce@assays$RNA@counts
library(future)
options(future.globals.maxSize = 1000 * 1024^2)#将其设置为1GB
data.input<- NormalizeData(data.input)
identity = data.frame(group =sce$celltype, row.names = names(sce$celltype)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels
#创建一个Cell Chat对象
cellchat <- createCellChat(object  = data.input)
cellchat
summary(cellchat)
#把metadata信息加到CellChat对象中
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
#导入配受体数据库
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# set the used database in the object
cellchat@DB <- CellChatDB.use
#预处理用于细胞通信分析的表达数据
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
#save(cellchat, file = "cellchat_input_step1.rdata")
#load(file = '/home/ug1243/Gofree/GSE183904/10cellchat/cellchat_input_step1.rdata')
future::plan("multisession", workers = 20) # do parallel
#future::plan("multiprocess", workers = 20)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
#########第二部分：细胞通信网络的推断#######
#计算通信概率并推断cellchat网络
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#在信号通路级别推断细胞-细胞通信
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞通信网络
cellchat <- aggregateNet(cellchat)
##########################################################
#使用圆图显示任意两个细胞组之间的相互作用次数或总交互强度（比重）
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf(file="1Circle plot.pdf",width=10,height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
#控制参数edge.weight.max，以便我们可以比较不同网络之间的边缘权重。
mat <- cellchat@net$weight
par(mfrow = c(1,1), xpd=TRUE)
pdf(file="2Single_Circle.pdf",width=10,height=6)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  #pdf(file="2Single Circle plot.pdf",width=10,height=6)
  #win.graph(width=6, height=5,pointsize=9) 
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
##########第三部分：细胞通信网络的可视化###
#所有显示重要通信的信号通路均可通过cellchat@netP$pathways获取。
cellchat@netP$pathways
pathways.show <- c("CSF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
pdf(file="3netVisual_CSF.pdf",width=10,height=6)
#win.graph(width=6, height=5,pointsize=9)
netVisual_aggregate(cellchat, signaling = pathways.show,  
                    vertex.receiver = vertex.receiver)
dev.off()
# Circle plot
par(mfrow=c(1,1))
pdf(file="4Circle_TGFb.pdf",width=10,height=6)
#win.graph(width=6, height=5,pointsize=9)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
# Chord diagram
par(mfrow=c(1,1))
pdf(file="5Chord_TGFb.pdf",width=8,height=10)
#win.graph(width=5, height=4,pointsize=8)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
# Heatmap
par(mfrow=c(1,1))
pdf(file="6Heatmap_TGFb.pdf",width=8,height=10)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
#> Do heatmap based on a single object
#计算每个配体受体对整体信号通路的贡献，并可视化由单个配体受体对调节的细胞通信
pdf(file="7netVisual_contribution.pdf",width=10,height=6)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
#我们还可以可视化由单个配体受体对调节的细胞-细胞通信。我们提供一个函数extractEnrichedLR来提取给定信号通路的所有重要相互作用（L-R对）和相关信号基因。
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
pdf(file="8Hierarchy_TGFb_LR.pdf",width=10,height=6)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
dev.off()
# Circle plot
pdf(file="9Circle_TGFb_LR.pdf",width=10,height=6)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
dev.off()
# Chord diagram
par(mfrow=c(2,2))
pdf(file="10Chord_TGFb_LR.pdf",width=8,height=10)
#win.graph(width=10, height=8,pointsize=5)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
############################################
#找各组之间差异的通路
meta = sce@meta.data
#提取分组
cell.PT = rownames(meta)[meta$group == "CONTROL"] # extract the cell names from disease data
cell.M = rownames(meta)[meta$group == "TAA"] # extract the cell names from disease data

gene_exp = sce@assays$RNA@counts
library(future)
options(future.globals.maxSize = 100000 * 1024^2)#将其设置为1GB
gene_exp<- NormalizeData(gene_exp)

data.input.M = gene_exp[, cell.M]
meta.M = meta[cell.M, ]

data.input.PT = gene_exp[, cell.PT]
meta.PT = meta[cell.PT, ]

sce.M <- CreateSeuratObject(counts = data.input.M,
                            meta.data = meta.M,
                            
                            project = "sce")

sce.PT <- CreateSeuratObject(counts = data.input.PT,
                             meta.data = meta.PT,
                             project = "sce")  # min.cells = 50

# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$celltype) # check the cell labels

### 创建cellchat对象
cco.M <- createCellChat(sce.M@assays$RNA@data, meta = sce.M@meta.data, group.by = "celltype")
cco.PT <- createCellChat(sce.PT@assays$RNA@data, meta = sce.PT@meta.data, group.by = "celltype")

cco.M@data[1:20,1:20]
#save(cco.M, cco.PT, file = "cellchat_coo_M_PT.Rdata")

# 2. 细胞通讯网络分析-------------------------
# 2.1 数据准备和路径切换
# dir.create("./Cellchat/")
# setwd("Cellchat/")
# load("../cco.rda")
cco.PT <- setIdent(cco.PT, ident.use = "celltype")
cco.M <- setIdent(cco.M, ident.use = "celltype")
table(cco.PT@idents)
table(cco.M@idents)

# 2.2 分析样本cco.normal的细胞通讯网络
# ⚠️：cellchat不太稳定，identifyOverExpressedGenes常出错（不出现进度条就是出错了），重启Rstudio后再运算。
cellchat <- cco.PT
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#########################################################
# 信号网络的流行学习与分类
# 把共同起作用的信号通路归纳在一起，分为基于功能的归纳和基于拓扑结构的归纳
# 基于功能相似性的流行学习与分类
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
# 基于拓扑相似性的流行学习与分类
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
#######################################################
cco.PT <- cellchat
#save(cco.PT, file = "cco.PT.rdata")
# 2.3 分析样本cco.tumor的细胞通讯网络
cellchat <- cco.M
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#################################################################
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
#################################################################
cco.M <- cellchat
#save(cco.M, file = "cco.M.rdata")
# load(file = 'cco.M.rdata')
# load(file = 'cco.PT.rdata')
# 2.4 合并cellchat对象
cco.list <- list(PT=cco.PT, M=cco.M)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)
save(cco.list, cellchat,file = 'merge_cellchat_list.rdata')
save(cellchat, file = "cellchat_input_step1.rdata")
#############################################################
# 3. 将整合的cellchat可视化---------------------
# 3.1 所有细胞群总体观：通讯数量与强度对比------------------
# 比较交互总数和交互强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
p
ggsave("1Overview_number_strength.pdf", p, width = 8, height = 6)
#   结果解释：左图展示通讯数量之间的差异，右图展示通讯强度之间的差异。本例中信号通路强度weight值过低，导致显示时
# 均为0（实际上有数值的，只是过小，显示为0）
# 数量与强度差异网络图
par(mfrow = c(1,2), xpd=TRUE)
pdf(file="2netVisual_diffInteraction.pdf",width=10,height=6)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight") 
dev.off()
#手动保存，8 * 12
# 两个数据集之间细胞通信网络中交互或交互强度的差异数可以使用圆图可视化， 
# 与第一个数据集相比，[红色]上调，[蓝色]边下调
# 数量与强度差异热图
#   我们还可以使用热图在更大的细节中显示交互的差异数或交互强度。顶部彩色条形图表示热图（传入信号）
# 中显示的列值的总和。右边的彩色条形图表示一行值（传出信号）的总和。在色条中红色或蓝色表示第二个
# 数据集中与第一个数据集相比增加或[减少]信号。
par(mfrow = c(1,1), xpd=TRUE)
pdf(file="3netVisual_heatmap.pdf",width=10,height=6)
h1 = netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
p = h1+h2
p
dev.off()
# save as Diff_number_strength_heatmap.pdf  6 * 11
# case和control对比，红色是上调，蓝色是下调。
# 细胞互作数量对比网络图
par(mfrow = c(1,2))
pdf(file="4netVisual_circle.pdf",width=10,height=6)
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}
dev.off()
# save as Counts_Compare_net.pdf  6 * 9
# 比较 2D 空间中的主要来源和目标
# 比较2D空间中的incoming和outgoing交互强度，可以识别不同数据集之间显著变化的发送或接收信号的细胞群。
num.link <- sapply(cco.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(cco.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cco.list[[i]], title = names(cco.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf(file="5wrap_plots.pdf",width=10,height=6)
patchwork::wrap_plots(plots = gg)
dev.off()
# 3.2 指定细胞互作数量对比网络图----------------
table(cellchat@idents$PT)
table(cellchat@idents$M)
par(mfrow = c(1,1))
s.cell <- c("HSC", "Fibroblast","T_cell","Kupffer_cell")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
pdf(file="6netVisual_circle_MAC.pdf",width=10,height=6)
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[2]))
dev.off()
# 3.3 保守和特异性信号通路的识别与可视化-------------------
## 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c('#6ABDEF','#F38080'))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c('#6ABDEF','#F38080'))
p <- gg1 + gg2
p
ggsave("7Compare_pathway_strengh1.pdf", p, width = 14, height = 10, dpi = 300)
## 特定细胞通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, 
               sources.use = c("HSC"),
               targets.use = c("T_cell"),color.use = c('#6ABDEF','#F38080'))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, 
               sources.use = c("HSC"),
               targets.use = c("Kupffer_cell"),color.use = c('#6ABDEF','#F38080'))
p <- gg1 + gg2
p
ggsave("8Compare_pathway_strengh_Endo2MAC.pdf", p, width = 10, height = 6,dpi = 300)
# 3.4 流行学习识别差异信号通路-----------------
# 这里function的图出不来，只有structural的图可以出来
########################################################################
# cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")
# netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
# netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
# cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural")
# netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
# netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
# p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
# ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)
# 
# save(cellchat, file = "merge_cellchat_riskscore_0923.rdata")
######################################################################3
# 3.5 细胞信号模式对比------------------------
library(ComplexHeatmap)
# 总体信号模式对比
cco.list[[1]]@netP$pathways
cco.list[[2]]@netP$pathways
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)  ##⚠️️可定制通路⚠️（根据ranknet）
pdf(file="9netAnalysis_signalingRole_heatmap.pdf",width=10,height=12)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.list)[1],width = 6, height = 21)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.list)[2],width = 6, height = 21)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
write.csv(ht1@matrix,'signalingRole_heatmap_all_PT.csv')
write.csv(ht2@matrix,'signalingRole_heatmap_all_M.csv')
# save as Compare_signal_pattern_all.pdf  10*6
# ⚠️输出信号模式对比
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
pdf(file="10out_signalingRole_heatmap.pdf",width=10,height=12)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(cco.list)[1],width = 5, height = 21)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cco.list)[2],width = 5, height = 21)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
# save as Compare_signal_pattern_outgoing.pdf  15*15
write.csv(ht1@matrix,'signalingRole_heatmap_outgoing_PT.csv')
write.csv(ht2@matrix,'signalingRole_heatmap_outgoing_M.csv')

# ⚠️输入信号模式对比
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
pdf(file="11in_signalingRole_heatmap.pdf",width=10,height=12)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cco.list)[1],width = 5, height = 21)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cco.list)[2],width = 5, height = 21)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
write.csv(ht1@matrix,'signalingRole_heatmap_incoming_PT.csv')
write.csv(ht2@matrix,'signalingRole_heatmap_incoming_M.csv')
# save as Compare_signal_pattern_incoming.pdf  15*15
# 3.6 特定信号通路的对比--------------------
# 网络图
pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show) 
pdf(file="12netVisual_THBS_netplot.pdf",width=10,height=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(cco.list)[i]))
}
dev.off()
# save as Compare_IL16_net.pdf  10*6.5
# 热图
par(mfrow = c(1,2), xpd=TRUE)
pdf(file="13netVisual_heatmap_Endothelial_cell.pdf",width=10,height=6)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds", 
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()
write.csv(ht[[2]]@matrix,'prob_cco_high-MK.csv')
write.csv(ht[[1]]@matrix,'prob_cco_low-MK.csv')
# save as Compare_IL16_heatmap.pdf  12*6.5
# 和弦图
par(mfrow = c(1,1), xpd=TRUE)
pdf(file="14netVisual_JAM_aggregate2.pdf",width=12,height=10)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "chord", pt.title = 3, title.space = 0.05,
                      vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(cco.list)[i]))
}
dev.off()
# save as Compare_IL16_chord.pdf  10*6.5
# 3.7 配体-受体对比分析----------------
# 气泡图展示所有配体受体对的差异
# 我们可以比较由某些细胞群到其他细胞组的配体受体对调节的通信概率。这可以通过设置comparison在函数netVisual_bubble中来完成。
levels(cellchat@idents$joint)
#这里的source和target必须填2个
p <- netVisual_bubble(cellchat, sources.use = c("HSC","T_cell"), targets.use = c("HSC","T_cell"), 
                      comparison = c(1, 2), angle.x = 45,color.text = c('#277694','#BD1919'))
p
ggsave("15Compare_LR_bubble_CD4_epi.pdf", p, width = 6, height = 25,dpi = 300)
# 此外，我们可以在一个数据集中识别与另一个数据集相比，上升（增加）和下降调节（减少）信号配体受体对。
# 这可以通过指定max.dataset和min.dataset在函数netVisual_bubble中完成。信号增加意味着这些信号在一个
# 数据集中与其他数据集相比具有更高的通信概率（强度）。
# 气泡图展示上调或下调的配体受体对
pdf(file="16netVisual_LR_bubble_M.pdf",width=20,height=20)
p1 <- netVisual_bubble(cellchat, sources.use = c(1,4), targets.use = c(1,4), comparison = c(1, 2),color.text = c('#277694','#BD1919'), 
                       max.dataset = 2, title.name = "Increased signaling in M_group", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(1,4), targets.use = c(1,4), comparison = c(1, 2),color.text = c('#277694','#BD1919'), 
                       max.dataset = 1, title.name = "Decreased signaling in M_group", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
pc
dev.off()
write.csv(p1$data,'increased_signaling_in_M.csv')
write.csv(p2$data,'decreased_signaling_in_M.csv')
write.csv(pc$data,'merged_signaling_in_M.csv')

ggsave("17Compare_LR_regulated_M.pdf", pc, width = 15, height = 20,dpi = 300)


# NB：气泡图中显示的配体受体对可以通过⚠️ signaling.LSIncreased = gg1$data ⚠️访问。

#   通过比较每个 L-R 对和每对细胞组的两个数据集之间的通信概率，可以采用上述方法来识别上调和下调的信号。
# 另外，我们可以根据微分基因表达分析来识别上调和下调的信号配体对。具体来说，我们对每个细胞组执行两种
# 生物条件（即RS_Low和RS_High）之间的微分表达分析，然后根据发送者细胞中配体和接收器细胞中受体的折叠变化获得上调
# 和下调的信号。此类分析可如下所示。
#############################################################################
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
#RS_High <-read.csv(file="merged_signaling_in_RS_High.csv", header=TRUE, sep=",")
pos.dataset = "M"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.M <- subsetCommunication(cellchat, net = net, datasets = "M",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.PT <- subsetCommunication(cellchat, net = net, datasets = "PT",ligand.logFC = -0.1, receptor.logFC = -0.1)
# 由于信号基因在多亚单位中可能很复杂，我们可以使用net.upnet.down进一步的来获得单个信号基因
gene.M <- extractGeneSubsetFromPair(net.M, cellchat)
gene.PT <- extractGeneSubsetFromPair(net.PT, cellchat)
# 然后，我们使用气泡图或和弦图可视化上调和向下调的信号配体对。
pairRS.use.M = net.M[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairRS.use.M, sources.use = c(4), 
                        targets.use = c(1:10), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(cco.list)[2]))
#> Comparing communications on a merged object
pairRS.use.PT = net.PT[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairRS.use.PT, sources.use = c(4), 
                        targets.use = c(1:10), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(cco.list)[2]))
#> Comparing communications on a merged object
pc = gg1 + gg2
pc
ggsave("17Compare_LR_regulated_M_PT.pdf", pc, width = 15, height = 15,dpi = 300)
write.csv(gg1$data,'increased_signaling_in_M-gg1.csv')
write.csv(gg2$data,'decreased_signaling_in_M-gg2.csv')
write.csv(pc$data,'merged_signaling_in_M-gg.csv')
#pc
#ggsave("Compare_LR_regulated_signaling.png", pc, width = 12, height = 8,dpi = 240)
##############################################################################
# 使用和弦图可视化上调和下调的信号配体对
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf(file="18netVisual_chord_LR.pdf",width=20,height=15)
netVisual_chord_gene(cco.list[[2]], sources.use = c(1), targets.use = c(4), slot.name = 'net', net = net.M, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(cco.list)[2]))
#> Note: The first link end is drawn out of sector 'SPP1'.
netVisual_chord_gene(cco.list[[1]], sources.use = c(1), targets.use = c(4), slot.name = 'net', net = net.PT, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(cco.list)[2]))
dev.off()
##############################################################
# 第四部分：使用层次结构图、圆图或和弦图可视比较细胞-细胞通信----------------------
# 与单个数据集的 CellChat 分析类似，我们可以使用层次结构图、圆图或和弦图可视化细胞通信网络。
# 边缘颜色/重量、节点颜色/大小/形状：在所有可视化图中，边缘颜色与发送者源一致，边缘权重与交互强度成正比。较厚的边缘线表示信号更强。在层次结构图和圆图中，圆的大小与每个细胞组中的细胞数量成正比。在层次图中，实心和开放的圆分别代表源和目标。在和弦图中，内条颜色表示从相应的外条接收信号的目标。内条大小与目标接收的信号强度成正比。这种内条有助于解释复杂的和弦图。请注意，有一些内条没有任何和弦的一些细胞组，请忽略他们，因为这是一个本包尚未解决的问题。
cellchat@netP$PT[["pathways"]] #查看通路
pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
pdf(file="19netVisual_aggregate_MK.pdf",width=10,height=6)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}  ## 4.5* 9 
dev.off()

pathways.show <- c("SPP1") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}  
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
pdf(file="20ComplexHeatmap_MK.pdf",width=10,height=6)
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))  ##  4.5 * 9
dev.off()
# Chord diagram
pathways.show <- c("SPP1") 
par(mfrow = c(1,2), xpd=TRUE)
pdf(file="21netVisual_aggregate_MK.pdf",width=10,height=6)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cco.list)[i]))
}
dev.off()
#> Note: The first link end is drawn out of sector 'Inflam. FIB'. ## 50 * 100
# netVisual_chord_cell对于和弦图，CellChat 具有独立函数，通过调整circlize包中的不同参数来灵活可视化信号网络。例如，我们可以定义一个group命名的字符矢量，以创建多组和弦图，将细胞群集分组到不同的细胞类型。
########################################################################
# Chord diagram
group.CellType <- c(rep("1", 10), rep("2", 10), rep("3", 10)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.CellType) <- levels(cco.list[[1]]@idents)
pathways.show <- c("SPP1") 
par(mfrow = c(1,2), xpd=TRUE)
pdf(file="22netVisual_chord_cell.pdf",width=10,height=6)
for (i in 1:length(cco.list)) {
  netVisual_chord_cell(cco.list[[i]], signaling = pathways.show, group = group.CellType, title.name = paste0(pathways.show, " signaling network - ", names(cco.list)[i]))
}
dev.off()
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
#####################################################################
#   使用和弦图，CellChat 提供两个函数netVisual_chord_cell和netVisual_chord_gene，可视化具有不同目的和
# 不同级别的细胞-细胞通信。 netVisual_chord_cell用于可视化不同细胞群之间的细胞-细胞通信（和弦图中的
# 每个部分是细胞组），netVisual_chord_gene用于可视化由多个配体受体或信号通路调解的细胞-细胞通信（和弦图中的每个部分都是配体、受体或信号通路）。
table(cco.list$PT@idents)
par(mfrow = c(1, 1), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
pdf(file="23netVisual_chord_Fib_CD8.pdf",width=15,height=10)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]], sources.use = 1, targets.use = c(3), lab.cex = 0.5, title.name = paste0("Signaling from Crypt_cell - ", names(cco.list)[i]))
}
dev.off()
# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 1), xpd=TRUE)
pdf(file="24netVisual_chord_multi.pdf",width=20,height=15)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]], sources.use = c(1,3, 4), targets.use = c(1,3,4),  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(cco.list)[i]), legend.pos.x = 10)
}
dev.off()
# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 1), xpd=TRUE)
pdf(file="25netVisual_chord_multi.pdf",width=20,height=15)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]], sources.use = c(1,2,3,4), targets.use = c(5:11),slot.name = "netP", title.name = paste0("Signaling pathways sending from multi - ", names(cco.list)[i]), legend.pos.x = 10)
}
dev.off()
#> Note: The second link end is drawn out of sector ' '.
#> Note: The first link end is drawn out of sector 'SPP1'.
#> Note: The second link end is drawn out of sector ' '.
#> Note: The first link end is drawn out of sector 'SPP1 '.
# NB：在生成绘图时，请忽略注释，例如"Note: The first link end is drawn out of sector ‘SPP1’"。如果基因名称重叠，您可以通过降低small.gap值来调整参数。
# 第五部分：比较不同数据集之间的信号基因表达分布-------------------
# 我们可以利用seurat包装的函数plotGeneExpression绘制与L-R对或信号通路相关的信号基因的基因表达分布图。
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("PT", "M")) # set factor level
pdf(file="26pathway_colplot.pdf",width=10,height=6)
plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", colors.ggplot = T, color.use = c('#68CCCA','#F38888'))
dev.off()
###########################################