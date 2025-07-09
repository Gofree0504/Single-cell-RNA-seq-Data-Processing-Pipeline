## 计算出每个样本中每个细胞类型的百分比--------------------
#install.packages("ggstatsplot")
#install.packages("ratantools")
#install.packages("afex")
library(ggstatsplot)
#查看各组细胞数
table(sce.all$orig.ident)
prop.table(table(sce.all@meta.data$celltype))
table(sce.all@meta.data$celltype, sce.all$treated)
#绘制堆叠柱状图
Cellratio <- prop.table(table((sce.all@meta.data$celltype), sce.all$treated), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
pdf(file="14.stackplot.pdf",width=10,height=6)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
dev.off()
################################################
#画cluster中高表达基因
# 计算每个聚类的平均表达量
sce.all.markers <- FindAllMarkers(sce.all, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
# 找到每个聚类中表达量最高的前10个基因
top10 <- sce.all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf(file="15.1.Doheatmap.pdf",width=10,height=6)
DoHeatmap(sce.all, features = top10$gene)
dev.off()