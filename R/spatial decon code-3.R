##  https://github.com/Nanostring-Biostats/SpatialDecon

library(SpatialDecon)
library(GeomxTools)
library(Seurat)
cell_count_matrix <- read.csv("SCLC_cell_count_matrix.csv", row.names = 1, header = TRUE)
cell_info <- read.csv("SCLC_cell_info.csv", row.names = 1, header = TRUE)
cell_info <- cell_info[colnames(cell_count_matrix),]
seurat_object <- CreateSeuratObject(counts = cell_count_matrix)
seurat_object <- AddMetaData(seurat_object, metadata = cell_info)
annots <- data.frame(CellID = rownames(cell_info),
                     LabeledCellType = cell_info$cell_type)
custom_mtx <- create_profile_matrix(mtx = cell_count_matrix,            # cell x gene count matrix
                                    cellAnnots = annots,  # cell annotations with cell type and cell name as columns 
                                    cellTypeCol = "LabeledCellType",  # column containing cell type
                                    cellNameCol = "CellID",           # column containing cell ID/name
                                    matrixName = "custom_mini_colon", # name of final profile matrix
                                    outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
                                    normalize = FALSE,                # Should data be normalized? 
                                    minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
                                    minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                    scalingFactor = 5,                # what should all values be multiplied by for final matrix
                                    discardCellTypes = TRUE)          # should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.

data <- read.csv("WTA_multi_ROIdata.csv",stringsAsFactors = F,row.names = 1)
name <- intersect(rownames(custom_mtx),rownames(data))
data <- as.matrix(data[name,]) ##  CIBERSORTx  input
custom_mtx <- as.matrix(custom_mtx[name,])  ##  CIBERSORTx  input
#write.csv(data,"wta_CIBERSORTx.csv",quote = F)
#write.csv(custom_mtx,"single_CIBERSORTx.csv",quote = F)


###  CIBERSORTx result  ###
library(ggplot2)
library(ggpubr)
library(gghalves)
data <- read.csv("CIBERSORTx_Job5_Results.csv",
                 stringsAsFactors = F,row.names = 1)
cluster <- read.csv("cluster_roi.csv",
                    stringsAsFactors = F,row.names = 1)
cluster <- cluster[rownames(data),]
data$cluster <- cluster$cluster_new
data$cluster <- factor(data$cluster,
                       levels = c("highest-heterogeneity",
                                  "middle-heterogeneity",
                                  "low-heterogeneity"))
a1 <- ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=cluster, y=Tcell, fill = cluster),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=cluster, y=Tcell, color=cluster), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  scale_color_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  stat_compare_means(data = data, 
                     aes(x=cluster, y=Tcell),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="Tcell")
a2 <- ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=cluster, y=Bcell, fill = cluster),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=cluster, y=Bcell, color=cluster), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  scale_color_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  stat_compare_means(data = data, 
                     aes(x=cluster, y=Bcell),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="Bcell")
a3 <- ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=cluster, y=Mast, fill = cluster),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=cluster, y=Mast, color=cluster), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  scale_color_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  stat_compare_means(data = data, 
                     aes(x=cluster, y=Mast),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="Mast")
a4 <- ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=cluster, y=Mye, fill = cluster),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=cluster, y=Mye, color=cluster), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  scale_color_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  stat_compare_means(data = data, 
                     aes(x=cluster, y=Mye),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="Mye")
a5 <- ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=cluster, y=tumor, fill = cluster),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=cluster, y=tumor, color=cluster), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  scale_color_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  stat_compare_means(data = data, 
                     aes(x=cluster, y=tumor),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="tumor")
a6 <- ggplot() +
  geom_half_boxplot(data = data, 
                    aes(x=cluster, y=normal, fill = cluster),
                    width = 0.4, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = data, 
                  aes(x=cluster, y=normal, color=cluster), 
                  shape = 16, alpha=1, size = 0.8)+
  scale_fill_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  scale_color_manual(values = c('#E59572','#D3BD74','#9BAEC8'))+
  stat_compare_means(data = data, 
                     aes(x=cluster, y=normal),
                     method = "kruskal.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.6),
        legend.position="none")+ labs(x='', y="normal")

data_new1 <- data[which(data$cluster != "low-heterogeneity"),]
data_new2 <- data[which(data$cluster != "middle-heterogeneity"),]
data_new3 <- data[which(data$cluster != "highest-heterogeneity"),]
wilcox.test(data_new1$Tcell~data_new1$cluster)
wilcox.test(data_new2$Tcell~data_new2$cluster)
wilcox.test(data_new3$Tcell~data_new3$cluster)








