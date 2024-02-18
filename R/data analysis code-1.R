###  data log  ###
data <- read.csv("WTA_data.csv",stringsAsFactors = F,row.names = 1)
data <- log2(data+1)
#write.csv(data,"WTA_datalog2.csv",quote = F)



#######################  ANPY cluster   #########################
library(pheatmap)
gene <- c('ASCL1','NEUROD1','POU2F3','YAP1')
data <- read.csv("WTA_datalog2.csv",
                 stringsAsFactors = F,row.names = 1)
data <- data[gene,]
red <- rgb(255,0,0,maxColorValue = 255)
blue <- "#13334c";
ye <- "#F5E0A3";
linshi <- apply(data,2,scale)
rownames(linshi) <- rownames(data)
hist(linshi)
out <- pheatmap(linshi,fontsize=6,cutree_cols = 3,
                color  = colorRampPalette(c(blue,ye,red))(100),
                #annotation_col = annotation_col,
                #annotation_colors = ann_colors,
                clustering_method = "ward.D",
                border_color = "grey60",
                cluster_cols = T, cluster_rows = F,
                show_rownames = T, show_colnames = T)


###  ROI NE subtype  ###
gene <- read.csv("NE_gene.csv",stringsAsFactors = F,row.names = 1)
data <- read.csv("WTA_multi_ROIdatalog2.csv",stringsAsFactors = F,row.names = 1)
intergene <- intersect(rownames(gene),rownames(data))
data <- data[intergene,]
gene <- gene[intergene,]
score <- c()
for (i in 1:dim(data)[2]){
  ne <- cor.test(data[,i],gene$NE.class.mean.expression)
  non_ne <- cor.test(data[,i],gene$Non.NE.class.mean.expression)
  score <- c(score,(ne$estimate - non_ne$estimate)/2)
}
names(score) <- colnames(data)

###  NE score  heatmap   ### 
library(ggplot2)
library(pheatmap)
gene <- read.csv("NE_gene.csv",stringsAsFactors = F,row.names = 1)
data <- read.csv("WTA_multi_ROIdatalog2.csv",stringsAsFactors = F,row.names = 1)
intergene <- intersect(rownames(gene),rownames(data))
data <- data[intergene,]
ne <- read.csv("NE_subtype.csv",stringsAsFactors = F,row.names = 1)
ne <- ne[order(ne$label),]
data <- data[,rownames(ne)]

annotation_col = data.frame( Cluster = as.factor(ne$label))
rownames(annotation_col) = colnames(data)
ann_colors = list(Cluster = c('NE_high' = "#F38181","NE_low" = "#95E1D3"))  
red <- rgb(255,0,0,maxColorValue = 255)
blue <- rgb(0,0,255,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)
linshi <- apply(data,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(data)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(-2)
out <- pheatmap(linshi,fontsize=6,
                annotation_col = annotation_col,
                annotation_colors = ann_colors,
                clustering_method = "ward.D2",
                border_color = "grey60",
                cluster_cols = F, cluster_rows = F,
                show_rownames = T, show_colnames = T)



















########################   垃圾  ！！！  ############################
###  NE 基因 和  ANPY 基因聚类  ###
setwd("D:\\北京_小细胞肺癌Yap1\\尝试一")
data <- read.csv("D:\\北京_小细胞肺癌\\zzc处理后的数据\\29例WTA_log2原始数据.csv",
                 stringsAsFactors = F,row.names = 1)
colnames(data) <- substr(colnames(data),2,19)
gene <- read.csv("D:\\北京_小细胞肺癌Yap1\\zzc处理后的数据\\NE_gene.csv",
                 stringsAsFactors = F,row.names = 1)
gene <- rownames(gene)
gene <- c(gene,'ASCL1','NEUROD1','POU2F3','YAP1','INSM1')
intergene <- intersect(gene,rownames(data))
data <- data[intergene,]

label <- read.csv("病人亚型label.csv",
                  stringsAsFactors = F,row.names = 1)
label <- label[colnames(data),]
annotation_col = data.frame( Cluster = as.factor(label$multi_label))
rownames(annotation_col) = rownames(label)
ann_colors = list(Cluster = c('NEhigh-A' = "#de9974","NEhigh-ANPmixed" = "#d74b4b",
                              "NElow" = "#293462","NEmixed" = "#1f640a",
                              "NEhigh-N" = "#39BAE8"))  

red <- rgb(255,0,0,maxColorValue = 255)
blue <- rgb(0,0,255,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)
linshi <- apply(data,2,scale)
rownames(linshi) <- rownames(data)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(2)
out <- pheatmap(linshi,fontsize=6,cutree_cols = 2,cutree_rows = 2,
                color  = colorRampPalette(c(blue,white,red))(100),
                annotation_col = annotation_col,
                annotation_colors = ann_colors,
                clustering_method = "ward.D",
                border_color = "grey60",
                cluster_cols = T, cluster_rows = T,
                show_rownames = T, show_colnames = T
)  #聚类热图
out.dist=dist(t(linshi)) 
out.hclust=hclust(out.dist,method="ward.D")
out.id=cutree(out.hclust,k=2)  
#write.csv(out.id,  "111.csv",quote = F)
clinical <- read.csv("D:\\北京_小细胞肺癌Yap1\\zzc处理后的数据\\29临床病理信息表.csv",
                     stringsAsFactors = F,row.names = 1)
label <- read.csv("111.csv",stringsAsFactors = F,row.names = 1)
rownames(label) <- substr(rownames(label),2,22)
clinical <- clinical[rownames(label),]
clinical$label <- label$x
surv <- survfit(Surv(as.numeric(DFS),DFSState)~label,data = clinical)
surv
survdiff(Surv(as.numeric(DFS),DFSState)~label,data = clinical)
summary(surv)
summary(coxph(Surv(as.numeric(DFS),DFSState)~label,data = clinical))
ggsurvplot(surv,
           pval = TRUE,#添加P值
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette = c("#DE9974","#D74B4B",'#293462',
                       "#1f640a"),
           legend = c(2,0.5), # 指定图例位置
           legend.title = "", # 设置图例标题
           conf.int = F,#添加置信区间
           #legend.labs = c("A-NEhigh","Mix-NEhigh","NElow","NEmixed"), # 指定图例分组标签
           xlab = "Time (months)",
           ylab = "Disease-free survival"
)
surv <- survfit(Surv(as.numeric(OS),OSState)~label,data = clinical)
surv
survdiff(Surv(as.numeric(OS),OSState)~label,data = clinical)
summary(surv)
summary(coxph(Surv(as.numeric(OS),OSState)~label,data = clinical))
ggsurvplot(surv,
           pval = TRUE,#添加P值
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette = c("#DE9974","#D74B4B",'#293462',
                       "#1f640a"),
           legend = c(2,0.5), # 指定图例位置
           legend.title = "", # 设置图例标题
           conf.int = F,#添加置信区间
           #legend.labs = c("A-NEhigh","Mix-NEhigh","NElow","NEmixed"), # 指定图例分组标签
           xlab = "Time (months)",
           ylab = "Overall survival"
)

a <- data.frame(y = as.numeric(data["INSM1",]),
                x = as.factor(label$x))
ggplot(data = a , aes(x = x , y = y,fill=x)) + 
  stat_compare_means(label.y = 0.6, label.x = 1.5)  +  ### 添加总体的检验
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_jitter(aes(fill=x),width =0.2,shape = 21,size=3.5)+
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF",alpha = 0.4)+ #绘制箱线图
  scale_fill_manual(values = c("#D94E48","#5175A4"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 0), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=10,colour="black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 10,colour="black"), #设置y轴标题的字体属性
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=10),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("INSM1")

























