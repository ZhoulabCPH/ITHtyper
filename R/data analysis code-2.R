
######  cv score and cluster   ######
library(dplyr)
data <- read.csv("WTA_multi_ROIdata.csv",stringsAsFactors = F,row.names = 1)
# calculate  CV
cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}
cv_result1 <- apply(data, 1, cal_cv) 
cv_result1 <- cv_result1[cv_result1 %>% order %>% rev]
#write.csv(cv_result1,"gene_cv.csv",quote = F)
cv_result2 <- apply(data, 2, cal_cv) # 
cv_result2 <- cv_result2[cv_result2 %>% order %>% rev]
#write.csv(cv_result2,"ROI_cv.csv",quote = F)


#####  
library(pheatmap)
data <- read.csv("WTA_multi_ROIdata.csv",stringsAsFactors = F,row.names = 1)
clinical <- read.csv("clinical.csv",stringsAsFactors = F,row.names = 1)
roi <- read.csv("roi_code.csv",stringsAsFactors = F,row.names = 1)
gene <- read.csv("gene_cv.csv",stringsAsFactors = F,row.names = 1)
gene_top200 <- rownames(gene)[1:200]
red <- "#D94E48";
blue <- "#5175A4";
white <- rgb(255,255,255,maxColorValue = 255)
linshi <- apply(data,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(data)
rownames(linshi) <- rownames(data)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(-2)
roi1 <- roi[colnames(linshi),]
names(roi1) <- colnames(linshi)
clinical <- clinical[roi1,]
annotation_col <- data.frame(patient = as.factor(roi1),
                             Gender = clinical$Gender,
                             SmokingHistory = clinical$SmokingHistory,
                             TumorLocation = clinical$TumorLocation,
                             PCI = clinical$PCI,
                             Tstage = as.factor(clinical$Tstage),
                             Nstage = as.factor(clinical$Nstage),
                             Mstage = as.factor(clinical$Mstage))
rownames(annotation_col) <- names(roi1)
ann_colors = list(patient = c("Pt-100"="#996699","Pt-115"="#006699","Pt-125"="#58B4AB","Pt-162"="#E96463",
                              "Pt-177"="#7A9D96","Pt-181"="#8b220d","Pt-185"="#c89c0e","Pt-191"="#B89582",
                              "Pt-192"="#492711","Pt-214"="#FFCC99","Pt-227"="#cbe86b","Pt-237"="#bf5704",
                              "Pt-241"="#E39183","Pt-245"="#5B6A27","Pt-269"="#EEDCC6","Pt-274"="#ff85cb",
                              "Pt-278"="#53bbf4","Pt-289"="#C5C6B6","Pt-294"="#ffc33c","Pt-323"="#dee2d1",
                              "Pt-332"="#f9fbba","Pt-78"="#d4edf4","Pt-80"="#57D1C9","Pt-88"="#AC8697",
                              "Pt-97"="#9baec8"),
                  Gender = c("Male" = "#75c0b8", "Female" = "#e09a7d"),
                  SmokingHistory = c("No" = "#45a3cb","Yes" = "#edab4c"),
                  TumorLocation = c("Left" = "#f3cfcb","Right" = "#9d8fba"),
                  PCI = c("No" = "#eae1e1","Yes" = "#778bbe"),
                  Tstage = c("1"="#e8d2c5","2"="#a1918f","3"="#716c75"),
                  Nstage = c("0"="#cccccc","1"="#91c0c0",
                             "2"="#c29e2e","3"="#647370"),
                  Mstage = c("0"="#d0e7ee"))
linshi200 <- linshi[gene_top200,]
out <- pheatmap(linshi200,fontsize=6,cutree_cols = 3,
                color  = colorRampPalette(c(blue,white,red))(100),
                annotation_col = annotation_col,
                annotation_colors = ann_colors,
                clustering_method = "ward.D2",
                border_color = "grey60",
                cluster_cols = T, cluster_rows = T,
                show_rownames = F, show_colnames = T
)  
out.dist=dist(t(linshi200)) 
out.hclust=hclust(out.dist,method="ward.D2")
name <- out.hclust$labels[out.hclust$order]
out.id <- cutree(out.hclust,k=3)
#write.csv(out.id,"cluster_roi.csv",quote = F)
###  ITH score PMID: 32917965 ###
#if (!requireNamespace("devtools", quietly = TRUE))
#  install.packages("devtools")
#devtools::install_github("WangX-Lab/DEPTH")
library(DEPTH)
library(ggplot2)
library(ggpubr)
data <- read.csv("WTA_multi_ROIdata.csv",stringsAsFactors = F,row.names = 1)
type <- data.frame(State = colnames(data),
                   Identification = rep("Tumor",dim(data)[2]))
ith_score <- DEPTH(data, type) 
rownames(ith_score) <- ith_score[,1]
ith_score <- as.data.frame(ith_score)
ith_score$`ITH score` <- as.numeric(ith_score$`ITH score`)
ith_score <- ith_score[name,]
barplot(ith_score$`ITH score`)
out.id <- out.id[rownames(ith_score)]
ith_score$cluster <- factor(out.id,levels = c(2,3,1)) ## 2:H-H, 3:M-H, 1:L=H
ith_score$cluster <- factor(ith_score$cluster,
                            levels = c("2","3","1"))
my_com <- list(c("2","3"),c("2","1"),c("3","1"))
ggplot(data = ith_score, aes(x = cluster , y = `ITH score`,fill=cluster)) + 
  stat_compare_means(label.y = 3, label.x = 1.5)  +  ###
  stat_compare_means(comparisons = my_com,label = "p.signif")  +
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_jitter(aes(fill=cluster),width =0.2,shape = 21,size=3.5)+
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF",alpha = 0.4)+ 
  scale_fill_manual(values = c("#9baec8","#F2BDD0","#D3BD74"))+ 
  theme_bw()+
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 0), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("ith_score")
###  ConsensusClusterPlus validation  ###
library(ConsensusClusterPlus)
library(ggplot2)
library(ggpubr)
data <- read.csv("WTA_multi_ROIdata.csv",stringsAsFactors = F,row.names = 1)
clinical <- read.csv("clinical.csv",stringsAsFactors = F,row.names = 1)
roi <- read.csv("roi_id.csv",stringsAsFactors = F,row.names = 1)
gene <- read.csv("gene_cv.csv",stringsAsFactors = F,row.names = 1)
gene_top200 <- rownames(gene)[1:200]
red <- "#D94E48";
blue <- "#5175A4";
white <- rgb(255,255,255,maxColorValue = 255)
linshi <- apply(data,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(data)
rownames(linshi) <- rownames(data)
linshi200 <- linshi[gene_top200,]
results <- ConsensusClusterPlus(linshi200, maxK = 10,
                                reps = 1000, pItem = 1,
                                pFeature = 1,  
                                clusterAlg = "km", 
                                distance="euclidean",
                                plot = "pdf")  
## best k
rivzi2016_PAC=c()  
for(i in 2:10){
  cM=results[[i]]$consensusMatrix
  Fn=ecdf(cM[lower.tri(cM)])  #
  rivzi2016_PAC=c(rivzi2016_PAC, Fn(0.9) - Fn(0.1))  #
}#end for i# The optimal K
rivzi2016_optK=c(2:10)[which.min(rivzi2016_PAC)]
label <- results[[3]][['consensusClass']]
#write.csv(label,"consensusCluster_label.csv",quote = F)
PAC_data <- data.frame(x = 2:10,
                       pac = rivzi2016_PAC)
ggplot(PAC_data, aes(x = x, y = pac,color = "#4F6F52")) + 
  geom_point() + geom_line() + theme_bw() + 
  scale_color_identity() + labs(title = "Consensus Cluster")
###  compare  ###
cluster <- read.csv("cluster_roi.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
ith_score <- read.csv("ROI_ITH_score.csv",
                      stringsAsFactors = F,row.names = 1,check.names = F)
ith_score <- ith_score[rownames(ith_score),]
cluster$ith_score <- ith_score$`ITH score`
ggplot(data = cluster, aes(x = consensusCluster , y = ith_score,
                           fill=consensusCluster)) + 
  stat_compare_means(label.y = 3, label.x = 1.5)  +  
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_jitter(aes(fill=consensusCluster),width =0.2,shape = 21,size=3.5)+
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF",alpha = 0.4)+ 
  scale_fill_manual(values = c("#9baec8","#F2BDD0","#D3BD74"))+ 
  theme_bw()+
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 0), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("ith_score")
table(paste(cluster$consensusCluster,cluster$cluster_new))



###  ROI   TSNE  ###
library(factoextra)
gene <- read.csv("gene_cv.csv",stringsAsFactors = F,row.names = 1)
gene_top200 <- rownames(gene)[1:200]
data <- read.csv("WTA_multi_ROIdata.csv",stringsAsFactors = F,row.names = 1)
label <- read.csv("cluster_roi.csv",stringsAsFactors = F,row.names = 1)
data <- t(data)
data_new <- apply(data, 2, scale)
rownames(data_new) <- rownames(data)
label <- label[rownames(data_new),]
library(Rtsne)
tsne_out = Rtsne(data_new,dims = 2,
                 pca = T,max_iter = 1000,
                 theta = 0.4,perplexity = 20,
                 verbose = F )
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) = c("tSNE1","tSNE2")
tsne_result$patient <- label$patient
tsne_result$cluster <- label$cluster_new
ggplot(tsne_result,aes(tSNE1,tSNE2,color=patient)) +
  geom_point(size = 3)+
  scale_color_manual(values = c("#996699","#006699","#58B4AB","#E96463",
                                "#7A9D96","#8b220d","#c89c0e","#B89582",
                                "#492711","#FFCC99","#cbe86b","#bf5704",
                                "#E39183","#5B6A27","#EEDCC6","#ff85cb",
                                "#53bbf4","#C5C6B6","#ffc33c","#dee2d1",
                                "#f9fbba","#d4edf4","#57D1C9","#AC8697",
                                "#9baec8"))
ggplot(tsne_result,aes(tSNE1,tSNE2,color=cluster)) +
  geom_point(size = 3)+
  scale_color_manual(values = c("#e59572", "#9baec8", "#d3bd74"))


###   ROI cv boxplot  ###
library(ggplot2)
library(ggpubr)
cv <- read.csv("ROI_cv.csv",stringsAsFactors = F,row.names = 1)
cluster <- read.csv("cluster_roi.csv",stringsAsFactors = F,row.names = 1)
cv <- cv[rownames(cluster),]
cv$cluster <- cluster$cluster_new
cv$cluster <- factor(cv$cluster,
                     levels = c("highest-heterogeneity",
                                "middle-heterogeneity",
                                "low-heterogeneity"))
my_com <- list(c("highest-heterogeneity","middle-heterogeneity"),
               c("highest-heterogeneity","low-heterogeneity"),
               c("middle-heterogeneity","low-heterogeneity"))
ggplot(data = cv, aes(x = cluster , y = CV,fill=cluster)) + 
  stat_compare_means(label.y = 0.3, label.x = 1.5)  +  
  stat_compare_means(comparisons = my_com,label = "p.signif")  +
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_jitter(aes(fill=cluster),width =0.2,shape = 21,size=3.5)+ #绘制箱线图
  scale_fill_manual(values = c("#F2BDD0","#9baec8","#D3BD74"))+ #
  theme_bw()+ #背景变为白色
  #ylim(0,0.01) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 0), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black",
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("CV")


###  ROI SPD 
library(ggplot2)
library(ggpubr)
data <- read.csv("WTA_ROI_x_y.csv",stringsAsFactors = F,row.names = 1)
label <- read.csv("cluster_roi.csv",stringsAsFactors = F,row.names = 1)
data <- data[rownames(label),]
patient <- unique(label$patient)
result <- matrix(,,4)
for (i in 1:length(patient)){
  weizhi <- which(patient[i] == label$patient)
  new_data <- data[rownames(label)[weizhi],]
  if (dim(new_data)[1] > 2){
    for (j in 1:(dim(new_data)[1]-1)){
      for (z in (j+1):(dim(new_data)[1])){
        a <- sqrt((new_data[j,1]-new_data[z,1])^2 + (new_data[j,2]-new_data[z,2])^2)
        b <- c(a,patient[i],label$cluster_new[weizhi][1],paste(rownames(new_data)[j],rownames(new_data)[z],sep = "_"))
        result <- rbind(result,b)
      }
    }
  }
  else if (dim(new_data)[1] == 2) {
    a <- sqrt((new_data[1,1]-new_data[2,1])^2 + (new_data[1,2]-new_data[2,2])^2)
    b <- c(a,patient[i],label$cluster_new[weizhi][1],paste(rownames(new_data)[1],rownames(new_data)[2],sep = "_"))
    result <- rbind(result,b)
  }
}
result <- result[-1,]
result <- as.data.frame(result)
colnames(result) <- c("Physical_distance","patient","label","ROI")
result$Physical_distance <- as.numeric(result$Physical_distance)
#write.csv(result,"ROI_SPD.csv",quote = F)
ggplot(data = result, aes(x = label , y = Physical_distance,fill=label)) + 
  stat_compare_means(label.y = 4000, label.x = 1.5)  +  
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_jitter(aes(fill=label),width =0.2,shape = 21,size=3.5)+
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF",alpha = 0.4)+ #绘制箱线图
  scale_fill_manual(values = c("#e59572","#9baec8","#D3BD74"))+ 
  theme_bw()+ 
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 0), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  #
                                 size=10),
        legend.title=element_text(colour="black", #
                                  size=10),
        panel.grid.major = element_blank(),   #
        panel.grid.minor = element_blank())+  #
  ylab("Physical_distance")
data <- read.csv("ROI_SPD.csv",
                 stringsAsFactors = F,row.names = 1)
ggplot(data, aes(Physical_distance)) +
  geom_histogram()


### C-score and SPD ###
library(ggplot2)
library(ggpubr)
cor <- read.csv("C-score.csv",stringsAsFactors = F,
                row.names = 1,check.names = F)
distan <- read.csv("ROI_SPD.csv",stringsAsFactors = F,
                   row.names = 1,check.names = F)
cor <- cor[rownames(distan),]
data_new <- data.frame(cor = cor$rho,
                       distan = distan$Physical_distance)
ggplot(data = data_new, aes(x = cor, y = distan)) + 
  geom_point(colour = "#a696c8", size = 1) +  
  geom_smooth(method = lm,colour='#fd5f00',fill='#E7E1D7') + 
  ylab("Physical_distance") + xlab("C-score") + theme_bw()+ 
  stat_cor(method = "spearman", label.x = median(data_new$cv),
           label.y = median(data_new$ith)) + theme_bw()



### specific gene in ROI  ###
library(DESeq2)
data <- read.csv("WTA_multi_ROIdata.csv",
                 stringsAsFactors = F,row.names = 1,check.names = F)
label <- read.csv("cluster_roi.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
data <- as.data.frame(data[,rownames(label)])
##  l-ITH
cluster1 <- as.factor(label$label1)
coldata <- data.frame(condition = factor(label$label1))
dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design= ~condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('condition', 'LH', 'HH_MH'))
res <- as.data.frame(res)
data1 <- res
data1 <- data1[order(data1$log2FoldChange),]
data1$padj <- p.adjust(data1$pvalue,method = "fdr")
#write.csv(data1,"l-ITH-gene.csv",quote = F)
## m-ITH
cluster1 <- as.factor(label$label2)
coldata <- data.frame(condition = factor(label$label2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design= ~condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('condition', 'MH', 'LH_HH'))
res <- as.data.frame(res)
data2 <- res
data2 <- data2[order(data2$log2FoldChange),]
data2$padj <- p.adjust(data2$pvalue,method = "fdr")
#write.csv(data2,"m-ITH-gene.csv",quote = F)
## h-ITH
cluster1 <- as.factor(label$label3)
coldata <- data.frame(condition = factor(label$label3))
dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design= ~condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('condition', 'HH', 'LH_MH'))
res <- as.data.frame(res)
data3 <- res
data3 <- data3[order(data3$log2FoldChange),]
data3$padj <- p.adjust(data3$pvalue,method = "fdr")
#write.csv(data3,"m-ITH-gene.csv",quote = F)

###  vocal plot  ###
library(ggplot2)
library(ggrepel)
library(ggpubr)
hh <- read.csv("h-ITH-gene.csv",stringsAsFactors = F,row.names = 1)
mh <- read.csv("m-ITH-gene.csv",stringsAsFactors = F,row.names = 1)
lh <- read.csv("l-ITH-gene.csv",stringsAsFactors = F,row.names = 1)
FC = 1.5
P.Value = 0.05
hh$sig[(-1*log10(hh$pvalue) < -1*log10(P.Value)|hh$pvalue=="NA")|(hh$log2FoldChange < log2(FC))& hh$log2FoldChange  > -log2(FC)] <- "NotSig"
hh$sig[-1*log10(hh$pvalue) >= -1*log10(P.Value) & hh$log2FoldChange >= log2(FC)] <- "Up"
hh$sig[-1*log10(hh$pvalue) >= -1*log10(P.Value) & hh$log2FoldChange <= -log2(FC)] <- "Down"
hh$label <- rownames(hh)
hh$label[11:18666] <- ""
hh$logp <- (-1*log10(hh$pvalue))
ggplot(hh,aes(logp,hh$log2FoldChange)) +    
  geom_point(aes(color = sig)) +                           
  labs(title="H-H",x="", y="log2(Fold change)") + 
  scale_color_manual(values = c("#e59572","#E6E6E6","#e59572")) + 
  geom_hline(yintercept=c(-log2(FC),log2(FC)),linetype=2) +
  ggrepel::geom_label_repel(aes(label = label),data = hh,color="black")
mh$sig[(-1*log10(mh$pvalue) < -1*log10(P.Value)|mh$pvalue=="NA")|(mh$log2FoldChange < log2(FC))& mh$log2FoldChange  > -log2(FC)] <- "NotSig"
mh$sig[-1*log10(mh$pvalue) >= -1*log10(P.Value) & mh$log2FoldChange >= log2(FC)] <- "Up"
mh$sig[-1*log10(mh$pvalue) >= -1*log10(P.Value) & mh$log2FoldChange <= -log2(FC)] <- "Down"
mh$label <- rownames(mh)
mh$label[11:18666] <- ""
mh$logp <- (-1*log10(mh$pvalue))
ggplot(mh,aes(logp,mh$log2FoldChange)) +   
  geom_point(aes(color = sig)) +                          
  labs(title="M-H",x="", y="log2(Fold change)") + 
  scale_color_manual(values = c("#d3bd74","#E6E6E6","#d3bd74")) + 
  geom_hline(yintercept=c(-log2(FC),log2(FC)),linetype=2) +
  ggrepel::geom_label_repel(aes(label = label),data = mh,color="black",size = 3)
lh$sig[(-1*log10(lh$pvalue) < -1*log10(P.Value)|lh$pvalue=="NA")|(lh$log2FoldChange < log2(FC))& lh$log2FoldChange  > -log2(FC)] <- "NotSig"
lh$sig[-1*log10(lh$pvalue) >= -1*log10(P.Value) & lh$log2FoldChange >= log2(FC)] <- "Up"
lh$sig[-1*log10(lh$pvalue) >= -1*log10(P.Value) & lh$log2FoldChange <= -log2(FC)] <- "Down"
lh$label <- rownames(lh)
lh$label[11:18666] <- ""
lh$logp <- (-1*log10(lh$pvalue))
ggplot(lh,aes(logp,lh$log2FoldChange)) +    
  geom_point(aes(color = sig)) +                          
  labs(title="L-H",x="", y="log2(Fold change)") + 
  scale_color_manual(values = c("#9baec8","#E6E6E6","#9baec8")) + 
  geom_hline(yintercept=c(-log2(FC),log2(FC)),linetype=2) +
  ggrepel::geom_label_repel(aes(label = label),data = lh,color="black",size = 3)


###  gene network  ###
library(Hmisc)
library(igraph)
library(ggraph)
data <- read.csv("WTA_multi_ROIdata.csv",stringsAsFactors = F,row.names = 1)
hh <- read.csv("h-ITH-gene.csv",stringsAsFactors = F,row.names = 1)
mh <- read.csv("m-ITH-gene.csv",stringsAsFactors = F,row.names = 1)
lh <- read.csv("l-ITH-gene.csv",stringsAsFactors = F,row.names = 1)
all_up_fc <- rbind(hh[which(hh$sig == "Up"),],
                   mh[which(mh$sig == "Up"),],
                   lh[which(lh$sig == "Up"),])
data.hh <- data[rownames(hh)[which(hh$sig == "Up")],]
data.mh <- data[rownames(mh)[which(mh$sig == "Up")],]
data.lh <- data[rownames(lh)[which(lh$sig == "Up")],]
data.all <- t(rbind(data.hh,data.mh,data.lh))
res <- rcorr(data.all,type=c("spearman"))
res.r <- res$r
res.p <- res$P
mat.r <- matrix(,30450,3)
mat.p <- matrix(,30450,3)
for (i in 1:175){
  for (j in 1:175){
    if (i != j){
      a <- c(rownames(res.r)[i],colnames(res.r)[j],res.r[i,j])
      mat.r <- rbind(mat.r,a)
    }
  }
}
mat.r <- mat.r[-1,]
for (i in 1:175){
  for (j in 1:175){
    if (i != j){
      a <- c(rownames(res.p)[i],colnames(res.p)[j],res.p[i,j])
      mat.p <- rbind(mat.p,a)
    }
  }
}

mat.p <- mat.p[-1,]
colnames(mat.r) <- c("node1","node2","link")
colnames(mat.p) <- c("node1","node2","p")
mat.r <- mat.r[which(mat.p[,3] < 0.01),]  ## 
mat.r <- mat.r[which(mat.r[,3] > 0.5),]  ##  
mat.r <- as.data.frame(mat.r)
label <- read.csv("cluster_roi.csv",stringsAsFactors = F,row.names = 1)
data.hh.roi <- data.all[rownames(label)[which(label$cluster_new == "highest-heterogeneity")],]
data.mh.roi <- data.all[rownames(label)[which(label$cluster_new == "middle-heterogeneity")],]
data.lh.roi <- data.all[rownames(label)[which(label$cluster_new == "low-heterogeneity")],]
gene.report <- data.frame(name = c(rownames(data.hh),
                                   rownames(data.mh),
                                   rownames(data.lh)),
                          type = c(rep("hh",dim(data.hh)[1]),
                                   rep("mh",dim(data.mh)[1]),
                                   rep("lh",dim(data.lh)[1])),
                          value = all_up_fc$log2FoldChange)
img <- graph_from_data_frame(mat.r, directed = F)
deg<-degree(img,mode="all")
gene.report <- gene.report[match(names(deg),gene.report$name),]
img <- graph_from_data_frame(mat.r, directed = F, vertices = gene.report)
deg<-degree(img,mode="all")
deg.fc <- all_up_fc[names(deg),2]
names(deg.fc) <- names(deg)
ggraph(img,layout = 'kk') +   
  geom_edge_link(color="grey",size = 0.1)+ 
  geom_node_point(aes(color = type,
                      size = deg.fc)) +   ## 
  scale_color_manual(values = c("#e59572","#d3bd74","#9baec8"))+
  geom_node_text(aes(label = name)) +
  theme_graph()
## h-ITH
gene <- hh[names(deg),2]
gene.report.hh <- data.frame(name = names(deg),
                             type = c(rep("hh",36),
                                      rep("mh",48),
                                      rep("lh",26)),
                             value = all_up_fc[names(deg),2],
                             gene = gene)
img.hh <- graph_from_data_frame(mat.r, directed = F, vertices = gene.report.hh)
deg.hh <- degree(img.hh,mode="all")
ggraph(img.hh,layout = 'kk') +   
  geom_edge_link(color="grey",size = 0.1)+ 
  geom_node_point(aes(color = gene,
                      size = deg.fc)) +   
  #scale_color_manual(values = c("#e59572","#d3bd74","#9baec8"))+
  theme_graph() + 
  scale_colour_gradient2(low = "#afb4ff", high = "#fc3a52", midpoint = 0)
## m-ITH
gene <- mh[names(deg),2]
gene.report.mh <- data.frame(name = names(deg),
                             type = c(rep("hh",36),
                                      rep("mh",48),
                                      rep("lh",26)),
                             value = all_up_fc[names(deg),2],
                             gene = gene)
img.mh <- graph_from_data_frame(mat.r, directed = F, vertices = gene.report.mh)
deg.mh <- degree(img.mh,mode="all")
ggraph(img.mh,layout = 'kk') +   
  geom_edge_link(color="grey",size = 0.1)+ 
  geom_node_point(aes(color = gene,
                      size = deg.fc)) +  
  #scale_color_manual(values = c("#e59572","#d3bd74","#9baec8"))+
  theme_graph() + 
  scale_colour_gradient2(low = "#afb4ff", high = "#fc3a52", midpoint = 0)
## l-ITH
gene <- lh[names(deg),2]
gene.report.lh <- data.frame(name = names(deg),
                             type = c(rep("hh",36),
                                      rep("mh",48),
                                      rep("lh",26)),
                             value = all_up_fc[names(deg),2],
                             gene = gene)
img.lh <- graph_from_data_frame(mat.r, directed = F, vertices = gene.report.lh)
deg.lh <- degree(img.lh,mode="all")
ggraph(img.lh,layout = 'kk') +   
  geom_edge_link(color="grey",size = 0.1)+ 
  geom_node_point(aes(color = gene,
                      size = deg.fc)) +   
  #scale_color_manual(values = c("#e59572","#d3bd74","#9baec8"))+
  theme_graph() + 
  scale_colour_gradient2(low = "#afb4ff", high = "#fc3a52", midpoint = 0)


###  gene network enrichment  ###
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
hh1 <- c('ATP5MG','ATRX','KCNJ9','SUPT7L','CCDC13','TMEM30A','OBSCN','CYP2W1',
         'ACTA2','POSTN','TAGLN','SPARC','COL4A1','COL1A1','COL1A2','IGFBP4',
         'COL6A2','VIM','IGFBP7','TIMP1','IGHA1','JCHAIN','IGKC','IGHG3',
         'IGHG4','IGHG2','IGHG1')
hh1 <- bitr(hh1,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
hh1 <- dplyr::distinct(hh1,SYMBOL,.keep_all=TRUE) 

mh1 <- c('FGF3','ASCL2','ZIC1','ISM1','MYH7B','AXIN2','FGF8','MAPK13',
         'MME','KRT5','NECTIN1','EIF4EBP1','ATP1A1','DKK4','DEFA6','DKK1',
         'IGFBPL1','NKD1','DEFA5','FGF19','BSX','CST1')
mh1 <- bitr(mh1,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")  
mh1 <- dplyr::distinct(mh1,SYMBOL,.keep_all=TRUE)

lh1 <- c('IL32','STAT1','NREP','CDKAL1','RPL19','RPS28','RPS21','RPL37A','LAPTM4B',
         'TAP1','DLL1','RPL17','RPS5','RBM38','CD274','RPS12','NCOA7','GBP5',
         'GSTP1','BCL2','SRSF6','TGIF2','RIC1','COLCA2','UHRF2','SDC4','FGL1',
         'KRT7','FGF14','PLGRKT','KCNMB4','SMARCA2','ERMP1','ANXA1','ACADSB',
         'IRX2','FOXC1','CHSY3','NEUROG2','KAZN','MYC','SOX9','AZGP1','EYA1',
         'L3MBTL4')
lh1 <- bitr(lh1,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
lh1 <- dplyr::distinct(lh1,SYMBOL,.keep_all=TRUE) 

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

hh1.res <- enricher(hh1$ENTREZID, TERM2GENE=m_t2g,
                    pvalueCutoff = 1,minGSSize = 1,
                    maxGSSize = 500,qvalueCutoff = 1)
hh1.res <- as.data.frame(hh1.res)
mh1.res <- enricher(mh1$ENTREZID, TERM2GENE=m_t2g,
                    pvalueCutoff = 1,minGSSize = 1,
                    maxGSSize = 500,qvalueCutoff = 1)
mh1.res <- as.data.frame(mh1.res)
lh1.res <- enricher(lh1$ENTREZID, TERM2GENE=m_t2g,
                    pvalueCutoff = 0.05,minGSSize = 1,
                    maxGSSize = 500,qvalueCutoff = 1)
lh1.res <- as.data.frame(lh1.res)


### clusterprofiler functional enrichment
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(dplyr)
# h-ITH
gene <- read.csv("m-ITH-gene.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
data <- data.frame(SYMBOL = rownames(gene),
                   logFC = gene$log2FoldChange)
gene <- data$SYMBOL
gene <- bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")  #
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE) #
data_all <- data %>% 
  inner_join(gene,by="SYMBOL")
data_all_sort <- data_all %>% 
  arrange(desc(logFC))
geneList = data_all_sort$logFC #
names(geneList) <- data_all_sort$ENTREZID #
kegg_gmt <- read.gmt("c2.cp.kegg.v2022.1.Hs.entrez.gmt") #
gsea1 <- GSEA(geneList,
              minGSSize = 1,
              pvalueCutoff = 0.05,
              pAdjustMethod = "fdr",
              TERM2GENE = kegg_gmt) #GSEA
gsea <- as.data.frame(gsea1)
gsea <- gsea[which(gsea$p.adjust < 0.05),]
# m-ITH
gene <- read.csv("m-ITH-gene.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
data <- data.frame(SYMBOL = rownames(gene),
                   logFC = gene$log2FoldChange)
gene <- data$SYMBOL
gene <- bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")  #
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE) #
data_all <- data %>% 
  inner_join(gene,by="SYMBOL")
data_all_sort <- data_all %>% 
  arrange(desc(logFC))
geneList = data_all_sort$logFC #
names(geneList) <- data_all_sort$ENTREZID #
kegg_gmt <- read.gmt("c2.cp.kegg.v2022.1.Hs.entrez.gmt") #
gsea1 <- GSEA(geneList,
              minGSSize = 1,
              pvalueCutoff = 0.05,
              pAdjustMethod = "fdr",
              TERM2GENE = kegg_gmt) #GSEA
gsea <- as.data.frame(gsea1)
gsea <- gsea[which(gsea$p.adjust < 0.05),]
# l-ITH
gene <- read.csv("l-ITH-gene.csv",stringsAsFactors = F,
                 row.names = 1,check.names = F)
data <- data.frame(SYMBOL = rownames(gene),
                   logFC = gene$log2FoldChange)
gene <- data$SYMBOL
gene <- bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")  #
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE) #
data_all <- data %>% 
  inner_join(gene,by="SYMBOL")
data_all_sort <- data_all %>% 
  arrange(desc(logFC))
geneList = data_all_sort$logFC #
names(geneList) <- data_all_sort$ENTREZID #
kegg_gmt <- read.gmt("c2.cp.kegg.v2022.1.Hs.entrez.gmt") #
gsea1 <- GSEA(geneList,
              minGSSize = 1,
              pvalueCutoff = 0.05,
              pAdjustMethod = "fdr",
              TERM2GENE = kegg_gmt) #GSEA
gsea <- as.data.frame(gsea1)
gsea <- gsea[which(gsea$p.adjust < 0.05),]



###  tumor immune microenvironment  ###
library(pheatmap)
library(ggplot2)
library(ggpubr)
box.fun <- function(x,y){
  ggplot(data = x, aes(x = cluster , y = x[,y],fill=cluster)) + 
    stat_compare_means(label.y = median(x[,y]), label.x = 1.5)  +  ###
    stat_compare_means(comparisons = my_com,label = "p.signif")  +
    geom_jitter(aes(fill=cluster),width =0.2,shape = 21,size=2.5)+
    geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,
                 width = 0.8,alpha = 0.4)+ #
    scale_fill_manual(values = c("#e59572","#9BAEC8","#D3BD74"))+ #
    theme_bw()+ 
    theme(legend.text=element_text(colour="black",size=10),
          legend.title=element_text(colour="black",size=10),
          panel.grid.major = element_blank(),legend.position="none")+  
    ylab(colnames(x)[y])
}
timer <- read.csv("TIMER_ROI.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
mcp <- read.csv("MCPcounter_ROI.csv",stringsAsFactors = F,
                row.names = 1,check.names = F)
ciber <- read.csv("Cibersort_ROI.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
cluster <- read.csv("cluster_roi.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
cluster <- cluster[order(cluster$cluster_new),]
timer <- timer[,rownames(cluster)]
mcp <- mcp[,rownames(cluster)]
ciber <- ciber[,rownames(cluster)]
barplot(as.matrix(ciber),col = rep(c('#8bbfdf','#0063a4','#005a6b','#e80b7d','#f9b66e',
                                     '#e11a2a','#ee6b92','#00417c','#518f36','#f2a41b',
                                     '#9d6a13','#98c788','#f2938c','#149e4b','#fcc214',
                                     '#fcc214','#df1a30','#bda7c5','#ec4334','#7b2874',
                                     '#f28237','#007dba'),22),
        legend.text = c(rownames(ciber)))
# cibrsort boxplot
ciber_new <- as.data.frame(t(ciber)[rownames(cluster),])
ciber_new$cluster <- cluster$cluster_new
ciber_new$cluster <- factor(ciber_new$cluster,
                            levels = c("highest-heterogeneity",
                                       "middle-heterogeneity",
                                       "low-heterogeneity"))
my_com <- list(c("highest-heterogeneity","middle-heterogeneity"),
               c("highest-heterogeneity","low-heterogeneity"),
               c("middle-heterogeneity","low-heterogeneity"))
for (i in 1:22){
  a <- kruskal.test(ciber_new[,i]~ciber_new$cluster)
  if (a$p.value < 0.05){
    print(i)
  }
}
a1 <- box.fun(ciber_new,3)
# cibrsort boxplot
ciber_new <- as.data.frame(t(ciber)[rownames(cluster),])
ciber_new$cluster <- cluster$cluster_new
ciber_new$cluster <- factor(ciber_new$cluster,
                            levels = c("highest-heterogeneity",
                                       "middle-heterogeneity",
                                       "low-heterogeneity"))
my_com <- list(c("highest-heterogeneity","middle-heterogeneity"),
               c("highest-heterogeneity","low-heterogeneity"),
               c("middle-heterogeneity","low-heterogeneity"))
for (i in 1:22){
  a <- kruskal.test(ciber_new[,i]~ciber_new$cluster)
  if (a$p.value < 0.05){
    print(i)
  }
}
a1 <- box.fun(ciber_new,3)
a2 <- box.fun(ciber_new,4)
a3 <- box.fun(ciber_new,5)
a4 <- box.fun(ciber_new,17)
ggarrange(a1,a2,a3,a4,nrow = 2,ncol = 2)

mat <- data.frame(fraction = as.numeric(as.matrix(ciber_new[,c(-3,-4,-5,-17,-23)])),
                  cluster = rep(ciber_new$cluster,18),
                  cell = c(rep(colnames(ciber_new)[1],79),rep(colnames(ciber_new)[2],79),
                           rep(colnames(ciber_new)[6],79),rep(colnames(ciber_new)[7],79),
                           rep(colnames(ciber_new)[8],79),rep(colnames(ciber_new)[9],79),
                           rep(colnames(ciber_new)[10],79),rep(colnames(ciber_new)[11],79),
                           rep(colnames(ciber_new)[12],79),rep(colnames(ciber_new)[13],79),
                           rep(colnames(ciber_new)[14],79),rep(colnames(ciber_new)[15],79),
                           rep(colnames(ciber_new)[16],79),rep(colnames(ciber_new)[18],79),
                           rep(colnames(ciber_new)[19],79),rep(colnames(ciber_new)[20],79),
                           rep(colnames(ciber_new)[21],79),rep(colnames(ciber_new)[22],79)))
ggplot(data = mat, 
       aes(x = cell , y = fraction,fill=factor(cluster))) + 
  #geom_jitter(aes(fill=cluster),width =0.2,shape = 21,size=2.5)+
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,
               width = 0.8,alpha = 0.4)+ #
  scale_fill_manual(values = c("#e59572","#D3BD74","#9BAEC8"))+ #
  theme_bw()+ #
  theme(legend.text=element_text(colour="black",size=10),
        legend.title=element_text(colour="black",size=10),
        panel.grid.major = element_blank())+  
  ylab("Fraction")
# timer
timer_new <- as.data.frame(t(timer)[rownames(cluster),])
timer_new$cluster <- cluster$cluster_new
timer_new$cluster <- factor(timer_new$cluster,
                            levels = c("highest-heterogeneity",
                                       "middle-heterogeneity",
                                       "low-heterogeneity"))
my_com <- list(c("highest-heterogeneity","middle-heterogeneity"),
               c("highest-heterogeneity","low-heterogeneity"),
               c("middle-heterogeneity","low-heterogeneity"))
for (i in 1:22){
  a <- kruskal.test(timer_new[,i]~timer_new$cluster)
  if (a$p.value < 0.05){
    print(i)
  }
}
b1 <- box.fun(timer_new,1)
b2 <- box.fun(timer_new,2)
b3 <- box.fun(timer_new,3)
b4 <- box.fun(timer_new,6)
# MCPCOUNTER
mcp_new <- as.data.frame(t(mcp)[rownames(cluster),])
mcp_new$cluster <- cluster$cluster_new
mcp_new$cluster <- factor(mcp_new$cluster,
                          levels = c("highest-heterogeneity",
                                     "middle-heterogeneity",
                                     "low-heterogeneity"))
my_com <- list(c("highest-heterogeneity","middle-heterogeneity"),
               c("highest-heterogeneity","low-heterogeneity"),
               c("middle-heterogeneity","low-heterogeneity"))
for (i in 1:11){
  a <- kruskal.test(mcp_new[,i]~mcp_new$cluster)
  if (a$p.value < 0.05){
    print(i)
  }
}
c1 <- box.fun(mcp_new,1)
c2 <- box.fun(mcp_new,2)
c3 <- box.fun(mcp_new,3)
c4 <- box.fun(mcp_new,5)
c5 <- box.fun(mcp_new,8)
ggarrange(b1,b2,b3,c1,c2,c3,c4,nrow = 1,ncol = 7)


###  ROI  IHC  validation  ###
cluster <- read.csv("cluster_roi.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
cluster <- cluster[order(cluster$cluster_new),]
cd8 <- read.csv("CD8_Hscore.csv",
                stringsAsFactors = F,row.names = 1,check.names = F)
cd8 <- cd8[rownames(cluster),]
cd8_new <- as.data.frame(cd8[rownames(cluster),25])
cd8_new$cluster <- cluster$cluster_new
cd8_new$cluster <- factor(cd8_new$cluster,
                          levels = c("highest-heterogeneity",
                                     "middle-heterogeneity",
                                     "low-heterogeneity"))
my_com <- list(c("highest-heterogeneity","middle-heterogeneity"),
               c("highest-heterogeneity","low-heterogeneity"),
               c("middle-heterogeneity","low-heterogeneity"))
ebtop<-function(x){
  return(mean(x)+sd(x)/sqrt(length(x)))
}
ebbottom<-function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
ggplot(data=cd8_new,aes(x=cluster,y=cd8_new[,1],fill=cluster))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",fun.min = ebbottom,
               fun.max = ebtop,position = position_dodge(0.9),width=0.2)+
  stat_compare_means(label.x = 1.5)  +  ### 
  stat_compare_means(comparisons = my_com,label = "p.signif")  + 
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  geom_jitter(aes(fill=cluster),width =0.3,shape = 21,size=1.5) +
  theme_bw()+ theme(panel.grid = element_blank())+ ylab("CD8 H-score")+
  scale_fill_manual(values = c("#e59572","#D3BD74","#9BAEC8"),name="")



###  relation with ROI on patient  ###
library(plyr)
library(ggplot2)
label <- read.csv("cluster_roi.csv",stringsAsFactors = F,row.names = 1)
label$number <- 1
label1 <- ddply(label,'cluster_patient',transform,percent = 1/sum(number)*100)
ggplot(label1)+
  scale_fill_manual(values = c("#e59572","#D3BD74","#9BAEC8"))+ 
  geom_bar(aes(x=patient,fill=cluster_new),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 90), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())



###  ITH score boxplot on patient ###
library(ggplot2)
library(ggpubr)
ith <- read.csv("ITH_score.csv",stringsAsFactors = F,row.names = 1)
cluster <- read.csv("cluster_roi.csv",stringsAsFactors = F,row.names = 1)
ith <- ith[rownames(cluster),]
ith$cluster <- cluster$patient_label
ggplot(data = ith, aes(x = cluster , y = ITH.score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1.5)  +  ### 
  geom_violin(trim=F,color="white",width = 0.8) +
  geom_jitter(aes(fill=cluster),width =0.2,shape = 21,size=3.5)+
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.8,fill = "#FFFFFF",alpha = 0.4)+ #
  scale_fill_manual(values = c("#f06966","#a39391"))+ 
  theme_bw()+ 
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=10,angle = 0), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("ITH score")


###  patient KM-curv  ###
library(survival)
library(survminer)
cluster <- read.csv("cluster_patient.csv",stringsAsFactors = F,row.names = 1)
clinical <- read.csv("clinical.csv",
                     stringsAsFactors = F,row.names = 1)
clinical <- clinical[rownames(cluster),]
clinical$cluster <- cluster$cluster
# OS
surv <- survfit(Surv(as.numeric(OS),OSState)~cluster,data = clinical)
surv
survdiff(Surv(as.numeric(OS),OSState)~cluster,data = clinical)
summary(surv)
summary(coxph(Surv(as.numeric(OS),OSState)~cluster,data = clinical))
ggsurvplot(surv,
           pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", 
           legend = c(2,0.5), #
           legend.title = "", 
           conf.int = F,
           xlab = "Time (months)",
           ylab = "Overall Survival"
)
# DFS
surv <- survfit(Surv(as.numeric(DFS),DFSState)~cluster,data = clinical)
surv
survdiff(Surv(as.numeric(DFS),DFSState)~cluster,data = clinical)
summary(surv)
summary(coxph(Surv(as.numeric(DFS),DFSState)~cluster,data = clinical))
ggsurvplot(surv,
           pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           legend = c(2,0.5), # 
           legend.title = "", # 
           conf.int = F,
           xlab = "Time (months)",
           ylab = "Disease-free Survival"
)

###  patient survival event  ###
library(plyr)
library(ggplot2)
cluster <- read.csv("cluster_patient.csv",stringsAsFactors = F,row.names = 1)
clinical <- read.csv("clinical.csv",
                     stringsAsFactors = F,row.names = 1)
clinical <- clinical[rownames(cluster),]
clinical$cluster <- cluster$cluster_new1
clinical$OSState <- as.factor(clinical$OSState)
clinical$DFSState <- as.factor(clinical$DFSState)
# DFS
clinical$number <- 1
clinical <- ddply(clinical,'DFSState',transform,percent = 1/sum(number)*100)
ggplot(clinical)+
  scale_fill_manual(values = c("#9DD3A8","#F0B775"))+ 
  geom_bar(aes(x=cluster,fill=DFSState),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL)) 



###  relation between ITH subtype and traditional subtype ###
library(tidyverse)
library(viridis)
library(patchwork)
library(networkD3)
MisLinks <- read.csv("link_matrix.csv",stringsAsFactors = F)
MisNodes <- read.csv("node_matrix.csv",stringsAsFactors = F)
Node2index = list()
Node2index[MisNodes$name] = 0:length(MisNodes$name)
MisLinks = MisLinks %>%
  mutate(source2 = unlist(Node2index[source])) %>%
  mutate(target2 = unlist(Node2index[target]))
# color
color2project = paste(unique(MisNodes$group_color),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',color2project,'"]).range(["',color2project,'"])')
sankeyNetwork(Links = MisLinks, 
              Nodes = MisNodes,
              Source = "source2", 
              Target = "target2",
              Value ="value",
              NodeID = "name",
              NodeGroup = "group_color", 
              colourScale = JS(my_color),
              fontSize = 0
)



###  immune molecular among subypes  ###
library(ggplot2)
library(ggpubr)
data <- read.csv("WTA_multi_ROIdata.csv",stringsAsFactors = F,row.names = 1)
cluster <- read.csv("cluster_roi.csv",stringsAsFactors = F,
                    row.names = 1,check.names = F)
cluster <- cluster[order(cluster$cluster_new),]
data <- t(data)[rownames(cluster),]
mhc1 <- c("HLA-A","HLA-B","HLA-C","TAP1","TAP2","B2M")
mhc2 <- c("HLA-DPA1","HLA-DPB1","HLA-DPB2","HLA-DQA1","HLA-DQA2",
          "HLA-DQB1","HLA-DQB2","HLA-DRB1","HLA-DRB5","HLA-DRB6")
mhc.other <- c("HLA-E","HLA-F","HLA-G","HLA-H")
coinhi <- c('CTLA4','TIGIT','BTLA','CD48','PDCD1','LAG3',
            'CD274','HACR2','BTN2A2','LAIR1','BTN3A1',
            'PDCD1LG2','BTN1A1','VTCN1','BTNL2')
costi <- c('ICOS','TNFRSF9','CD70','CD80','TNFRSF13C','TMIGD2',
           'TNFRSF13B','CD27','SLAMF1','TNFSF13B','CD86','TNFSF4',
           'TNFRSF18','CD28','CD226','TNFSF9','TNFSF8','TNFSF18',
           'HAVCR1','TNFRSF4','TNFSF15','TNFSF13','CD58','ICOSLG',
           'TNFRSF8','TNFRSF14','BTNLB')
data.new <- data[,intersect(c(mhc1,mhc2,mhc.other,coinhi,costi),colnames(data))]
data.new.h <- data.new[rownames(cluster)[which(cluster$cluster_new == "highest-heterogeneity")],]
data.new.m <- data.new[rownames(cluster)[which(cluster$cluster_new == "middle-heterogeneity")],]
data.new.l <- data.new[rownames(cluster)[which(cluster$cluster_new == "low-heterogeneity")],]
fc.h.log2 <- log2(colMeans(data.new.h)/colMeans(rbind(data.new.m,data.new.l)))
fc.m.log2 <- log2(colMeans(data.new.m)/colMeans(rbind(data.new.h,data.new.l)))
fc.l.log2 <- log2(colMeans(data.new.l)/colMeans(rbind(data.new.m,data.new.h)))
p <- c()
for (i in 1:56){
  a <- kruskal.test(c(data.new.h[,i],data.new.m[,i],data.new.l[,i])~cluster$cluster_new)
  p <- c(p,a$p.value)
}
names(p) <- colnames(data.new)
heat.data <- rbind(fc.h.log2,fc.m.log2,fc.l.log2)

red <- "#D94E48";
blue <- "#5175A4";
white <- rgb(255,255,255,maxColorValue = 255)
out <- pheatmap(heat.data,fontsize=6,cutree_cols = 3,
                color  = colorRampPalette(c(blue,white,red))(100),
                clustering_method = "ward.D2",
                border_color = "grey60",
                cluster_cols = F, cluster_rows = F,
                show_rownames = T, show_colnames = T)


















