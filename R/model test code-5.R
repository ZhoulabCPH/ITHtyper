######    Extrenal and inhouse data prediction   #######   
gene <- c("NKX1-2","TLE2","TPBG","SRSF6","DAZ4", "GPR31","CD274" ,"LYZ", "PCP4","ZIC1" )
###  WTA 
roi <- read.csv("cluster_roi.csv",stringsAsFactors = F,
                row.names = 1,check.names = F)
train_name <- roi_duiying[names(XG_pred_train1),]
vali_name <-roi_duiying[names(XG_pred_vali1),]
test_name <- roi_duiying[names(XG_pred_test1),]
data_wta <- read.csv("WTA_log2_rawdata.csv",
                     stringsAsFactors = F,row.names = 1,check.names = F)
data_wta <- data_wta[gene,]
clinical_wta <- read.csv("clinical.csv",stringsAsFactors = F,
                         row.names = 1,check.names = F)
name <- intersect(colnames(data_wta),rownames(clinical_wta))
data_wta <- data_wta[,name]
clinical_wta <- clinical_wta[name,]
clinical_train <- clinical_wta[train_name,]
clinical_train$predict_label <- XG_pred_train1
clinical_vali <- clinical_wta[vali_name,]
clinical_vali$predict_label <- XG_pred_vali1
clinical_test <- clinical_wta[test_name,]
clinical_test$predict_label <- XG_pred_test1

### george cohort
data_george <- read.table("george_81tumor_18816gene_data_fpkm.txt",
                          stringsAsFactors = F,row.names = 1,check.names = F)
data_george <- data_george[gene,]
clinical_george <- read.csv("clinical.csv",stringsAsFactors = F,
                            row.names = 1,check.names = F)
name <- intersect(colnames(data_george),rownames(clinical_george))
data_george <- data_george[,name]
clinical_george <- clinical_george[name,]
### GSE60052
data_GSE60052 <- read.table("GSE60052_78tumor_51131gene_count.txt",
                            stringsAsFactors = F,row.names = 1,check.names = F)
data_GSE60052 <- data_GSE60052[gene,]
clinical_GSE60052 <- read.csv("clinical.csv",stringsAsFactors = F,
                              row.names = 1,check.names = F)
data_GSE60052 <- data_GSE60052[,rownames(clinical_GSE60052)]
### Roper
data_Roper <- read.csv("Roper_data_rpkm.csv",stringsAsFactors = F,
                       row.names = 1,check.names = F)
data_Roper <- data_Roper[gene,]
rownames(data_Roper)[5] <- gene[5]
data_Roper[5,] <- rep(0,dim(data_Roper)[2])
clinical_Roper <- read.csv("clinical.csv",stringsAsFactors = F,
                           row.names = 1,check.names = F)
data_Roper <- data_Roper[,rownames(clinical_Roper)]
###  input data type construct 
data_wta_input <- apply(data_wta, 1, scale)
rownames(data_wta_input) <- colnames(data_wta)
data_Roper_input <- apply(data_Roper, 1, scale)
rownames(data_Roper_input) <- colnames(data_Roper)
data_george_input <- t(data_george)
data_GSE60052_input <- t(data_GSE60052)
XG_pred_wta_score <- predict(mxgb4m,data_wta_input,type = "response")
XG_pred_george_score <- predict(mxgb4m,data_george_input,type = "response")
names(XG_pred_george_score) <- rownames(data_george_input)
XG_pred_GSE60052_score <- predict(mxgb4m,data_GSE60052_input,type = "response")
names(XG_pred_GSE60052_score) <- rownames(data_GSE60052_input)
XG_pred_Roper_score <- predict(mxgb4m,data_Roper_input,type = "response")
names(XG_pred_Roper_score) <- rownames(data_Roper_input)
###  label prediction 
XG_pred_wta_label <- XG_pred_wta_score
XG_pred_wta_label[XG_pred_wta_score>=0.45] <- 1  ### ITHtyper high
XG_pred_wta_label[XG_pred_wta_score<0.45] <- 0  ### ITHtyper low
XG_pred_george_label <- XG_pred_george_score
XG_pred_george_label[XG_pred_george_score>=0.45] <- 1  ### ITHtyper high
XG_pred_george_label[XG_pred_george_score<0.45] <- 0  ### ITHtyper low
XG_pred_GSE60052_label <- XG_pred_GSE60052_score
XG_pred_GSE60052_label[XG_pred_GSE60052_score>=0.45] <- 1  ###  ITHtyper high
XG_pred_GSE60052_label[XG_pred_GSE60052_score<0.45] <- 0  ### ITHtyper low
XG_pred_Roper_label <- XG_pred_Roper_score
XG_pred_Roper_label[XG_pred_Roper_score>=0.45] <- 1  ###  ITHtyper high
XG_pred_Roper_label[XG_pred_Roper_score<0.45] <- 0  ### ITHtyper low
###  merge clinical information and ITHtyper label
clinical_wta$score <- XG_pred_wta_score
XG_pred_george_score <- XG_pred_george_score[name]
clinical_george <- clinical_george[name,]
clinical_george$score <- XG_pred_george_score
XG_pred_GSE60052_score <- XG_pred_GSE60052_score[rownames(clinical_GSE60052)]
clinical_GSE60052$score <- XG_pred_GSE60052_score
XG_pred_Roper_score <- XG_pred_Roper_score[rownames(clinical_Roper)]
clinical_Roper$score <- XG_pred_Roper_score
clinical_wta$hlabel <- as.factor(XG_pred_wta_label)
XG_pred_george_score <- XG_pred_george_score[name]
clinical_george <- clinical_george[name,]
clinical_george$hlabel <- as.factor(XG_pred_george_label)
XG_pred_GSE60052_label <- XG_pred_GSE60052_label[rownames(clinical_GSE60052)]
clinical_GSE60052$hlabel <- as.factor(XG_pred_GSE60052_label)
XG_pred_Roper_label <- XG_pred_Roper_label[rownames(clinical_Roper)]
clinical_Roper$hlabel <- as.factor(XG_pred_Roper_label)

###  ITHtyper in survival analysis  ###
library(survival)
library(survminer)
library(plyr)
library(ggplot2)
###  training
surv <- survfit(Surv(as.numeric(OS),OSState)~predict_label,data = clinical_train)
surv
survdiff(Surv(as.numeric(OS),OSState)~predict_label,data = clinical_train)
summary(surv)
summary(coxph(Surv(as.numeric(OS),OSState)~predict_label,data = clinical_train))
ggsurvplot(surv,
           pval = TRUE,#
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette = c("#1c1e24","#A39391","#F06966"),
           legend = c(2,0.5), #
           legend.title = "", # 
           conf.int = F,#
           add.all = TRUE,
           legend.labs = c("All","ITHtyper low","ITHtyper high"), # 指定图例分组标签
           xlab = "Time (months)",
           ylab = "Overall Survival"
)
clinical_train$predict_label <- as.factor(clinical_train$predict_label)
clinical_train$OSState <- as.factor(clinical_train$OSState)
clinical_train$number <- 1
clinical_train <- ddply(clinical_train,'OSState',transform,percent = 1/sum(number)*100)
ggplot(clinical_train)+
  scale_fill_manual(values = c("#475F77","#EBB481"))+ #
  geom_bar(aes(x=predict_label,fill=OSState),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))
surv <- survfit(Surv(as.numeric(DFS),DFSState)~predict_label,data = clinical_train)
surv
survdiff(Surv(as.numeric(DFS),DFSState)~predict_label,data = clinical_train)
summary(surv)
summary(coxph(Surv(as.numeric(DFS),DFSState)~predict_label,data = clinical_train))
ggsurvplot(surv,
           pval = TRUE,#
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette = c("#1c1e24","#A39391","#F06966"),
           legend = c(2,0.5), # 
           legend.title = "", # 
           conf.int = F,#
           add.all = TRUE,
           legend.labs = c("All","ITHtyper low","ITHtyper high"), #
           xlab = "Time (months)",
           ylab = "Disease−free Survival")
clinical_train$predict_label <- as.factor(clinical_train$predict_label)
clinical_train$DFSState <- as.factor(clinical_train$DFSState)
clinical_train$number <- 1
clinical_train <- ddply(clinical_train,'DFSState',transform,percent = 1/sum(number)*100)
ggplot(clinical_train)+
  scale_fill_manual(values = c("#475F77","#EBB481"))+ #
  geom_bar(aes(x=predict_label,fill=DFSState),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))


###  testing
clinical_vali_test <- rbind(clinical_vali,clinical_test)
clinical_vali_test <- clinical_vali_test[nchar(rownames(clinical_vali_test))<7,]
surv <- survfit(Surv(as.numeric(OS),OSState)~predict_label,data = clinical_vali_test)
surv
survdiff(Surv(as.numeric(OS),OSState)~predict_label,data = clinical_vali_test)
summary(surv)
summary(coxph(Surv(as.numeric(OS),OSState)~predict_label,data = clinical_vali_test))
ggsurvplot(surv,
           pval = TRUE,#
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette = c("#1c1e24","#A39391","#F06966"),
           legend = c(2,0.5), # 
           legend.title = "", # 
           conf.int = F,#
           add.all = TRUE,
           legend.labs = c("All","ITHtyper low","ITHtyper high"), # 
           xlab = "Time (months)",
           ylab = "Overall Survival"
)
clinical_vali_test$predict_label <- as.factor(clinical_vali_test$predict_label)
clinical_vali_test$OSState <- as.factor(clinical_vali_test$OSState)
clinical_vali_test$number <- 1
clinical_vali_test <- ddply(clinical_vali_test,'OSState',transform,percent = 1/sum(number)*100)
ggplot(clinical_vali_test)+
  scale_fill_manual(values = c("#475F77","#EBB481"))+ 
  geom_bar(aes(x=predict_label,fill=OSState),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))

surv <- survfit(Surv(as.numeric(DFS),DFSState)~predict_label,data = clinical_vali_test)
surv
survdiff(Surv(as.numeric(DFS),DFSState)~predict_label,data = clinical_vali_test)
summary(surv)
summary(coxph(Surv(as.numeric(DFS),DFSState)~predict_label,data = clinical_vali_test))
ggsurvplot(surv,
           pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette = c("#1c1e24","#A39391","#F06966"),
           legend = c(2,0.5), 
           legend.title = "", 
           conf.int = F,
           add.all = TRUE,
           legend.labs = c("All","ITHtyper low","ITHtyper high"), 
           xlab = "Time (months)",
           ylab = "Disease−free Survival"
)
clinical_vali_test <- clinical_vali_test[which(is.na(clinical_vali_test$OSState) == F),]
clinical_vali_test$predict_label <- as.factor(clinical_vali_test$predict_label)
clinical_vali_test$DFSState <- as.factor(clinical_vali_test$DFSState)
clinical_vali_test$number <- 1
clinical_vali_test <- ddply(clinical_vali_test,'OSState',transform,percent = 1/sum(number)*100)
ggplot(clinical_vali_test)+
  scale_fill_manual(values = c("#475F77","#EBB481"))+ 
  geom_bar(aes(x=predict_label,fill=DFSState),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))
## Roper  ICB response
roc1 <- roc(clinical_Roper$ICB_response,clinical_Roper$score,
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "#2694ab")
roc2 <- roc(clinical_Roper$ICB_response,as.numeric(clinical_Roper$hlabel),
            plot=TRUE, print.thres=TRUE, ci=TRUE,
            print.auc=TRUE,legacy.axes = TRUE,col = "#2694ab")
roc(clinical_Roper$ICB_response,clinical_Roper$score)
roc(clinical_Roper$ICB_response,as.numeric(clinical_Roper$hlabel))

data_Roper <- read.csv("data_rpkm.csv",stringsAsFactors = F,row.names = 1,check.names = F)
data_Roper <- as.data.frame(t(data_Roper)[rownames(clinical_Roper),])
pd1 <- data_Roper$PDCD1
pd1_label <- pd1
pd1_label[pd1>median(pd1)] <- 1  ### ITHtyper high
pd1_label[pd1<median(pd1)] <- 0  ### ITHtyper low
pdl1 <- data_Roper$CD274
pdl1_label <- pdl1
pdl1_label[pdl1>median(pdl1)] <- 1  ### ITHtyper high
pdl1_label[pdl1<median(pdl1)] <- 0  ### ITHtyper low

roc_pd1 <- roc(clinical_Roper$ICB_response,pd1)
roc_pd1_label <- roc(clinical_Roper$ICB_response,as.numeric(pd1_label))
roc(clinical_Roper$ICB_response,pd1, plot=TRUE, print.thres=TRUE, 
    ci=TRUE,print.auc=TRUE,legacy.axes = TRUE,col = "#E73A38")
roc(clinical_Roper$ICB_response,as.numeric(pd1_label), plot=TRUE,
    print.thres=TRUE, ci=TRUE,print.auc=TRUE,legacy.axes = TRUE,col = "#E73A38")
roc_pdl1 <- roc(clinical_Roper$ICB_response,pdl1)
roc_pdl1_label <- roc(clinical_Roper$ICB_response,as.numeric(pdl1_label))
roc(clinical_Roper$ICB_response,pdl1, plot=TRUE, print.thres=TRUE, 
    ci=TRUE,print.auc=TRUE,legacy.axes = TRUE,col = "#E73A38")
roc(clinical_Roper$ICB_response,as.numeric(pdl1_label), plot=TRUE,
    print.thres=TRUE, ci=TRUE,print.auc=TRUE,legacy.axes = TRUE,col = "#E73A38")

clinical_Roper <- clinical_Roper[which(is.na(clinical_Roper$ICB_response) == F),]
clinical_Roper$ICB_response <- as.factor(clinical_Roper$ICB_response)
wilcox.test(clinical_Roper$score~as.factor(clinical_Roper$ICB_response))
ggplot(clinical_Roper,aes(ICB_response,score))+
  stat_boxplot(geom = "errorbar",width=0.15)+
  geom_boxplot(aes(fill=ICB_response)) +
  scale_fill_manual(values = c("#7F95D1","#D87575")) +
  geom_jitter(aes(fill=ICB_response),width =0.2,shape = 21,size=3.5) +
  stat_compare_means(cluster = 0.6, label.x = 1.5)
table(paste(clinical_Roper$ICB_response,clinical_Roper$hlabel,sep="_"))

clinical_Roper$hlabel <- as.factor(clinical_Roper$hlabel)
clinical_Roper$ICB_response <- as.factor(clinical_Roper$ICB_response)
clinical_Roper$number <- 1
clinical_Roper <- ddply(clinical_Roper,'ICB_response',transform,percent = 1/sum(number)*100)
ggplot(clinical_Roper)+
  scale_fill_manual(values = c("#475F77","#EBB481"))+ 
  geom_bar(aes(x=hlabel,fill=ICB_response),position = "fill")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))

##  External cohort survival
clinical_all <- data.frame(os_time = c(clinical_george$overall_survival..months.,
                                       clinical_GSE60052$os_time),
                           os_event = c(clinical_george$os,
                                        clinical_GSE60052$os),
                           hlabel= c(clinical_george$hlabel,
                                     clinical_GSE60052$hlabel))
surv <- survfit(Surv(as.numeric(os_time),os_event)~hlabel,data = clinical_all)
surv
survdiff(Surv(as.numeric(os_time),os_event)~hlabel,data = clinical_all)
summary(surv)
summary(coxph(Surv(as.numeric(os_time),os_event)~hlabel,data = clinical_all))
ggsurvplot(surv,
           pval = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette = c("#1c1e24","#A39391","#F06966"),
           legend = c(2,0.5), 
           legend.title = "", 
           conf.int = F,
           add.all = TRUE,
           legend.labs = c("All","ITHtyper low","ITHtyper high"), 
           xlab = "Time (months)",
           ylab = "Overall Survival"
)






















