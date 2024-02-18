#####  ROI WTA process  #####
library(ggridges)
library(ggplot2)
library(reshape2)
data <- read.csv("WTA_multi_ROIdata.csv",stringsAsFactors = F,row.names = 1)
data <- as.matrix(t(data))
roi <- read.csv("roi_id.csv",stringsAsFactors = F,row.names = 1)
roi <- roi[rownames(data),]
data <- data[rownames(roi),]
roi$all <-rep(1,dim(roi)[1])
data_pat <- data
rownames(data_pat) <- roi$code
data_sex <- data
rownames(data_sex) <- roi$Gender
data_loca <- data
rownames(data_loca) <- roi$TumorLocation
data_all <- data
rownames(data_all) <- roi$all
test_wide <- melt(data, varnames = c("ROI","gene"),value.name="exp")
test_wide <- as.data.frame(test_wide)
test_wide1 <- melt(data_pat, varnames = c("patient","gene"),value.name="exp")
test_wide1 <- as.data.frame(test_wide1)
test_wide2 <- melt(data_sex, varnames = c("sex","gene"),value.name="exp")
test_wide2 <- as.data.frame(test_wide2)
test_wide3 <- melt(data_loca, varnames = c("location","gene"),value.name="exp")
test_wide3 <- as.data.frame(test_wide3)
test_wide4 <- melt(data_all, varnames = c("all","gene"),value.name="exp")
test_wide4 <- as.data.frame(test_wide4)
test_wide1$patient <- factor(test_wide1$patient,
                             levels = c("Pt-100","Pt-115","Pt-125","Pt-162",
                                        "Pt-177","Pt-181","Pt-185","Pt-191",
                                        "Pt-192","Pt-214","Pt-227","Pt-237",
                                        "Pt-241","Pt-245","Pt-269","Pt-274",
                                        "Pt-278","Pt-289","Pt-294","Pt-323",
                                        "Pt-332","Pt-78","Pt-80","Pt-88",
                                        "Pt-97"))
ggplot(test_wide1, aes(x = exp, y = patient,fill = patient)) +
  geom_density_ridges(alpha=0.7,stat="binline", bins=33) + theme_ridges() +
  scale_x_continuous(breaks=c(0,4,8,12,16)) + 
  scale_fill_manual(values = c("#996699","#006699","#58B4AB","#E96463",
                               "#7A9D96","#8b220d","#c89c0e","#B89582",
                               "#492711","#FFCC99","#cbe86b","#bf5704",
                               "#E39183","#5B6A27","#EEDCC6","#ff85cb",
                               "#53bbf4","#C5C6B6","#ffc33c","#dee2d1",
                               "#f9fbba","#d4edf4","#57D1C9","#AC8697",
                               "#9baec8"))
ggplot(test_wide2, aes(x = exp, y = sex,fill = sex)) +
  geom_density_ridges(alpha=0.7,stat="binline", bins=33) + theme_ridges() +
  scale_x_continuous(breaks=c(0,4,8,12,16)) + 
  scale_fill_manual(values = c("#75c0b8","#e09a7d"))
ggplot(test_wide3, aes(x = exp, y = location,fill = location)) +
  geom_density_ridges(alpha=0.7,stat="binline", bins=33) + theme_ridges() +
  scale_x_continuous(breaks=c(0,4,8,12,16)) + 
  scale_fill_manual(values = c("#f3cfcb","#9d8fba"))
ggplot(test_wide, aes(x = exp, y = ROI,fill = ROI)) +
  geom_density_ridges(alpha=0.7,stat="binline", bins=50) + 
  theme_ridges() +
  scale_x_continuous(breaks=c(0,4,8,12,16))
