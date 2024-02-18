# ITHtyper
![image](https://github.com/ZhoulabCPH/ITHtyper/assets/143063392/5d791a48-bca6-4424-8c57-d12b002d0fe2)

# All codes depend on the R 4.3.0 platform with windows 11

# Overview
# The manuscript corresponding to the code is under review. We will update the summary as soon as the review is completed.


# Use of ITHtyper
# HCs-TME is represented as 1 in model code
# LCs-TME is represented as 0 in model code
# Key gene was "NKX1-2","TLE2","TPBG","SRSF6","DAZ4", "GPR31","CD274" ,"LYZ", "PCP4","ZIC1" ; 
# Please provide the expression profile of these ten genes. The rows are samples and the columns are genes.
# Please refer to "data/Training/train_data.csv" for the specific format.

setwd("Your/path")
Source(model training code-4.R)   
predict(mxgb4m, Yourdata ,type = "response")


# If you have any questions, please contact zhoumeng@wmu.edu.cn
