train <- read.csv("train_data.csv",stringsAsFactors = F,
                  row.names = 1,check.names = F)
###  XGboost  ###
library(xgboost)
library(Matrix)
library(verification)
traindata1 <- data.matrix(...)  ## 
traindata2 <- Matrix(as.matrix(...,sparse=T)  ## 
traindata3 <- factor(label,levels = c(0,1))   ### 
traindata4 <- list(data=traindata2,label=traindata3)  ### candidate training data
dtrain <- xgb.DMatrix(data = traindata4$data, label = as.character(traindata4$label))
mxgb4m <- xgboost(data = dtrain,   
                  objective='binary:logistic',
                  nround=100, 
                  nfold = 5,  
                  max_depth=5,  
                  subsample = 0.8,
                  colsample_bytree = 0.8,
                  eta=0.8,
                  eval_metric = "error")









  
