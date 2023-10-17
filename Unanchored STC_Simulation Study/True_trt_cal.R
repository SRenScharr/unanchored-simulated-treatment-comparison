###########run############
library(dplyr)
library(lmtest)
library(copula)
library(randtoolbox)
library(boot)
setwd("X:/HAR_WG/WG/Single-arm_Trial/ESTC/R code/ESTC Code Github")

source("All_functions.R")


#1. simulation setting and hyperparameters
for (ipt1 in 1) {#ipt1 in 1:10

true_trt <- rep(NA,12)

for (ipt2 in 1:12){
  #2. model input
  model_input_raw<-read.csv(paste("model_input_S",ipt1,"_BB.csv",sep=''))
  row.names(model_input_raw) <-model_input_raw[,1]
  model_input<-model_input_raw[,ipt2+1, drop=FALSE]
  
  #3.calculate the true treatment effects
  true_trt[ipt2]<- true_trt_function(N=5e6,
                               char_cov = model_input["char_cov",],
                               cor_input_AgD = as.numeric(model_input["cor_input_AgD",]),
                               cor_input_IPD = as.numeric(model_input["cor_input_IPD",]),
                               p1_AgD = as.numeric(model_input["p1_AgD",]),
                               p2_AgD = as.numeric(model_input["p2_AgD",]),
                               p1_IPD = as.numeric(model_input["p1_IPD",]),
                               p2_IPD = as.numeric(model_input["p2_IPD",]),
                               interaction_X2 = model_input["interaction_X2",],
                               b_0 = as.numeric(model_input["b_0",]),
                               b_1 = as.numeric(model_input["b_1",]),
                               b_2 = as.numeric(model_input["b_2",]),
                               b_trt_B = as.numeric(model_input["b_trt_B",]),
                               b_X2_trt = as.numeric(model_input["b_X2_trt",]))
  
}

write.csv(true_trt, file = paste("Results/true_trt_S",ipt1,".csv", sep=""), row.names = F)

}
