###########run############
library(dplyr)
library(lmtest)
library(copula)
library(randtoolbox)
library(boot)

setwd("X:/HAR_WG/WG/Single-arm_Trial/ESTC/R code/Github/unanchored-simulated-treament-comparison/Unanchored STC_Simulation Study")
source("All_functions.R")

#######################################################################################
#1. We take the first scenario as an example#
Task_ID <- 1 #take value from 1 to 120, the scenario to run


#######################################################################################
#2. model input
n_mc <- 10 #number of monte carlo simulations
n_bt <- 100 #number of bootstrap
ipt1 <- ceiling(Task_ID/12) #Scenario number
ipt2 <- Task_ID-(ipt1-1)*12 #Study/Column number
model_input_raw<- read.csv(paste("Model_input_S",ipt1,"_BB.csv",sep=''))
row.names(model_input_raw) <-model_input_raw[,1]
model_input<-model_input_raw[,ipt2+1, drop=FALSE]

#######################################################################################
#3. run simulations, mean and SE saved 
final_result <- matrix(NA, nrow = n_mc, ncol = 2)
colnames(final_result) <- c("Mean","SE")
i <- 0
while (i < n_mc) {
  i <-i+1
   final_result[i, ] <-
     my_sim_fun(
       N = as.numeric(model_input["N",]),
       n_bt=n_bt,
       covariates = c("X1", "X2"),
       char_cov = as.character(model_input["char_cov",]),
       cor_input_AgD = as.numeric(model_input["cor_input_AgD",]),
       cor_input_IPD = as.numeric(model_input["cor_input_IPD",]),
       p1_AgD = as.numeric(model_input["p1_AgD",]),
       p2_AgD = as.numeric(model_input["p2_AgD",]),
       p1_IPD = as.numeric(model_input["p1_IPD",]),
      p2_IPD = as.numeric(model_input["p2_IPD",]),
      b_0 = as.numeric(model_input["b_0",]),
      b_1 = as.numeric(model_input["b_1",]),
      b_2 = as.numeric(model_input["b_2",]),
      b_trt_B = as.numeric(model_input["b_trt_B",]),
      interaction_X2 = model_input["interaction_X2",],
      b_X2_trt = as.numeric(model_input["b_X2_trt",])
    )
 }

tag<- 2 # to distinguish the output
write.csv(final_result, file = paste("Output/n",model_input["N",],"_S",ipt1,"_C",ipt2,"_mc",n_mc,"_",tag,".csv",sep=""), row.names = F)

#######################################################################################
