# To obtain the results of the performance assessment of the simulation study
library(dplyr)
library(copula)

#setwd("~") # set an appropriate working direction
source("Output_functions.R")

#1. simulation setting
for (ipt1 in 1:10) {
  n_mc <- rep(1000, 2)  # 2000 Monte Carlo simulations
  n_mc_l <- c(0, cumsum(n_mc))
  n_mc_all <- sum(n_mc)
  
  true_trt <-
    as.matrix(read.csv(paste(
      "Results/true_trt_S", ipt1, ".csv", sep = ""
    ), header = T))
  
  Output_S <- matrix(NA, nrow = 12, ncol = 10)
  for (ipt2 in 1:12) {
    #2. model input
    model_input_raw <-
      read.csv(paste("model_input_S", ipt1, "_BB.csv", sep = ''))
    row.names(model_input_raw) <- model_input_raw[, 1]
    model_input <- model_input_raw[, ipt2 + 1, drop = FALSE]
    
    #3. combine all the results
    final_result <- matrix(NA, nrow = n_mc_all, ncol = 4)
    colnames(final_result) <- c("true_trt", "STC", "sd", "coverage")
    final_result[, 1] <- as.numeric(true_trt[ipt2, ])
    
    for (i in 1:length(n_mc)) {
      final_result[(n_mc_l[i] + 1):(n_mc_l[i + 1]), 2:3] <-
        as.matrix(read.csv(
          paste(
            "Output/n",
            model_input["N", ],
            "_S",
            ipt1,
            "_C",
            ipt2,
            "_mc",
            n_mc[i],
            "_",
            i,
            ".csv",
            sep = ""
          )
        ))
    }
    
    
    for (i in 1:n_mc_all)
      final_result[i, 4] <-
      ifelse(
        final_result[i, 1] > final_result[i, 2] + 1.96 * final_result[i, 3]
        ||
          final_result[i, 1] < final_result[i, 2] - 1.96 * final_result[i, 3],
        0,
        1
      )
    
    #4. final output
    Output_S[ipt2, ] <- output_fun(final_result, M = n_mc_all)
    
  }
  
  colnames(Output_S) <-
    c(
      "bias",
      "bias_mcse",
      "MSE",
      "MSE_mcse",
      "empSE",
      "empSE_mcse",
      "modSE",
      "modSE_mcse",
      "coverage",
      "coverage_mcse"
    )
  row.names(Output_S) <- colnames(model_input_raw)[-1]
  
  # save the results
  write.csv(Output_S, file = paste("Results/Output_S", ipt1, ".csv", sep = ""))
  
}
