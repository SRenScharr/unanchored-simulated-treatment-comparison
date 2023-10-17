# function to save output from simulation
# run_output_function<-function(SA_name,N, M=1000){
#   mydata <- read.csv(paste("ESTC/Result/2 covs/n", N,"_", SA_name, ".csv",sep=""))
#   
#   output<-output_fun(mydata, M)
#   output
# }

output_fun <- function(final_result, M) {
  # assessing the performance for the treatment effect
  # obtain mean bias
  bias <- mean(final_result[, 2]-final_result[, 1])
  bias_mcse <-
    sqrt(1 / (M * (M - 1)) * sum((final_result[, 2] - mean(final_result[, 2])) ^
                                   2))
  
  #MSE
  MSE <- mean((final_result[, 2]-final_result[, 1])^2)
  MSE_mcse <- sqrt(sum(((final_result[,2]-final_result[,1])^2-MSE)^2)/(M*(M-1)))
    
    
  # empirical SE
  empSE <- sd(final_result[, 2])
  empSE_mcse <- empSE / sqrt(2 * (M - 1))
  
  # model SE
  modSE <- sqrt(mean((final_result[, 3]) ^ 2))
  modSE_mcse <- sqrt(var(final_result[, 3] ^ 2) / (4 * M * modSE ^ 2))
  
  # use the true trt effect obtained using the mean of each true_trt n=200
  coverage <- mean(final_result[, 4])
  coverage_mcse <- sqrt(coverage * (1 - coverage) / M)
  
  
  # bias in jaccard coef
  #bias_jacc <- mean(final_result[,5]-final_result[,4])
  
  # output
  output_trt <-
    c(
      bias,
      bias_mcse,
	  MSE,
	  MSE_mcse,
      empSE,
      empSE_mcse,
      modSE,
      modSE_mcse,
      coverage,
      coverage_mcse
    )
  
  output_trt
  
}