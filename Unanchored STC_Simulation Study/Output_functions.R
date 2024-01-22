##########Output function#############
# A function to calculate bias, model and empirical SE and covareage from simulations
output_fun <- function(final_result, M) {
  # assessing the performance for the treatment effect
  # obtain mean bias
  bias <- mean(final_result[, 2] - final_result[, 1])
  bias_mcse <-
    sqrt(1 / (M * (M - 1)) * sum((final_result[, 2] - mean(final_result[, 2])) ^
                                   2))
  
  #MSE
  MSE <- mean((final_result[, 2] - final_result[, 1]) ^ 2)
  MSE_mcse <-
    sqrt(sum(((
      final_result[, 2] - final_result[, 1]
    ) ^ 2 - MSE) ^ 2) / (M * (M - 1)))
  
  
  # empirical SE
  empSE <- sd(final_result[, 2])
  empSE_mcse <- empSE / sqrt(2 * (M - 1))
  
  # model SE
  modSE <- sqrt(mean((final_result[, 3]) ^ 2))
  modSE_mcse <-
    sqrt(var(final_result[, 3] ^ 2) / (4 * M * modSE ^ 2))
  
  # use the true trt effect obtained using the mean of each true_trt n=200
  coverage <- mean(final_result[, 4])
  coverage_mcse <- sqrt(coverage * (1 - coverage) / M)
  
  
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