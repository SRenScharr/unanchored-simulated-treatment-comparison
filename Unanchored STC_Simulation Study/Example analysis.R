# An example analysis based on one scenario setting from the simulation study
# S1: two covariates, both prognostic factors
#     small correlation between the two covariates within the study:0.2
#     large overlap both p_AgD=0.1, both p_IPD=0.2
#     lage treatment effect: -0.45

##The example analysis code has two components
##1. Data prep: set up the data used for this example analysis
##2. Unanchored STC with bootstrap: demonstrate the use of unanchored STC in practice


library(dplyr)
library(copula)

setwd("~") # set an appropriate working direction
source("All_functions.R")

###################################################################################
##1. Data prep
###################################################################################
N <- 200 #number of patients
covariates <- c("X1", "X2") #name of the covaraites
char_cov <- "bin" #binary covariates
cor_input_AgD <-
  0.2 #correlation between covariates in AgD population
cor_input_IPD <-
  0.2 #correlation between covariates in IPD population
p1_AgD <- 0.1
p2_AgD <- 0.1 #covariates in AgD population
p1_IPD <- 0.2
p2_IPD <- 0.2 #covariates in IPD population
b_0 <- -0.25
b_1 <- 0.09
b_2 <- 0.15 #model parameters
b_trt_B <-
  -0.45
interaction_X2 <- "FALSE"
b_X2_trt <- "NA" #model parameters
n_bt <- 10000 #number of bootstrap

# Sample covariate for IPD and AgD using cov_fun.R #
cov_all <- cov_fun(
  N = N,
  char_cov = char_cov,
  cor_input_AgD = cor_input_AgD,
  cor_input_IPD = cor_input_IPD,
  m1_AgD = m1_AgD,
  sd1_AgD = sd1_AgD,
  m2_AgD = m2_AgD,
  sd2_AgD = sd2_AgD,
  p1_AgD = p1_AgD,
  p2_AgD = p2_AgD,
  m1_IPD = m1_IPD,
  sd1_IPD = sd1_IPD,
  m2_IPD = m2_IPD,
  sd2_IPD = sd2_IPD,
  p1_IPD = p1_IPD,
  p2_IPD = p2_IPD
)

if (char_cov == "bin") {
  if (interaction_X2 == TRUE) {
    AB_IPD_all <-
      cov_all %>%
      # Generate outcomes using logistic model
      mutate(yprob = 1 / (1 + exp(-(
        b_0 + b_1 * X1 + b_2 * X2  +
          if_else(trt == "B", b_trt_B + b_X2_trt * X2 , 0)
      ))),
      y = rbinom(N * 4, 1, yprob)) %>%
      select(-yprob) # Drop the yprob column
  }
  
  
  if (interaction_X2 == FALSE) {
    AB_IPD_all <-
      cov_all %>%
      # Generate outcomes using logistic model
      mutate(yprob = 1 / (1 + exp(-(
        b_0 + b_1 * X1  + b_2 * X2   +
          if_else(trt == "B", b_trt_B , 0)
      ))),
      y = rbinom(N * 4, 1, yprob)) %>%
      select(-yprob) # Drop the yprob column
  }
}

# construct the data set IPD B arm and AgD A arm
# non-treated A arm of AgD population
A_AgD <- AB_IPD_all %>%
  filter(trt == "A", pop == "AgD")

r_AgD <- sum(A_AgD$y) #number of events in AgD population
p_x1_AgD <- mean(A_AgD$X1) #covariates X1 in AgD population
p_x2_AgD <- mean(A_AgD$X2) #covariates X2 in AgD population
summary_AgD <- c(r_AgD, N, p_x1_AgD, p_x2_AgD)
names(summary_AgD) <-
  c("number of events", "number of patients", "X1", "X2")

# data used for analysis related to AgD population
write.csv(summary_AgD, "A_AgDSummary.csv")

# treated B arm of IPD population
B_IPD <- AB_IPD_all %>%
  filter(trt == "B", pop == "IPD")

# data used for analysis related to IPD population
write.csv(B_IPD, "B_IPD.csv", row.names = FALSE)


###################################################################################
##2. Unanchored STC with bootstrap
###################################################################################
boot_fun <- function(i) {
  # bootstrap samples
  index <- sample(nrow(B_IPD), N, replace = TRUE)
  boot_IPD <- B_IPD[index,]
  
  # obtain the correlation structure for the IPD study
  ipd_cor <- cor(boot_IPD[, c("X1", "X2")])
  
  # simulate covariates for the AgD study
  ipd_copula <-
    copula::normalCopula(copula::P2p(ipd_cor),
                         dim = ncol(ipd_cor),
                         dispstr = "un") #construct a normal copula for n_X random variables with correlation structure as ipd_cor.
  
  
  if (char_cov == "bin") {
    Mvd <-
      mvdc(
        copula = ipd_copula,
        margins = c("binom", "binom"),
        paramMargins = list(list(size = 1, prob = p_x1_AgD),
                            list(size = 1, prob = p_x2_AgD))
      )
    
  }
  
  AgD_cov <- rMvdc(10000, Mvd)
  AgD_cov <- data.frame(AgD_cov)
  
  full_cov <- c("X1", "X2")
  colnames(AgD_cov) <- full_cov
  
  # outcome regression model
  m_STC1 <-
    glm(as.formula(paste("y ~ ",
                         paste(
                           covariates, collapse = "+"
                         ))), data = boot_IPD, family = binomial)
  
  # prediction
  p_STC1 <- predict(m_STC1,
                    newdata = data.frame(AgD_cov),
                    type = "response")
  
  # average treatment effect for AgD population treating with trt B
  p_B <- mean(p_STC1)
  t11 <- log(p_B / (1 - p_B))
  t11
}

# estimate mean and SE using bootstrap
boot.t11 <- sapply(1:n_bt, boot_fun)
t12 <- log(r_AgD / N / (1 - r_AgD / N))
var_A <- N / (r_AgD * (N - r_AgD)) #variance of A arm

#Mean and SE of the marginal treatment effect
c(mean(boot.t11) - t12, sqrt(var(boot.t11) + var_A)) 