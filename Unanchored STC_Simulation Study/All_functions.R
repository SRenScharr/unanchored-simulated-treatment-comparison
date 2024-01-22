# All functions used in the simulation study

##########cov function#############
# A function to sample covariates using NORTA algorithm (also called Normal copula)
cov_fun <-
  function(N ,
           char_cov,
           cor_input_AgD,
           cor_input_IPD,
           m1_AgD,
           sd1_AgD,
           m2_AgD,
           sd2_AgD,
           p1_AgD,
           p2_AgD,
           m1_IPD,
           sd1_IPD,
           m2_IPD,
           sd2_IPD,
           p1_IPD,
           p2_IPD) {
    if (char_cov == "bin") {
      # for AgD study
      myCop_AgD <-
        copula::normalCopula(cor_input_AgD, dim = 2, dispstr = "un")
      myMvd_AgD <-
        mvdc(
          copula = myCop_AgD,
          margins = c("binom", "binom"),
          paramMargins = list(list(size = 1, prob = p1_AgD),
                              list(size = 1, prob = p2_AgD))
        )
      
      # for IPD study
      myCop_IPD <-
        normalCopula(cor_input_IPD, dim = 2, dispstr = "un")
      myMvd_IPD <-
        mvdc(
          copula = myCop_IPD,
          margins = c("binom",  "binom"),
          paramMargins = list(list(size = 1, prob = p1_IPD),
                              list(size = 1, prob = p2_IPD))
        )
      
    }
    
    # AgD study
    mycov_AgD <- data.frame(rMvdc(N, myMvd_AgD))
    cov_AgD_A <- mycov_AgD
    cov_AgD_A$trt <- "A"
    
    cov_AgD_B <- mycov_AgD
    cov_AgD_B$trt <- "B"
    
    cov_AgD <- rbind(cov_AgD_A, cov_AgD_B)
    cov_AgD$pop <- "AgD"
    
    # IPD study
    mycov_IPD <- data.frame(rMvdc(N, myMvd_IPD))
    cov_IPD_A <- mycov_IPD
    cov_IPD_A$trt <- "A"
    
    cov_IPD_B <- mycov_IPD
    cov_IPD_B$trt <- "B"
    
    cov_IPD <- rbind(cov_IPD_A, cov_IPD_B)
    cov_IPD$pop <- "IPD"
    
    # combine AgD study and IPD study
    cov_all <- rbind(cov_AgD, cov_IPD)
    cov_all
  }

#########true treatment function###########
# A function to calcualte the true marginal treatment effect
true_trt_function <-  function(N ,
                               char_cov,
                               cor_input_AgD,
                               cor_input_IPD,
                               m1_AgD,
                               sd1_AgD,
                               m2_AgD,
                               sd2_AgD,
                               p1_AgD,
                               p2_AgD,
                               m1_IPD,
                               sd1_IPD,
                               m2_IPD,
                               sd2_IPD,
                               p1_IPD,
                               p2_IPD,
                               interaction_X2,
                               b_0,
                               b_1,
                               b_2,
                               b_trt_B,
                               b_X2_trt) {
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
  
  
  # simulate AgD and IPD data
  if (char_cov == "bin") {
    if (interaction_X2 == TRUE) {
      # with interaction
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
      # without interaction
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
  
  
  # simulated treated and non-treated AgD population
  AgD_all <- AB_IPD_all %>%
    filter(pop == "AgD")
  
  
  # obtain the marginal trt effect within AgD treated group
  fit_trueM_A <- glm(y ~ trt , data = AgD_all, family = binomial)
  true_trt <- coef(fit_trueM_A)[2]
  true_trt
}


#########Unanchored STC function#############
# A function for obtain ITC using the proposed unanchored STC
my_sim_fun <-
  function(N,
           n_bt,
           cov_all,
           covariates = c("X1", "X2"),
           char_cov,
           cor_input_AgD,
           cor_input_IPD,
           m1_AgD,
           sd1_AgD,
           m2_AgD,
           sd2_AgD,
           p1_AgD,
           p2_AgD,
           m1_IPD,
           sd1_IPD,
           m2_IPD,
           sd2_IPD,
           p1_IPD,
           p2_IPD,
           b_0,
           b_1,
           b_2,
           b_trt_B,
           interaction_X2,
           b_X2_trt) {
    # Sample covariates for IPD and AgD using Cov_function.R
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
    
    
    # simulate AgD and IPD data
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
    # non-treated arm of AgD population
    A_AgD <- AB_IPD_all %>%
      filter(trt == "A", pop == "AgD")
    
    # treated arm of IPD population
    B_IPD <- AB_IPD_all %>%
      filter(trt == "B", pop == "IPD")
    
    
    # Treated and non-treated AgD population
    AgD_all <- AB_IPD_all %>%
      filter(pop == "AgD")
    
    #bootstrap function
    test_boot_fun <- function(i) {
      index <- sample(nrow(B_IPD), N, replace = TRUE)
      boot_IPD <- B_IPD[index,]
      
      ipd_cor <-
        cor(boot_IPD[, c("X1", "X2")])
      
      #construct a normal copula for n_X random variables with correlation structure as ipd_cor.
      ipd_copula <-
        copula::normalCopula(copula::P2p(ipd_cor),
                             dim = ncol(ipd_cor),
                             dispstr = "un")
      
      if (char_cov == "bin") {
        Mvd <-
          mvdc(
            copula = ipd_copula,
            margins = c("binom", "binom"),
            paramMargins = list(list(
              size = 1, prob = mean(A_AgD$X1)
            ),
            list(
              size = 1, prob = mean(A_AgD$X2)
            ))
          )
        
      }
      
      # simulated covariates for AgD study
      AgD_cov <- rMvdc(10000, Mvd)
      AgD_cov <- data.frame(AgD_cov)
      
      full_cov <- c("X1", "X2")
      colnames(AgD_cov) <- full_cov
      
      # outcome regression model
      m_STC1 <-
        glm(as.formula(paste(
          "y ~ ",
          paste(covariates, collapse = "+")
        )), data = boot_IPD, family = binomial)
      
      # prediction
      p_STC1 <- predict(m_STC1,
                        newdata = data.frame(AgD_cov),
                        type = "response")
      
      # average treatment effect for AgD population treating with trt B
      p_B <- mean(p_STC1)
      
      t11 <- log(p_B / (1 - p_B))
      t11
    }
    
    # bootstrap to estimate mean and SE
    boot.t11 <- sapply(1:n_bt, test_boot_fun)
    
    t12 <- log(sum(A_AgD$y) / N / (1 - sum(A_AgD$y) / N))
    
    var_A <- N / (sum(A_AgD$y) * (N - sum(A_AgD$y)))
    
    c(mean(boot.t11) - t12, sqrt(var(boot.t11) + var_A))
    
  }