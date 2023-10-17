setwd("X:/HAR_WG/WG/Single-arm_Trial/ESTC/R code/ESTC Code Github")

library(SimCorrMix)
library(copula)
library(faux) # KR found a new package does NORTA and intermediate correlation transformatio

char_cov <- "bin"
N <- 5e6
p_AgD <- c(0.1, 0.9, 0.8, 0.5)
p_IPD <- c(0.2, 0.8, 0.2, 0.6)

true_corr <- function(N,p1, p2, cor_input) {
  myCop <- normalCopula(cor_input, dim = 2, dispstr = "un")
  myMvd_AgD <-
    mvdc(
      copula = myCop,
      margins = c("binom", "binom"),
      paramMargins = list(list(size = 1, prob = p1),
                          list(size = 1, prob = p2))
    )
  mycov_AgD <- data.frame(rMvdc(N, myMvd_AgD))
  cor(mycov_AgD)[1,2]
}

x1<-sapply(1:4, function(i){true_corr(N,p_AgD[i], p_AgD[i],0.1)})
x2<-sapply(1:4, function(i){true_corr(N,p_IPD[i], p_IPD[i],0.1)})
x3<-sapply(1:4, function(i){true_corr(N,p_AgD[i], p_AgD[i],0.9)})
x4<-sapply(1:4, function(i){true_corr(N,p_IPD[i], p_IPD[i],0.9)})

output<-cbind(x1,x2,x3,x4)
colnames(output) <- c("low_AgD","low_IPD","high_AgD","high_IPD")

write.csv(output,"True correlation.csv", row.names = F)
