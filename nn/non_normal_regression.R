rm(list=ls())
library(grid)
library(mvtnorm)
library(stats)
library(MASS)
library(matrixcalc)
library(ggplot2)
printf <- function(...) invisible(print(sprintf(...)))

source("~/Research/ecology letters/rscripts/functions.R")
source("~/Research/ecology letters/rscripts/offsetmatching.R")

coeffsA1 <- NULL
coeffsA2 <- NULL
coeffsB1 <- NULL
coeffsB2 <- NULL
coeffsk  <- NULL
r2A1 <- NULL
r2A2 <- NULL
r2B1 <- NULL
r2B2 <- NULL
r2k  <- NULL

PTS <- 1000
OUTPUT <- NULL
INPUT <- NULL
pvals <- 10^(-3:-1)
for(P in pvals){
for(pts in 1:PTS){
  bool <- T # regenerate data until you get something non-normal, but doesn't explode
  while(bool){
    
    bb <- T
    while(bb){
      A <- runif(2,0,.01)
      B <- runif(2,0,.01)
      G <- rexp(2,1)
      n <- rexp(2,1/100)
      th <- rnorm(2,0,10)
      k <- rexp(1,0.1)
      m <- mu(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2])
      S <- SIGMA(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2],G[1],G[2])
      bb <- !is.positive.definite(S)
    }
    
    p <- 1
    while(p>.05){
      N <- rpois(1,20)
      while(N<3) N <- rpois(1,20)
      data <- rmvnorm(N,m,S)
      test1 <- shapiro.test(data[,1])
      test2 <- shapiro.test(data[,2])
      p1 <- test1$p.value
      p2 <- test2$p.value
      p <- min(p1,p2)
    }
    
    m <- colMeans(data)
    S <- (N-1)*var(data)/N
    MLS <- ml_sol(m[1],m[2],S[1,1],S[2,2],S[1,2],n[1],n[2],th[1],th[2],G[1],G[2])
    MLS <- c(MLS$A1,MLS$A2,MLS$B1,MLS$B2,MLS$k)
    ML_A <- MLS[1:2]
    ML_B <- MLS[3:4]
    ML_k <- MLS[5]
    bool <- F
    bool <- bool||OR(is.nan(MLS))||OR(is.infinite(MLS)) # screens for nans & infs
  }
  
  OUTPUT <- rbind(OUTPUT,c(ML_A,ML_B,sqrt(prod(ML_B))))
  INPUT  <- rbind(INPUT,c(A,B,sqrt(prod(B))))
  
}
  
  dfB1 <- data.frame(cbind(INPUT[,3],OUTPUT[,3]))
  dfB2 <- data.frame(cbind(INPUT[,4],OUTPUT[,4]))
  dfC <- data.frame(cbind(INPUT[,5],OUTPUT[,5]))
  wtsB1 <- 1/sqrt(dfB1$X1^2+dfB1$X2^2)
  wtsB2 <- 1/sqrt(dfB2$X1^2+dfB2$X2^2)
  wtsC <- 1/sqrt(dfC$X1^2+dfC$X2^2)
  lmB1 <- lm(X2~X1,data=dfB1,weights=wtsB1)
  lmB2 <- lm(X2~X1,data=dfB2,weights=wtsB2)
  
  # not general coeffs, but slopes, really  
  coeffsB1 <- c(coeffsB1,lmB1$coefficients[2])  
  coeffsB2 <- c(coeffsB2,lmB2$coefficients[2])  
  
  # r2
  r2B1 <- c(r2B1,summary(lmB1)$r.squared)  
  r2B2 <- c(r2B2,summary(lmB2)$r.squared)  
  
}


slpB1_nn <- data.frame(pvals,coeffsB1)
slpB2_nn <- data.frame(pvals,coeffsB2)
save(slpB1_nn,file="~/Research/ecology letters/data/nn/slopeB1_nn.Rda")
save(slpB2_nn,file="~/Research/ecology letters/data/nn/slopeB2_nn.Rda")
load("~/Research/ecology letters/data/nn/slopeB1_nn.Rda")

r2B1df_nn <- data.frame(pvals,r2B1)
r2B2df_nn <- data.frame(pvals,r2B2)
save(r2B1df_nn,file="~/Research/ecology letters/data/nn/r2B1_nn.Rda")
save(r2B2df_nn,file="~/Research/ecology letters/data/nn/r2B2_nn.Rda")
load("~/Research/ecology letters/data/nn/r2B1_nn.Rda")

slp_plt_nn <- ggplot(slpB1_nn,aes(x=pvals,y=coeffsB1))+
  geom_line(color="#00c970",size=2)+#geom_smooth(color="#00c970",method=lm,fullrange=T)+
  scale_x_continuous(trans = "log10")+
  xlab("Threshold of normality")+
  ylab("Slope of linear regression")+
  scale_y_continuous(breaks=seq(0,1.1,0.05),limits=c(1,1.1))+
  theme_bw()

slp_plt_nn
ggsave("~/Research/ecology letters/slp_nn_plt.eps",slp_plt_nn,width=5,height=2.5)

r2_plt_nn <- ggplot(r2B1df_nn,aes(x=pvals,y=r2B1))+
  geom_line(color="#1ec56f",size=2)+#geom_smooth(color="#1ec56f",method=lm,fullrange=T)+
  scale_x_continuous(trans = "log10")+
  xlab("Threshold of normality")+
  ylab(expression(paste(italic(R)^2," of linear regression")))+
  scale_y_continuous(breaks=seq(0,1,0.1),limits=c(.2,.4))+
  theme_bw()

r2_plt_nn
ggsave("~/Research/ecology letters/r2_nn_plt.eps",r2_plt_nn,width=5,height=2.5)
