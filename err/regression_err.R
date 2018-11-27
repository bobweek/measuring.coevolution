rm(list=ls())
library(grid)
library(mvtnorm)
library(stats)
library(MASS)
library(matrixcalc)
library(ggplot2)
printf <- function(...) invisible(print(sprintf(...)))

source("~/Research/ecology letters/rscripts/functions.R")

PTS <- 1000
coeffsA1 <- NULL
coeffsA2 <- NULL
coeffsB1 <- NULL
coeffsB2 <- NULL
coeffsC <- NULL
r2A1 <- NULL
r2A2 <- NULL
r2B1 <- NULL
r2B2 <- NULL
r2C <- NULL
SDth <- seq(0.1,1,.1) # proportion of metapopulation variance
for(SD in SDth){
  
  OUTPUT <- NULL
  INPUT <- NULL
  for(pt in 1:PTS){
    bool <- T
    while(bool){
      
      bb <- T
      while(bb){
        A <- runif(2,0,0.01)
        B <- runif(2,0,0.01)
        G <- rexp(2,1)
        n <- rexp(2,1/100)
        th <- rnorm(2,0,10)
        k <-  rexp(1,0.1)
        m <-    mu(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2])
        S <- SIGMA(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2],G[1],G[2])
        bb <- !is.positive.definite(S)
      }
      
      N <- rpois(1,20)
      while(N<3) N <- rpois(1,20)
      data <- rmvnorm(N,m,S)
      
      m <- colMeans(data)
      S <- (N-1)*var(data)/N
      th_hat <- c(0,0)
      th_hat[1] <- rnorm(1,th[1],SD*S[1,1])
      th_hat[2] <- rnorm(1,th[2],SD*S[2,2])
      MLS <- ml_sol(m[1],m[2],S[1,1],S[2,2],S[1,2],n[1],n[2],th_hat[1],th_hat[2],G[1],G[2])
      MLS <- c(MLS$A1,MLS$A2,MLS$B1,MLS$B2,MLS$k)
      bool <- F
      bool <- bool||OR(is.nan(MLS))||OR(is.infinite(MLS)) # screens for nans & infs
    }
    
    OUTPUT <- rbind(OUTPUT,c(MLS[1:4],sqrt(prod(MLS[3:4]))))
    INPUT  <- rbind(INPUT,c(A,B,sqrt(prod(B))))
    
  }
  
  dfA1 <- data.frame(cbind(INPUT[,1],OUTPUT[,1]))
  dfA2 <- data.frame(cbind(INPUT[,2],OUTPUT[,2]))
  dfB1 <- data.frame(cbind(INPUT[,3],OUTPUT[,3]))
  dfB2 <- data.frame(cbind(INPUT[,4],OUTPUT[,4]))
  dfC <- data.frame(cbind(INPUT[,5],OUTPUT[,5]))
  wtsA1 <- 1/sqrt(dfA1$X1^2+dfA1$X2^2)
  wtsA2 <- 1/sqrt(dfA2$X1^2+dfA2$X2^2)
  wtsB1 <- 1/sqrt(dfB1$X1^2+dfB1$X2^2)
  wtsB2 <- 1/sqrt(dfB2$X1^2+dfB2$X2^2)
  wtsC <- 1/sqrt(dfC$X1^2+dfC$X2^2)
  lmA1 <- lm(X2~X1,data=dfA1,weights=wtsA1)
  lmA2 <- lm(X2~X1,data=dfA2,weights=wtsA2)
  lmB1 <- lm(X2~X1,data=dfB1,weights=wtsB1)
  lmB2 <- lm(X2~X1,data=dfB2,weights=wtsB2)
  lmC <- lm(X2~X1,data=dfC,weights=wtsC)
  
  # extract exponents (or slopes for log-log regressions)
  coeffsA1 <- c(coeffsA1,lmA1$coefficients[2])  
  coeffsA2 <- c(coeffsA2,lmA2$coefficients[2])  
  coeffsB1 <- c(coeffsB1,lmB1$coefficients[2]) 
  coeffsB2 <- c(coeffsB2,lmB2$coefficients[2])
  coeffsC <- c(coeffsC,lmC$coefficients[2])
  
  # extract percent variance explained
  r2A1 <- c(r2A1,summary(lmA1)$r.squared)  
  r2A2 <- c(r2A2,summary(lmA2)$r.squared)  
  r2B1 <- c(r2B1,summary(lmB1)$r.squared)  
  r2B2 <- c(r2B2,summary(lmB2)$r.squared)  
  r2C <- c(r2C,summary(lmC)$r.squared)  
  
}

slpB1_err <- data.frame(SDth,coeffsB1)
slpB2_err <- data.frame(SDth,coeffsB2)
slpA1_err <- data.frame(SDth,coeffsA1)
slpA2_err <- data.frame(SDth,coeffsA2)
slpC_err <- data.frame(SDth,coeffsC)
save(slpB1_err,file="~/Research/ecology letters/data/err/slopeB1_sample_err.Rda")
save(slpB2_err,file="~/Research/ecology letters/data/err/slopeB2_sample_err.Rda")
save(slpA1_err,file="~/Research/ecology letters/data/err/slopeA1_sample_err.Rda")
save(slpA2_err,file="~/Research/ecology letters/data/err/slopeA2_sample_err.Rda")
save(slpC_err,file="~/Research/ecology letters/data/err/slopeC_sample_err.Rda")
load("~/Research/ecology letters/data/err/slopeC_sample_err.Rda")

r2B1df_err <- data.frame(SDth,r2B1)
r2B2df_err <- data.frame(SDth,r2B2)
r2A1df_err <- data.frame(SDth,r2A1)
r2A2df_err <- data.frame(SDth,r2A2)
r2Cdf_err <- data.frame(SDth,r2C)
save(r2B1df_err,file="~/Research/ecology letters/data/err/r2B1_sample_err.Rda")
save(r2B2df_err,file="~/Research/ecology letters/data/err/r2B2_sample_err.Rda")
save(r2A1df_err,file="~/Research/ecology letters/data/err/r2A1_sample_err.Rda")
save(r2A2df_err,file="~/Research/ecology letters/data/err/r2A2_sample_err.Rda")
save(r2Cdf_err,file="~/Research/ecology letters/data/err/r2C_sample_err.Rda")
load("~/Research/ecology letters/data/err/r2C_sample_err.Rda")

slp_plt_err <- ggplot(slpC_err,aes(x=SDth,y=coeffsC))+
  geom_line(color="#1ec56f",size=2)+#geom_smooth(color="#1ec56f",method=lm,fullrange=T)+
  xlab("Standard deviation of error")+
  scale_y_continuous(breaks=seq(0,1.1,0.1),limits=c(.9,1.1))+
  ylab("Slope of linear regression")+
  theme_bw()

slp_plt_err
ggsave("~/Research/ecology letters/slp_plt_err.eps",slp_plt_err,width=5,height=2.5)

r2_plt_err <- ggplot(r2Cdf_err,aes(x=SDth,y=r2C))+
  geom_line(color="#1ec56f",size=2)+#geom_smooth(color="#1ec56f",method=lm,fullrange=T)+
  xlab("Standard deviation of error")+
  ylab(expression(paste(italic(R)^2," of linear regression")))+
  scale_y_continuous(breaks=seq(0,1,0.1),limits=c(.3,.5))+
  theme_bw()

r2_plt_err
ggsave("~/Research/ecology letters/r2_plt_err.eps",r2_plt_err,width=5,height=2.5)
