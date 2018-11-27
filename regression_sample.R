rm(list=ls())
library(grid)
library(mvtnorm)
library(stats)
library(MASS)
library(matrixcalc)
library(ggplot2)
require(olsrr)
printf <- function(...) invisible(print(sprintf(...)))

source("~/Research/ecology letters/rscripts/functions.R")

PTS <- 10000
NN <- c(seq(3,19,1),seq(20,50,5))
intA1 <- NULL
intA2 <- NULL
intB1 <- NULL
intB2 <- NULL
intC <- NULL
intk <- NULL
slpA1 <- NULL
slpA2 <- NULL
slpB1 <- NULL
slpB2 <- NULL
slpC <- NULL
slpk <- NULL
r2A1 <- NULL
r2A2 <- NULL
r2B1 <- NULL
r2B2 <- NULL
r2C <- NULL
r2k <- NULL

for(N in NN){
  
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
      
      data <- rmvnorm(N,m,S)
      
      m <- colMeans(data)
      S <- (N-1)*var(data)/N
      MLS <- ml_sol(m[1],m[2],S[1,1],S[2,2],S[1,2],n[1],n[2],th[1],th[2],G[1],G[2])
      MLS <- c(MLS$A1,MLS$A2,MLS$B1,MLS$B2,MLS$k)
      bool <- F
      bool <- bool||OR(is.nan(MLS))||OR(is.infinite(MLS)) # screens for nans & infs
    }
    
    OUTPUT <- rbind(OUTPUT,c(MLS))
    INPUT  <- rbind(INPUT,c(A,B,k))
    
  }
  
  # wls
  dfA1 <- data.frame(I=INPUT[,1],O=OUTPUT[,1])
  dfA2 <- data.frame(I=INPUT[,2],O=OUTPUT[,2])
  dfB1 <- data.frame(I=INPUT[,3],O=OUTPUT[,3])
  dfB2 <- data.frame(I=INPUT[,4],O=OUTPUT[,4])
  dfk  <- data.frame(I=INPUT[,5],O=OUTPUT[,5])
  wtsA1 <- 1/sqrt(dfA1$I^2+dfA1$O^2)
  wtsA2 <- 1/sqrt(dfA2$I^2+dfA2$O^2)
  wtsB1 <- 1/sqrt(dfB1$I^2+dfB1$O^2)
  wtsB2 <- 1/sqrt(dfB2$I^2+dfB2$O^2)
  wtsk  <- 1/sqrt(dfk$I^2+dfk$O^2)
  lmA1 <- lm(O~I,data=dfA1,weights=wtsA1)
  lmA2 <- lm(O~I,data=dfA2,weights=wtsA2)
  lmB1 <- lm(O~I,data=dfB1,weights=wtsB1)
  lmB2 <- lm(O~I,data=dfB2,weights=wtsB2)
  lmk <- lm(O~I,data=dfk,weights=wtsk)
  
  Cs <- data.frame(I=sqrt(abs(dfB1$I*dfB2$I)),O=sqrt(abs(dfB1$O*dfB2$O)))
  wts <- 1/sqrt(Cs$I^2+Cs$O^2)
  lmC <- lm(O~I,data=Cs,weights = wts)
  
  intA1 <- c(intA1,lmA1$coefficients[1])  
  intA2 <- c(intA2,lmA2$coefficients[1])  
  intB1 <- c(intB1,lmB1$coefficients[1])  
  intB2 <- c(intB2,lmB2$coefficients[1])  
  intk  <- c(intk, lmk$coefficients[1])  
  intC  <- c(intC, lmC$coefficients[1])
  
  slpA1 <- c(slpA1,lmA1$coefficients[2])  
  slpA2 <- c(slpA2,lmA2$coefficients[2])
  slpB1 <- c(slpB1,lmB1$coefficients[2])  
  slpB2 <- c(slpB2,lmB2$coefficients[2])  
  slpk  <- c(slpk, lmk$coefficients[2]) 
  slpC  <- c(slpC, lmC$coefficients[2])
  
  r2A1 <- c(r2A1,summary(lmA1)$r.squared)  
  r2A2 <- c(r2A2,summary(lmA2)$r.squared)  
  r2B1 <- c(r2B1,summary(lmB1)$r.squared)  
  r2B2 <- c(r2B2,summary(lmB2)$r.squared)  
  r2k  <- c(r2k, summary(lmk)$r.squared)  
  r2C  <- c(r2C, summary(lmC)$r.squared)  
  
}

intB1df <- data.frame(NN,intB1)
intB2df <- data.frame(NN,intB2)
intA1df <- data.frame(NN,intA1)
intA2df <- data.frame(NN,intA2)
intkdf <- data.frame(NN,intk)
intCdf <- data.frame(NN,intC)
save(intB1df,file="~/Research/ecology letters/data/reg_sample/intB1_sample.Rda")
save(intB2df,file="~/Research/ecology letters/data/reg_sample/intB2_sample.Rda")
save(intA1df,file="~/Research/ecology letters/data/reg_sample/intA1_sample.Rda")
save(intA2df,file="~/Research/ecology letters/data/reg_sample/intA2_sample.Rda")
save(intCdf,file="~/Research/ecology letters/data/reg_sample/intC_sample.Rda")
save(intkdf,file="~/Research/ecology letters/data/reg_sample/intk_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/intB1_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/intB2_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/intA1_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/intA2_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/intC_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/intk_sample.Rda")

slpB1df <- data.frame(NN,slpB1)
slpB2df <- data.frame(NN,slpB2)
slpA1df <- data.frame(NN,slpA1)
slpA2df <- data.frame(NN,slpA2)
slpkdf <- data.frame(NN,slpk)
slpCdf <- data.frame(NN,slpC)
save(slpB1df,file="~/Research/ecology letters/data/reg_sample/slopeB1_sample.Rda")
save(slpB2df,file="~/Research/ecology letters/data/reg_sample/slopeB2_sample.Rda")
save(slpA1df,file="~/Research/ecology letters/data/reg_sample/slopeA1_sample.Rda")
save(slpA2df,file="~/Research/ecology letters/data/reg_sample/slopeA2_sample.Rda")
save(slpkdf,file="~/Research/ecology letters/data/reg_sample/slopek_sample.Rda")
save(slpCdf,file="~/Research/ecology letters/data/reg_sample/slopeC_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/slopeB1_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/slopeB2_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/slopeA1_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/slopeA2_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/slopek_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/slopeC_sample.Rda")

r2B1df <- data.frame(NN,r2B1)
r2B2df <- data.frame(NN,r2B2)
r2A1df <- data.frame(NN,r2A1)
r2A2df <- data.frame(NN,r2A2)
r2kdf <- data.frame(NN,r2k)
r2Cdf <- data.frame(NN,r2C)
save(r2B1df,file="~/Research/ecology letters/data/reg_sample/r2B1_sample.Rda")
save(r2B2df,file="~/Research/ecology letters/data/reg_sample/r2B2_sample.Rda")
save(r2A1df,file="~/Research/ecology letters/data/reg_sample/r2A1_sample.Rda")
save(r2A2df,file="~/Research/ecology letters/data/reg_sample/r2A2_sample.Rda")
save(r2Cdf,file="~/Research/ecology letters/data/reg_sample/r2C_sample.Rda")
save(r2kdf,file="~/Research/ecology letters/data/reg_sample/r2k_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/r2B1_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/r2B2_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/r2A1_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/r2A2_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/r2k_sample.Rda")
load("~/Research/ecology letters/data/reg_sample/r2C_sample.Rda")

int_plt <- ggplot(intCdf,aes(x=NN,y=intC))+
  geom_line(size=1,color="#f8766d")+#geom_smooth(color="#f8766d",fullrange=T)+
  xlab("Sample size")+
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(-.1,1.1))+
  ylab("Intercept of linear regression")+
  theme_bw()

int_plt
ggsave("~/Research/ecology letters/int_sample_plt.eps",int_plt,width=5,height=2.5,device = "eps")

slp_plt <- ggplot(slpCdf,aes(x=NN,y=slpC))+
  geom_line(size=1,color="#f8766d")+#geom_smooth(color="#f8766d",fullrange=T)+
  xlab("Sample size")+
  xlim(5,50)+
  scale_y_continuous(breaks=seq(1,1.5,0.15),limits=c(1,1.3))+
  ylab("Slope of linear regression")+
  theme_bw()
  
slp_plt
ggsave("~/Research/ecology letters/slp_sample_plt.eps",slp_plt,width=5,height=2.5,device="eps")

r2_plt <- ggplot(r2Cdf,aes(x=NN,y=r2C))+
  geom_line(size=1,color="#00c970")+#geom_smooth(color="#00c970",fullrange=T)+
  xlab("Sample size")+
  xlim(5,50)+
  scale_y_continuous(breaks=seq(0,0.75,0.3),limits=c(0,0.6))+
  ylab(expression(paste(italic(R)^2," of linear regression")))+
  theme_bw()

r2_plt
ggsave("~/Research/ecology letters/r2_sample_plt.eps",r2_plt,width=5,height=2.5,device="eps")

