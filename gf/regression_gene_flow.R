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
slpC <- NULL
r2A1 <- NULL
r2A2 <- NULL
r2B1 <- NULL
r2B2 <- NULL
r2C <- NULL
GFs <- 10^(-5:-2)
for(GF in GFs){
  
  OUTPUT <- NULL
  INPUT <- NULL
  for(pt in 1:PTS){
    bool <- T
    while(bool){
      
      bb <- T
      while(bb){
        A <- runif(2,0,0.01)
        B <- runif(2,0,0.01)
        d <- rep(GF,2) # d for dispersion
        G <- rexp(2,1)
        n <- rexp(2,1/100)
        th <- rnorm(2,0,10)
        k <- rexp(1,0.1)
        m <- mu(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2])
        S <- SIGMA_gf(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2],G[1],G[2],d[1],d[2])
        bb <- !is.positive.definite(S)
      }
     
      N <- rpois(1,20)
      while(N<3) N <- rpois(1,20)
      data <- rmvnorm(N,m,S)

      m <- colMeans(data)
      S <- (N-1)*var(data)/N
      ML_S <- ml_sol(m[1],m[2],S[1,1],S[2,2],S[1,2],n[1],n[2],th[1],th[2],G[1],G[2])
      ML_S <- c(ML_S$A1,ML_S$A2,ML_S$B1,ML_S$B2)
      bool <- F
      bool <- bool||OR(is.nan(ML_S))||OR(is.infinite(ML_S)) # screens for nans & infs
    }
    
    OUTPUT <- rbind(OUTPUT,ML_S)
    INPUT  <- rbind(INPUT,c(A,B))
    
  }

  dfB1 <- data.frame(cbind(INPUT[,3],OUTPUT[,3]))
  dfB2 <- data.frame(cbind(INPUT[,4],OUTPUT[,4]))
  Cs <- data.frame(I=sqrt(abs(dfB1$X1*dfB2$X1)),O=sqrt(abs(dfB1$X2*dfB2$X2)))
  wts <- 1/sqrt(Cs$I^2+Cs$O^2)
  lmC <- lm(O~I,data=Cs,weights = wts)
  
  slpC     <- c(slpC, lmC$coefficients[2])
  
  r2C  <- c(r2C, summary(lmC)$r.squared)  
  
}

slpC_gf <- data.frame(GFs,slpC)
save(slpC_gf, file="~/Research/ecology letters/data/gf/slopeC_sample_gf.Rda")
load("~/Research/ecology letters/data/gf/slopeC_sample_gf.Rda")

r2Cdf_gf <- data.frame(GFs,r2C)
save(r2Cdf_gf, file="~/Research/ecology letters/data/gf/r2C_sample_gf.Rda")
load("~/Research/ecology letters/data/gf/r2C_sample_gf.Rda")

slp_plt_gf <- ggplot(slpC_gf,aes(x=GFs,y=slpC))+
  geom_line(size=2,color="#1ec56f")+#geom_smooth(color="#1ec56f",method=lm,fullrange=T)+
  xlab("Gene flow")+
  scale_x_continuous(trans='log10')+
  ylab("Slope of linear regression")+
  scale_y_continuous(breaks=seq(0,1.5,0.2),limits=c(0.8,1.2))+
  theme_bw()

slp_plt_gf
ggsave("~/Research/ecology letters/slp_plt_gf.eps",slp_plt_gf,width=5,height=2.5,device="eps")

r2_plt_gf <- ggplot(r2Cdf_gf,aes(x=GFs,y=r2C))+
  geom_line(size=2,color="#1ec56f")+#geom_smooth(color="#1ec56f",method=lm,fullrange=T)+
  xlab("Gene flow")+
  scale_x_continuous(trans='log10')+
  ylab(expression(paste(italic(R)^2," of linear regression")))+
  scale_y_continuous(breaks=seq(0,1,0.25),limits=c(0,0.5))+
  theme_bw()

r2_plt_gf
ggsave("~/Research/ecology letters/r2_plt_gf.eps",r2_plt_gf,width=5,height=2.5,device="eps")
