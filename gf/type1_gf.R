rm(list=ls())
require(mvtnorm)
require(matrixcalc)

source("~/Research/ecology letters/rscripts/functions.R")

####################################################################

alpha <-0.05
type1s <- NULL
GFs <- 10^(-5:-2)
for(GF in GFs){
  
  type1 <- NULL
  
  while(length(type1)<10000){
    
    bool <- T
    while(bool){
      bb <- T
      while(bb){
        A  <- runif(2,0,.01)
        B1 <- runif(1,0,.01)
        B2 <- runif(1,0,.01)
        b1 <- rbinom(1,1,0.5)
        b2 <- rbinom(1,1,0.5)
        B1 <- b1*b2*B1
        B2 <- b1*(1-b2)*B2
        B  <- c(B1,B2)
        d  <- rep(GF,2)
        G  <- rexp(2,1)
        n  <- rexp(2,1/100)
        th <- rnorm(2,0,10)        
        k  <- rexp(1,0.1)
        m  <- mu(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2])
        S  <- SIGMA_gf(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2],G[1],G[2],d[1],d[2])
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
    
    mu_ml <- m
    S_ml  <- S
    lnL   <- likelihood(data,mu_ml,S_ml)
    
    mu_n1 <- c(th[1],m[2])
    S_n1  <- S
    lnL1  <- likelihood(data,mu_n1,S_n1)
    
    mu_n2 <- c(m[1],th[2])
    S_n2  <- S
    lnL2  <- likelihood(data,mu_n2,S_n2)
    
    # chisqrs
    p1 <- 1-pchisq(2*(lnL-lnL1),1)
    p2 <- 1-pchisq(2*(lnL-lnL2),1)
    
    detected <- 0
    if(!is.nan(p1)&&!is.nan(p2)){
      if(p1<alpha&&p2<alpha) detected <- 1
      type1 <- c(type1,detected)
    }
    
  }
  
  type1s <- c( type1s, sum(type1)/length(type1) )
  
}

t1df_gf <- data.frame(GFs,type1s)
save(t1df_gf,file="~/Research/ecology letters/data/gf/t1_sample_gf.Rda")
load("~/Research/ecology letters/data/gf/t1_sample_gf.Rda")

t1_plt_gf <- ggplot(t1df_gf,aes(x=GFs,y=type1s))+
  geom_line(size=2,color="#f8766d")+#geom_smooth(color="#f8766d",method=lm,fullrange=T)+
  xlab("Gene flow")+
  scale_x_continuous(trans='log10')+
  scale_y_continuous(breaks=seq(0,1,0.05),limits=c(0,.1))+
  ylab("False positive rate")+
  theme_bw()

t1_plt_gf
ggsave("~/Research/ecology letters/t1_plt_gf.eps",t1_plt_gf,width=5,height=2.5,device="eps")
