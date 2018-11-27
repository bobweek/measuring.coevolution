rm(list=ls())
require(mvtnorm)
require(matrixcalc)

source("~/Research/ecology letters/rscripts/functions.R")

####################################################################

alpha <-0.05
powers <- NULL
SDth <- seq(0.1,1,.1) # proportion of metapopulation variance
for(SD in SDth){
  
  power <- NULL
  
  while(length(power)<1000){
    
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
      power <- c(power,detected)
    }
    
  }
  
  powers <- c( powers, sum(power)/length(power) )
  
}

pwrdf_err <- data.frame(SDth,powers)
save(pwrdf_err,file="~/Research/ecology letters/data/err/pwr_sample_err.Rda")
load("~/Research/ecology letters/data/err/pwr_sample_err.Rda")

pwr_plt_err <- ggplot(pwrdf_err,aes(x=SDth,y=powers))+
  geom_line(color="#f8766d",size=2)+#geom_smooth(color="#f8766d",method=lm,fullrange=T)+
  xlab("Standard deviation of error")+
  scale_y_continuous(breaks=seq(0,1,0.1),limits=c(.8,1))+
  ylab("Power")+
  theme_bw()

pwr_plt_err
ggsave("~/Research/ecology letters/pwr_plt_err.eps",pwr_plt_err,width=5,height=2.5)
