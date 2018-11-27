rm(list=ls())
require(mvtnorm)
require(matrixcalc)
require(ggplot2)

source("~/Research/ecology letters/rscripts/functions.R")

####################################################################

alpha <-0.05
Ns <- c(seq(3,17,2),seq(20,50,5))
type1s <- NULL

for(N in Ns){
  
  type1 <- NULL
  
  while(length(type1)<10000){
    
    bool <- T
    while(bool){
      bb <- T
      while(bb){
        A <- runif(2,0,.01)
        B <- runif(2,0,.01)*(2*rbinom(2,1,.5)-1)
        b <- rbinom(2,1,0.5)
        B[1] <- b[1]*b[2]*B[1]
        B[2] <- b[1]*(1-b[2])*B[2]
        if(A[1] < -B[1]) B[1] <- -runif(1,0,A[1])
        if(A[2] < -B[2]) B[2] <- -runif(1,0,A[2])
        G <- rexp(2,1)
        n <- rexp(2,1/100)
        th <- rnorm(2,0,10)
        k <- rexp(1,0.1)
        m <- mu(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2])
        S <- SIGMA(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2],G[1],G[2])
        bb <- !is.positive.definite(S)
      }

      data <- rmvnorm(N,m,S)
      
      m <- colMeans(data)
      S <- (N-1)*var(data)/N
      MLS <- ml_sol(m[1],m[2],S[1,1],S[2,2],S[1,2],n[1],n[2],th[1],th[2],G[1],G[2])
      MLS <- c(MLS$A1,MLS$A2,MLS$B1,MLS$B2)
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
      type1 <- c(type1,detected)
    }
    
  }
  
  type1s <- c( type1s, sum(type1)/length(type1) )

}

t1df <- data.frame(Ns,type1s)
save(t1df,file="~/Research/ecology letters/data/t1/t1_sample.Rda")
load("~/Research/ecology letters/data/t1/t1_sample.Rda")
t1_plt <- ggplot(t1df,aes(x=Ns,y=type1s))+
  geom_line(size=1,color="#f8766d")+
  xlab("Sample size")+
  ylab("False positive rate")+
  xlim(5,50)+
  scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0,0.2))+
  theme_bw()

t1_plt
ggsave("~/Research/ecology letters/t1_sample_plt.eps",t1_plt,width=5,height=2.5,device="eps")
