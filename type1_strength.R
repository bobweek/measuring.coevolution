rm(list=ls())
require(mvtnorm)
require(matrixcalc)
require(ggplot2)

source("~/Research/ecology letters/rscripts/functions.R")

####################################################################

alpha <-0.05
Bs <- 10^(-6:-2)
type1s <- NULL

for(Bp in Bs){
  
  t1 <- NULL
  
  while(length(t1)<10000){
    bool <- T
    while(bool){
      
      bb <- T
      while(bb){
        bz <- rbinom(1,1,0.5)
        B <- (2*rbinom(2,1,0.5)-1)*Bp*c(bz,1-bz)
        if(Bp<0.01){
          A <- runif(2,0,.01)
          if(A[1] < -B[1]) A[1] <- runif(1,-B[1],.01) # we require Ai + Bi > 0
          if(A[2] < -B[2]) A[2] <- runif(1,-B[2],.01)
        }
        if(Bp==0.01){
          A <- runif(2,0,.01)
          if(A[1] == -B[1]) A[1] <- runif(1,-B[1],.011)
          if(A[2] == -B[2]) A[2] <- runif(1,-B[2],.011)
        } 
        G <- rexp(2,1)
        n <- rexp(2,1/100)
        th <- rnorm(2,0,10)
        k <- rexp(1,0.1)
        m <- mu(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2])
        S <- SIGMA(A[1],A[2],B[1],B[2],k,n[1],n[2],th[1],th[2],G[1],G[2])
        bb <- !is.positive.definite(S)
      }

      N <- rpois(1,20)
      while(N<3) N <- rpois(1,20)
      data <- rmvnorm(N,m,S)
      
      m <- colMeans(data)
      S <- (N-1)*var(data)/N
      MLS <- ml_sol(m[1],m[2],S[1,1],S[2,2],S[1,2],n[1],n[2],th[1],th[2],G[1],G[2])
      MLS <- c(MLS$A1,MLS$A2,MLS$B1,MLS$B2)
      ML_A <- MLS[1:2]
      ML_B <- MLS[3:4]
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
      t1 <- c(t1,detected)
    }
  } 
  
  type1s <- c( type1s, sum(t1)/length(t1) )
  
}

t1df <- data.frame(Bs,type1s)
save(t1df,file="~/Research/ecology letters/data/t1/t1_strength.Rda")
load("~/Research/ecology letters/data/t1/t1_strength.Rda")
t1_plt <- ggplot(t1df,aes(x=Bs,y=type1s))+
  geom_line(size=1,color="#1ec56f")+
  xlab("Strength of biotic selection")+
  ylab("False positive rate")+
  scale_x_continuous(trans='log10',breaks=c(1e-6,1e-4,1e-2))+
  scale_y_continuous(breaks=seq(0,0.25,0.125),limits=c(0,.25))+
  theme_bw()

t1_plt
ggsave("~/Research/ecology letters/t1_strength_plt.eps",t1_plt,width=5,height=2.5,device="eps")
