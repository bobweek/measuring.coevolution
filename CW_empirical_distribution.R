rm(list=ls())
library(stats)
library(MASS)
library(matrixcalc)
library(parallel)
library(foreach)
library(doParallel)
no_cores <- detectCores()-1
printf <- function(...) invisible(print(sprintf(...)))

#         1   2   3  4 5  6  7  8  9  10  11 12 13 14 15
# P <- c(mu1,mu2,V1,V2,C,m1,m2,G1,G2,th1,th2,n1,n2,k1,k2) DATA + BACKGROUND PARAMETERS

# this script does not include migration

Statistic <- function(D,mu_0,V_0,mu_1,V_1){
  N <- nrow(D)
  c <- N*log(det(V_0)/det(V_1))
  c0 <- 0
  for( i in 1:N){
    c0 <- c0 + mahalanobis(D[i,],mu_0,V_0)
  }
  c1 <- 0
  for( i in 1:N){
    c1 <- c1 + mahalanobis(D[i,],mu_1,V_1)
  }
  c <- c+c0[1]-c1[1]
  return(c)
}

Fhat <- function(x,L){
  N <- length(L)
  L0 <- NULL
  for(i in 1:N){
    if(L[i]<=x) L0 <- c(L0,L[i])
  }
  f <- length(L0)/N
  return(f)
}


SIGMAM <- function(X,P){
  y <- matrix(0,2,2);
  y[1,1] <- P[9]*P[13]*(X[2]+X[4])^2 + P[8]*( P[12]*X[1]^2 + P[13]*( X[1]*X[4]+X[3]*X[4]+X[2]*X[3] ) )
  y[2,2] <- P[8]*P[12]*(X[1]+X[3])^2 + P[9]*( P[13]*X[2]^2 + P[12]*( X[2]*X[3]+X[3]*X[4]+X[1]*X[4] ) )
  y[1,2] <- P[8]*P[12]*X[1]*(X[1]+X[3]) + P[9]*P[13]*X[2]*(X[2]+X[4])
  y[2,1] <- y[1,2];
  y <- y/( 2*P[12]*P[13]*( X[1]*X[4]+X[2]*X[3]+X[3]*X[4] )*( P[8]*(X[1]+X[3])+P[9]*(X[2]+X[4]) ) )
  return(y)
}

#         1   2   3  4 5  6  7  8  9  10  11 12 13 14 15
# P <- c(mu1,mu2,V1,V2,C,m1,m2,G1,G2,th1,th2,n1,n2,k1,k2) DATA + BACKGROUND PARAMETERS
# X <- c( B1, B2,A1,A2)

# with offsets set to differences
muM <- function(X,P){
  y <- numeric(2);
  y[1] <- (X[3]*X[4]*P[10] + X[2]*X[3]*P[10] + X[1]*X[4]*P[11] + sqrt(X[1])*X[4] + sqrt(X[1])*X[2] + X[1]*sqrt(X[2]) )
  y[2] <- (X[3]*X[4]*P[11] + X[1]*X[4]*P[11] + X[2]*X[3]*P[10] + sqrt(X[2])*X[3] + sqrt(X[1])*X[2] + X[1]*sqrt(X[2]) )
  y <- y/(X[3]*X[4] + X[2]*X[3] + X[1]*X[4])
  return(y)
}

female_snout <- c(9.13,10.42,9.89,10.31,10.05,10.66,9.12,9.21,9.61,13.63,10.06,12.98,11.68,11.48,14.54,
                  12.99,16.92,21.11,18.28,20.59)
camellia_coincidence <- c(3.88,6.07,6.30,6.13,6.66,7.73,6.76,6.42,7.52,11.65,7.77,12.80,11.89,11.13,12.49,
                          17.83,12.97,20.41,21.21,19.35)
data <- data.frame(cbind(female_snout,camellia_coincidence))
female_snout_var <- c(0.72,0.61,1.02,0.83,0.87,1.17,0.74,0.69,0.59,0.91,0.66,0.94,0.70,1.87,1.14,
                      1.87,1.55,1.85,1.61,0.46)^2
camellia_coincidence_var <- c(0.59,0.87,1.65,1.20,1.38,2.47,1.14,1.08,1.47,2.82,1.69,2.45,2.61,2.68,3.39,
                              1.80,3.60,3.99,2.38,2.68)^2
camellia_h <- 0.74 # the mean heritability reported
G1 <- 0.5*mean(female_snout_var)
G2 <- camellia_h*mean(camellia_coincidence_var)
m1 <- 10^-4
m2 <- 10^-4
camellia_absence <- c(5.22,5.27,5.48,7.54,4.68,5.19,4.91,4.99,5.55,7.65,6.82,6.42,7.63,7.69,7.28,7.56,6.61)
male_snout <- c(6.19,5.57,5.5,5.45,5.52,5.71,6.01,5.95,6.33,6.42,6.38)
k1 <- 0
k2 <- 0
th_hat_weevil <- mean(male_snout)+k1
th_hat_camellia <- mean(camellia_absence)+k2
th1 <- th_hat_weevil
th2 <- th_hat_camellia
n_weevil   <- c(26389,41944,19167,66389,44167,23333,30278)
n_camellia <- c(1178,1334,2155,1841,3112,2389,1770)
n1 <- 1/mean(1/n_weevil)
n2 <- 1/mean(1/n_camellia)
N <- nrow(data)
n <- N
P <- c(colMeans(data),(N-1)*var(female_snout)/N,(N-1)*var(camellia_coincidence)/N,(N-1)*var(data)[1,2]/N,m1,m2,G1,G2,th1,th2,n1,n2,k1,k2)
TRY <- 1000
seeds <- lapply(1:TRY, function(x) c(runif(4,-30,-5)))
cl <- makeCluster(no_cores,type="FORK")
registerDoParallel(cl)

# finding the ML solution
clusterEvalQ(cl,{
  Lrat <- function(X){
    X <- exp(X)
    N  <- nrow(data)
    mH <- colMeans(data)
    sH <- (N-1)*var(data)/N
    s  <- SIGMAM(X,P)
    out <- tryCatch(solve(s), error = function(e) e )
    while( is.singular.matrix(s) || any(class(out) == "error") ){
      X <- exp(runif(4,-30,-5))
      s  <- SIGMAM(X,P)
      out <- tryCatch(solve(s), error = function(e) e )
    }
    m  <- muM(X,P)
    c  <- N*log(det(s)/det(sH))
    for(i in 1:N){
      c <- c + mahalanobis(data[i,],m,s)
    }
    for(i in 1:N){
      c <- c - mahalanobis(data[i,],mH,sH)
    }
    return(c)
  }
})
solns <- parLapply(cl, seeds, function(x) c( optim(f=Lrat,par=x) ))
values <- NULL
for(i in 1:TRY){
  values <- c(values,solns[[i]]$value)
}
sol <- solns[[ which.min(values) ]]
xml <- muM(exp(sol$par),P)
Sml <- SIGMAM(exp(sol$par),P)
cat("\n************* Estimates of Selection Strengths *************\n")
B1 <- exp(sol$par)[1]
B2 <- exp(sol$par)[2]
A1 <- exp(sol$par)[3]
A2 <- exp(sol$par)[4]
printf("B1 = %0.3g,  B2 = %0.3g,  A1 = %0.3g,  A2 = %0.3g",B1,B2,A1,A2)

# find ml parameters for null1: B1=0
clusterEvalQ(cl,{
  Lrat <- function(X){
    X  <- exp(X)
    Y  <- c(0,X)
    N  <- nrow(data)
    mH <- colMeans(data)
    sH <- (N-1)*var(data)/N
    s  <- SIGMAM(Y,P)
    out <- tryCatch(solve(s), error = function(e) e )
    while( is.singular.matrix(s) || any(class(out) == "error") ){
      X <- exp(runif(3,-30,-5))
      Y  <- c(0,X)
      s  <- SIGMAM(Y,P)
      out <- tryCatch(solve(s), error = function(e) e )
    }
    m  <- muM(Y,P)
    c  <- N*log(det(s)/det(sH))
    for(i in 1:N){
      c <- c + mahalanobis(data[i,],m,s)
    }
    for(i in 1:N){
      c <- c - mahalanobis(data[i,],mH,sH)
    }
    return(c)
  }
})
solns <- parLapply(cl, seeds, function(x) c( optim(f=Lrat,par=x) ))
values <- NULL
for(i in 1:TRY){
  values <- c(values,solns[[i]]$value)
}
sol1 <- solns[[ which.min(values) ]]
Y <- exp(sol1$par)
Y <- c(0,Y)
mu1 <- muM(Y,P)
sigma1 <- SIGMAM(Y,P)
cat("\n************* Estimates of Selection Strengths Under Null 1 *************\n")
B1 <- 0
B2 <- exp(sol1$par)[1]
A1 <- exp(sol1$par)[2]
A2 <- exp(sol1$par)[3]
printf("B1 = %0.3g,  B2 = %0.3g,  A1 = %0.3g,  A2 = %0.3g",B1,B2,A1,A2)


# find ml parameters for null2: B2=0
clusterEvalQ(cl,{
  Lrat <- function(X){
    X  <- exp(X)
    Y <- c(X[1],0,X[2:3])
    N  <- nrow(data)
    mH <- colMeans(data)
    sH <- (N-1)*var(data)/N
    s  <- SIGMAM(Y,P)
    out <- tryCatch(solve(s), error = function(e) e )
    while( is.singular.matrix(s) || any(class(out) == "error") ){
      X <- exp(runif(3,-30,-5))
      Y <- c(X[1],0,X[2:3])
      s  <- SIGMAM(Y,P)
      out <- tryCatch(solve(s), error = function(e) e )
    }
    m  <- muM(Y,P)
    c  <- N*log(det(s)/det(sH))
    for(i in 1:N){
      c <- c + mahalanobis(data[i,],m,s)
    }
    for(i in 1:N){
      c <- c - mahalanobis(data[i,],mH,sH)
    }
    return(c)
  }
})
solns <- parLapply(cl, seeds, function(x) c( optim(f=Lrat,par=x) ))
values <- NULL
for(i in 1:TRY){
  values <- c(values,solns[[i]]$value)
}
sol2 <- solns[[ which.min(values) ]]
Y <- exp(sol2$par)
Y <- c(Y[1],0,Y[2:3])
mu2 <- muM(Y,P)
sigma2 <- SIGMAM(Y,P)
cat("\n************* Estimates of Selection Strengths Under Null 2 *************\n")
B1 <- exp(sol2$par)[1]
B2 <- 0
A1 <- exp(sol2$par)[2]
A2 <- exp(sol2$par)[3]
printf("B1 = %0.3g,  B2 = %0.3g,  A1 = %0.3g,  A2 = %0.3g",B1,B2,A1,A2)


# find ml parameters for null3: B1=B2=0
clusterEvalQ(cl,{
  Lrat <- function(X){
    X  <- exp(X)
    Y <- c(0,0,X)
    N  <- nrow(data)
    mH <- colMeans(data)
    sH <- (N-1)*var(data)/N
    s  <- SIGMAM(Y,P)
    out <- tryCatch(solve(s), error = function(e) e )
    while( is.singular.matrix(s) || any(class(out) == "error") ){
      X <- exp(runif(2,-30,-5))
      Y <- c(0,0,X)
      s  <- SIGMAM(Y,P)
      out <- tryCatch(solve(s), error = function(e) e )
    }
    m  <- muM(Y,P)
    c  <- N*log(det(s)/det(sH))
    for(i in 1:N){
      c <- c + mahalanobis(data[i,],m,s)
    }
    for(i in 1:N){
      c <- c - mahalanobis(data[i,],mH,sH)
    }
    return(c)
  }
})
solns <- parLapply(cl, seeds, function(x) c( optim(f=Lrat,par=x) ))
values <- NULL
for(i in 1:TRY){
  values <- c(values,solns[[i]]$value)
}
sol3 <- solns[[ which.min(values) ]]
Y <- exp(sol3$par)
Y <- c(0,0,Y)
mu3 <- muM(Y,P)
sigma3 <- SIGMAM(Y,P)
cat("\n************* Estimates of Selection Strengths Under Null 3 *************\n")
B1 <- 0
B2 <- 0
A1 <- exp(sol3$par)[1]
A2 <- exp(sol3$par)[2]
printf("B1 = %0.3g,  B2 = %0.3g,  A1 = %0.3g,  A2 = %0.3g",B1,B2,A1,A2)


#
# comparison with null1
#

cat("\n************* Inference Under Null 1 *************\n")
POINTS <- 250
TRY <- 30
seeds <- lapply(1:TRY, function(x) c(runif(4,-30,-5)))
x <- mu1
S <- sigma1
clusterExport(cl,c("x","S"))
lr1 <- Statistic(as.matrix(data),x,S,xml,Sml)
printf("LR-statistic for null1: %0.3f",lr1)
L1 <- NULL
while(length(L1) < POINTS){

  DATA <- mvrnorm(N,x,S)
  out <- tryCatch(solve( var(DATA) ), error = function(e) e )
  while( is.singular.matrix( var(DATA) ) || any(class(out) == "error") ){
    DATA <- mvrnorm(N,x,S)
    out <- tryCatch(solve( var(DATA) ), error = function(e) e )
  }
  clusterExport(cl,c("DATA"))

  clusterEvalQ(cl,{
    Lrat <- function(X){
      X <- exp(X)
      s  <- SIGMAM(X,P)
      m  <- muM(X,P)
      N  <- nrow(DATA)
      mH <- colMeans(DATA)
      sH <- (N-1)*var(DATA)/N
      c  <- N*log(det(s)/det(sH))
      for(i in 1:N){
        c <- c + mahalanobis(DATA[i,],m,s)
      }
      for(i in 1:N){
        c <- c - mahalanobis(DATA[i,],mH,sH)
      }
      return(c)
    }
  })
  solns <- foreach(i=seeds,.errorhandling = 'remove')%dopar%optim(f=Lrat,par=i)
  values <- NULL
  for(i in 1:length(solns)){
    values <- c(values,solns[[i]]$value)
  }
  min <- which.min(values)
  sol0 <- solns[[min]]
  x0 <- muM(exp(sol0$par),P)
  S0 <- SIGMAM(exp(sol0$par),P)
  L1 <- c(L1,Statistic(as.matrix(DATA),x,S,x0,S0))
}
L1 <- L1[L1 > 0]
Ln <- length(L1)
printf("Probability of rejecting null1: %0.3f",Fhat(lr1,L1))
printf("%.1f%% points lost from empirical distribution",100-100*Ln/POINTS)

#
# comparison with null2
#

cat("\n************* Inference Under Null 2 *************\n")
x <- mu2
S <- sigma2
clusterExport(cl,c("x","S"))
lr2 <- Statistic(as.matrix(data),x,S,xml,Sml)
printf("LR-statistic for null2: %0.3f",lr2)
L2 <- NULL
while(length(L2) < POINTS){

  DATA <- mvrnorm(N,x,S)
  out <- tryCatch(solve( var(DATA) ), error = function(e) e )
  while( is.singular.matrix( var(DATA) ) || any(class(out) == "error") ){
    DATA <- mvrnorm(N,x,S)
    out <- tryCatch(solve( var(DATA) ), error = function(e) e )
  }
  clusterExport(cl,c("DATA"))

  clusterEvalQ(cl,{
    Lrat <- function(X){
      X <- exp(X)
      s  <- SIGMAM(X,P)
      m  <- muM(X,P)
      N  <- nrow(DATA)
      mH <- colMeans(DATA)
      sH <- (N-1)*var(DATA)/N
      c  <- N*log(det(s)/det(sH))
      for(i in 1:N){
        c <- c + mahalanobis(DATA[i,],m,s)
      }
      for(i in 1:N){
        c <- c - mahalanobis(DATA[i,],mH,sH)
      }
      return(c)
    }
  })
  solns <- foreach(i=seeds,.errorhandling = 'remove')%dopar%optim(f=Lrat,par=i)
  values <- NULL
  for(i in 1:length(solns)){
    values <- c(values,solns[[i]]$value)
  }
  min <- which.min(values)
  sol0 <- solns[[min]]
  x0 <- muM(exp(sol0$par),P)
  S0 <- SIGMAM(exp(sol0$par),P)
  L2 <- c(L2,Statistic(as.matrix(DATA),x,S,x0,S0))
}
L2 <- L2[L2 > 0]
Ln <- length(L2)
printf("Probability of rejecting null2: %0.3f",Fhat(lr2,L2))
printf("%.1f%% points lost from empirical distribution",100-100*Ln/POINTS)

#
# comparison with null3
#

cat("\n************* Inference Under Null 3 *************\n")
x <- mu3
S <- sigma3
clusterExport(cl,c("x","S"))
lr3 <- Statistic(as.matrix(data),x,S,xml,Sml)
printf("LR-statistic for null3: %0.3f",lr3)
L3 <- NULL
while(length(L3) < POINTS){
  
  DATA <- mvrnorm(N,x,S)
  out <- tryCatch(solve( var(DATA) ), error = function(e) e )
  while( is.singular.matrix( var(DATA) ) || any(class(out) == "error") ){
    DATA <- mvrnorm(N,x,S)
    out <- tryCatch(solve( var(DATA) ), error = function(e) e )
  }
  clusterExport(cl,c("DATA"))
  
  clusterEvalQ(cl,{
    Lrat <- function(X){
      X <- exp(X)
      s  <- SIGMAM(X,P)
      m  <- muM(X,P)
      N  <- nrow(DATA)
      mH <- colMeans(DATA)
      sH <- (N-1)*var(DATA)/N
      c  <- N*log(det(s)/det(sH))
      for(i in 1:N){
        c <- c + mahalanobis(DATA[i,],m,s)
      }
      for(i in 1:N){
        c <- c - mahalanobis(DATA[i,],mH,sH)
      }
      return(c)
    }
  })
  solns <- foreach(i=seeds,.errorhandling = 'remove')%dopar%optim(f=Lrat,par=i)
  values <- NULL
  for(i in 1:length(solns)){
    values <- c(values,solns[[i]]$value)
  }
  min <- which.min(values)
  sol0 <- solns[[min]]
  x0 <- muM(exp(sol0$par),P)
  S0 <- SIGMAM(exp(sol0$par),P)
  L3 <- c(L3,Statistic(as.matrix(DATA),x,S,x0,S0))
}
stopCluster(cl)
L3 <- L3[L3 > 0]
Ln <- length(L3)
printf("Probability of rejecting null3: %0.3f",Fhat(lr3,L3))
printf("%.1f%% points lost from empirical distribution",100-100*Ln/POINTS)

#
# comparing null1 to null3
#

cat("\n************* Null1 vs Null3 *************\n")
lr13 <- Statistic(as.matrix(data),mu3,sigma3,mu1,sigma1)
printf("LR-statistic for null1 vs null3: %0.3f",lr13)
printf("Probability of rejecting null3 against null1: %0.3f",Fhat(lr13,L3))

#
# comparing null2 to null3
#

cat("\n************* Null2 vs Null3 *************\n")
lr23 <- Statistic(as.matrix(data),mu3,sigma3,mu2,sigma2)
printf("LR-statistic for null2 vs null3: %0.3f",lr23)
printf("Probability of rejecting null3 against null2: %0.3f",Fhat(lr23,L3))

