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

# with offsets
muM <- function(X,P){
  y <- numeric(2);
  y[1] <- (X[3]*X[4]*P[10] + X[2]*X[3]*P[10] + X[1]*X[4]*(P[11]+P[14]) + X[1]*X[2]*P[14]+X[1]*sqrt(X[2]) )
  y[2] <- (X[3]*X[4]*P[11] + X[1]*X[4]*P[11] + X[2]*X[3]*P[10]+X[3]*sqrt(X[2]) + X[1]*X[2]*P[14]+X[1]*sqrt(X[2]) )
  y <- y/(X[3]*X[4] + X[2]*X[3] + X[1]*X[4])
  return(y)
}

snout <- c(79.3,54.8,57.0,72.8,42.9,58.4,85.8,75.6)
tubes <- c(77.0,50.0,49.1,52.6,41.1,58.6,60.8,68.0)
data <- data.frame(cbind(snout,tubes))
m1 <- 10^-4
m2 <- 10^-4
snout_var <- mean(c(4.6,4.4,5.9,5.8,1.6,6.8,6.5,7.4)^2)
tubes_var <- mean(c(5.9,4.9,4.3,4.7,2.5,7.5,4.6,5.5)^2)
G1 <- snout_var
G2 <- tubes_var
k1 <- 16.97
other_pollinators <- c(27.8,43.3)
ancepsAsolo <- c(31.2,51.7)
k2 <- 4.69
other_flowers <- c(85.3,51.6,46.1,27.5,45.5,69.5,48.8,28.7,42.5,31.1,78.1,41.5,72.5,75.1,62.1,59.4,58.5,43.9,48.9,85.2,46.3,76.8,61.1,43.3,52.8,43.2,89.1,77.1)
th_hat_longir <- mean(other_flowers)+k1
th_hat_anceps <- mean(other_pollinators)+k2
#th_hat_anceps <- mean(ancepsAsolo)+k2
th1 <- th_hat_longir
th2 <- th_hat_anceps
n1 <- 500 # repeat for all permutations of (50,500,5000)
n2 <- 500
n <- nrow(data)
N <- n
P <- c(colMeans(data),(N-1)*var(snout)/N,(N-1)*var(tubes)/N,(N-1)*var(data)[1,2]/N,m1,m2,G1,G2,th1,th2,n1,n2,k1,k2)
TRY <- 500
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
#L1 <- L1[L1 > 0]
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
#L2 <- L2[L2 > 0]
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

