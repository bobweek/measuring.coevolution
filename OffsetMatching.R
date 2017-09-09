library("MASS")
library("matrixcalc")

#         1 2  3  4  5  6
# P <- c(mu,V,G1,G2,n1,n2) DATA + BACKGROUND PARAMETERS
#         1  2  3  4
# X <- c(M1,M2,D1,D2)

mu <- function(M1,M2,D1,D2,G1,G2){
  y <- ( D1*G1-D2*G2 )/( G1*M1+G2*M2 )
  return(y)
}

sigma <- function(M1,M2,G1,G2,n1,n2){
  y <- ( G1/n1 + G2/n2 )/(2*( G1*M1+G2*M2 ))
  return(y)
}

# has a great null but does not have a unique solution!
Lrat <- function(X){
  N  <- length(data)
  Mh <- mean(data)
  Sh <- (N-1)*var(data)/N
  M  <-    mu(X[1],X[2],X[3],X[4],P[3],P[4])
  S  <- sigma(X[1],X[2],P[4],P[5],P[5],P[6])
  c  <- N*log(S/Sh)
  for(i in 1:N){
    c <- c + (data[i] - M)^2/S
  }
  for(i in 1:N){
    c <- c - (data[i] - Mh)^2/Sh
  }
  return(c)
}
optim(f=Lrat, par = runif(4,0,1e-4))

# has a unique solution but no decent null!
Lrat <- function(X){
  N  <- length(data)
  Mh <- mean(data)
  Sh <- (N-1)*var(data)/N
  M  <-    mu(X[1],X[1],X[2],X[2],P[3],P[4])
  S  <- sigma(X[1],X[1],P[4],P[5],P[5],P[6])
  c  <- N*log(S/Sh)
  for(i in 1:N){
    c <- c + (data[i] - M)^2/S
  }
  for(i in 1:N){
    c <- c - (data[i] - Mh)^2/Sh
  }
  return(c)
}
optim(f=Lrat, par = runif(2,0,1e-4))

# a decent null but no unique solution!
Lrat <- function(X){
  N  <- length(data)
  Mh <- mean(data)
  Sh <- (N-1)*var(data)/N
  M  <-    mu(X[1],X[2],0,0,P[3],P[4])
  S  <- sigma(X[1],X[2],P[4],P[5],P[5],P[6])
  c  <- N*log(S/Sh)
  for(i in 1:N){
    c <- c + (data[i] - M)^2/S
  }
  for(i in 1:N){
    c <- c - (data[i] - Mh)^2/Sh
  }
  return(c)
}
optim(f=Lrat, par = runif(2,0,1e-4))


other_pollinators <- c(27.8,43.3)
other_flowers <- c(85.3,51.6,46.1,27.5,45.5,69.5,48.8,28.7,42.5,31.1,78.1,41.5,72.5,75.1,62.1,59.4,58.5,43.9,48.9,85.2,46.3,76.8,61.1,43.3,52.8,43.2,89.1,77.1)
th_hat_anceps <- mean(other_pollinators)
th_hat_longir <- mean(other_flowers)
ancepsAsolo <- c(31.2,51.7)
mean(ancepsAsolo)
th_hat_anceps <- mean(c(th_hat_anceps,mean(ancepsAsolo)))
snout <- c(79.3,54.8,57.0,72.8,42.9,58.4,85.8,75.6)
tubes <- c(77.0,50.0,49.1,52.6,41.1,58.6,60.8,68.0)
data <- snout-tubes
m1 <- 10^-4
m2 <- 10^-4
snout_var <- mean(c(4.6,4.4,5.9,5.8,1.6,6.8,6.5,7.4)^2)
tubes_var <- mean(c(5.9,4.9,4.3,4.7,2.5,7.5,4.6,5.5)^2)
G1 <- snout_var
G2 <- tubes_var
th1 <- th_hat_longir+25
th2 <- th_hat_anceps+15
n1 <- 100
n2 <- 100
N <- length(data)
#P <- c(colMeans(data),n*var(snout)/(n-1),n*var(tubes)/(n-1),n*var(data)[1,2]/(n-1),m1,m2,G1,G2,th1,th2,n1,n2)
P <- c(mean(data),(N-1)*var(data)/N,G1,G2,n1,n2)

# find ml parameters for null1: B1=0
# Y1 <- numeric(4);
# Y1[2] <- ( P[5]*P[9]*P[12]*P[3] + P[5]*P[8]*P[13]*P[4] ) / ( 2*P[12]*P[13]*P[3]*P[9]*( P[3]*P[4]-P[5]*P[5] ) )
# Y1[3] <- ( 2*P[12]*P[3] )^(-1) - (P[6]/P[8])
# Y1[4] <- ( P[5]*P[5]*P[8]*P[13] - P[5]*P[9]*P[12]*P[3] + P[9]*P[12]*P[3]*P[3] - P[5]*P[8]*P[13]*P[4]  ) / ( 2*P[12]*P[13]*P[3]*P[9]*( P[3]*P[4]-P[5]*P[5] ) ) - (P[7]/P[9])
# mu1 <- muM(Y1,P)
# sigma1 <- SIGMAM(Y1,P)
Lrat <- function(X){
  Y <- c(0,X)
  N <- nrow(data)
  mH <- colMeans(data)
  sH <- (N-1)*var(data)/N
  MOO <- muM(Y,P)
  SEEG <- SIGMAM(Y,P)
  c <- N*log(det(SEEG)/det(sH))
  for(i in 1:N){
    c <- c + as.numeric(data[i,] - MOO)%*%solve(SEEG)%*%as.numeric(data[i,] - MOO)
  }
  for(i in 1:N){
    c <- c - as.numeric(data[i,] - mH)%*%solve(sH)%*%as.numeric(data[i,] - mH)
  }
  return(c)
}
sol <- optim(f=Lrat, par = runif(3,0,1e-4))
Y <- sol$par
Y <- c(0,Y)
mu1 <- muM(Y,P)
sigma1 <- SIGMAM(Y,P)

# find ml parameters for null2: B2=0
# Y2 <- numeric(4)
# Y2[1] <- ( P[5]*P[8]*P[13]*P[4] + P[5]*P[9]*P[12]*P[3] ) / ( 2*P[12]*P[13]*P[4]*P[8]*( P[3]*P[4]-P[5]*P[5] ) )
# Y2[3] <- ( P[5]*P[5]*P[9]*P[12] - P[5]*P[8]*P[13]*P[4] + P[8]*P[13]*P[4]*P[4] - P[5]*P[9]*P[12]*P[3]  ) / ( 2*P[13]*P[12]*P[4]*P[8]*( P[3]*P[4]-P[5]*P[5] ) ) - (P[6]/P[8])
# Y2[4] <- ( 2*P[13]*P[4] )^(-1) - (P[7]/P[9])
# mu2 <- muM(Y2,P)
# sigma2 <- SIGMAM(Y2,P)
Lrat <- function(X){
  Y <- c(X[1],0,X[2],X[3])
  N <- nrow(data)
  mH <- colMeans(data)
  sH <- (N-1)*var(data)/N
  MOO <- muM(Y,P)
  SEEG <- SIGMAM(Y,P)
  c <- N*log(det(SEEG)/det(sH))
  for(i in 1:N){
    c <- c + as.numeric(data[i,] - MOO)%*%solve(SEEG)%*%as.numeric(data[i,] - MOO)
  }
  for(i in 1:N){
    c <- c - as.numeric(data[i,] - mH)%*%solve(sH)%*%as.numeric(data[i,] - mH)
  }
  return(c)
}
sol <- optim(f=Lrat, par = runif(3,0,1e-4))
Y <- sol$par
Y <- c(Y[1],0,Y[2],Y[3])
mu2 <- muM(Y,P)
sigma2 <- SIGMAM(Y,P)

# find ml parameters for null3: B1=B2=0
# Y3 <- numeric(4)
# Y3[3] <- ( 2*P[12]*P[3] )^(-1) - (P[6]/P[8])
# Y3[4] <- ( 2*P[13]*P[4] )^(-1) - (P[7]/P[9])
# mu3 <- muM(Y3,P)
# sigma3 <- SIGMAM(Y3,P)
Lrat <- function(X){
  Y <- c(0,0,X)
  N <- nrow(data)
  mH <- colMeans(data)
  sH <- (N-1)*var(data)/N
  MOO <- muM(Y,P)
  SEEG <- SIGMAM(Y,P)
  c <- N*log(det(SEEG)/det(sH))
  for(i in 1:N){
    c <- c + as.numeric(data[i,] - MOO)%*%solve(SEEG)%*%as.numeric(data[i,] - MOO)
  }
  for(i in 1:N){
    c <- c - as.numeric(data[i,] - mH)%*%solve(sH)%*%as.numeric(data[i,] - mH)
  }
  return(c)
}
sol <- optim(f=Lrat, par = runif(2,0,1e-4))
Y <- sol$par
Y <- c(0,0,Y)
mu3 <- muM(Y,P)
sigma3 <- SIGMAM(Y,P)

Statistic <- function(D,mu_0,V_0,mu_1,V_1){
  mu_t <- colMeans(D)
  n <- nrow(D)
  c <- n*log(det(V_0)/det(V_1))
  c0 <- 0
  for( i in 1:n){
    c0 <- c0 + t(D[i,]-mu_0)%*%solve(V_0)%*%(D[i,]-mu_0)
  }
  c1 <- 0
  for( i in 1:n){
    c1 <- c1 + t(D[i,]-mu_1)%*%solve(V_1)%*%(D[i,]-mu_1)
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



# FAKE SETUP
#muM(P[6],P[7],P[8],P[9],P[10],P[11],P[12],P[13],1e-5,1e-5,1e-6,1e-6)
#SIGMAM(P[6],P[7],P[8],P[9],P[10],P[11],P[12],P[13],1e-5,1e-5,1e-6,1e-6)
#data <- mvrnorm(n,x,S)

#xml <- colMeans(data)
#Sml <- (N-1)*var(data)/N

POINTSS <- 1000
N <- nrow(data)
x <- mu2
S <- sigma2
L <- NULL
for(i in 1:POINTSS){
  DATA <- mvrnorm(n,x,S)
  #x0 <- colMeans(DATA)
  #S0 <- (N-1)*var(DATA)/N
  
  Lrat <- function(X){
    N <- nrow(DATA)
    mH <- colMeans(DATA)
    sH <- (N-1)*var(DATA)/N
    MOO <- muM(X,P)
    SEEG <- SIGMAM(X,P)
    c <- N*log(det(SEEG)/det(sH))
    for(i in 1:N){
      c <- c + as.numeric(DATA[i,] - MOO)%*%solve(SEEG)%*%as.numeric(DATA[i,] - MOO)
    }
    for(i in 1:N){
      c <- c - as.numeric(DATA[i,] - mH)%*%solve(sH)%*%as.numeric(DATA[i,] - mH)
    }
    return(c)
  }
  sol <- optim(f=Lrat, par = runif(4,0,1e-4))
  x0  <- muM(sol$par,P)
  S0  <- SIGMAM(sol$par,P)
  
  L <- c(L,Statistic(as.matrix(DATA),x,S,x0,S0))
}

L <- L[L < 1e+02] # get rid of everything ridiculously large
L <- L[L > 0]     # get rid of everything ridiculously small

z <- seq(min(L),max(L),length.out=100)
Fh <- NULL
for(i in 1:100){
  Fh <- c(Fh,Fhat(z[i],L))
}
plot(z,Fh)

Lrat <- function(X){
  N <- nrow(data)
  mH <- colMeans(data)
  sH <- (N-1)*var(data)/N
  MOO <- muM(X,P)
  SEEG <- SIGMAM(X,P)
  c <- N*log(det(SEEG)/det(sH))
  for(i in 1:N){
    c <- c + as.numeric(data[i,] - MOO)%*%solve(SEEG)%*%as.numeric(data[i,] - MOO)
  }
  for(i in 1:N){
    c <- c - as.numeric(data[i,] - mH)%*%solve(sH)%*%as.numeric(data[i,] - mH)
  }
  return(c)
}
sol <- optim(f=Lrat, par = runif(4,0,1e-4))
xml <- muM(sol$par,P)
Sml <- SIGMAM(sol$par,P)

lr <- Statistic(as.matrix(data),x,S,xml,Sml)

print(sol$par)
print(Fhat(lr,L))
print(lr)