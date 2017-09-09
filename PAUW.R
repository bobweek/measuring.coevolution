library(ggplot2)
source('~/Research/DetectingCoevolution/ErrRates/likelihoodRatio.R')
source('~/Research/DetectingCoevolution/ErrRates/eq.R')

# the matching update functions for traits z1 and z2
updatez1M <- function(z11,z22,B11,A11,G11,n11,th11,m11,mu11){
  z11+rnorm(1,0,sqrt(G11/n11))+G11*(B11*(z22-z11)+A11*(th11-z11))+m11*(mu11-z11)
}
updatez2M <- function(z22,z11,B22,A22,G22,n22,th22,m22,mu22){
  z22+rnorm(1,0,sqrt(G22/n22))+G22*(B22*(z11-z22)+A22*(th22-z22))+m22*(mu22-z22)
}

snout <- c(79.3,54.8,57.0,72.8,42.9,58.4,85.8,75.6)
tubes <- c(77.0,50.0,49.1,52.6,41.1,58.6,60.8,68.0)
data <- data.frame(cbind(snout,tubes))
ggplot(data,aes(snout,tubes))+geom_point()+stat_smooth(method='lm',se=FALSE)+theme_bw()

other_pollinators <- c(27.8,43.3)
other_flowers <- c(85.3,51.6,46.1,27.5,45.5,69.5,48.8,28.7,42.5,31.1,78.1,41.5,72.5,75.1,62.1,59.4,58.5,43.9,48.9,85.2,46.3,76.8,61.1,43.3,52.8,43.2,89.1,77.1)
th_hat_anceps <- mean(other_pollinators)
th_hat_longir <- mean(other_flowers)

# values of anceps at locations without longirostris
ancepsAsolo <- c(31.2,51.7)
mean(ancepsAsolo)
th_hat_anceps <- mean(c(th_hat_anceps,mean(ancepsAsolo)))

# have determined estimates for abiotic optima, can use average sample variation as additive genetic variance
# still need estimates of population sizes and migration rates

# could start by assuming zero migration and making up estimates for local population size

# the methods section of the paper suggests that these populations are in the hundreds.
# considering how weakly dependent this method is on the value of population size, we can safely assume n=300

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
N <- nrow(data)
PARS <- c(colMeans(data),N*var(snout)/(N-1),N*var(tubes)/(N-1),N*var(data)[1,2]/(N-1),m1,m2,G1,G2,th1,th2,n1,n2);


muM <- function(X,P){
  y <- numeric(2);
  y[1] <- (X[4]*X[1]*P[11] + X[3]*(X[4] + X[2])*P[10])/(X[4]*X[1] + X[3]*(X[4] + X[2]));
  y[2] <- (X[3]*X[2]*P[10] + X[4]*(X[3] + X[1])*P[11])/(X[3]*X[2] + X[4]*(X[3] + X[1]));
  y
}

SIGMAM <- function(X,P){
  y <- matrix(0,2,2);
  y[1,1] <- (P[8]*((P[8]*P[9]*(X[1])^2)/P[13] + ((X[4]*P[9] + P[7])*(X[4]*P[9] + P[6] + P[7] + P[8]*(X[3] + X[1])) + P[9]*(X[3]*P[8] + 2*X[4]*P[9] + P[6] + 2*P[7])*X[2] + 
                                                   (P[9])^2*(X[2])^2)/P[12]))/
    (2.*((X[4]*P[9] + P[7])*(P[6] + P[8]*(X[3] + X[1])) + P[9]*(X[3]*P[8] + P[6])*X[2])*(P[6] + P[7] + P[8]*(X[3] + X[1]) + P[9]*(X[4] + X[2])));
  y[2,2] <- (P[9]*(P[12]*(P[6] + P[8]*(X[3] + X[1]))*(X[4]*P[9] + P[6] + P[7] + P[8]*(X[3] + X[1])) + P[9]*(X[3]*P[8] + P[6])*P[12]*X[2] + P[8]*P[9]*P[13]*(X[2])^2))/
    (2.*P[12]*P[13]*((X[4]*P[9] + P[7])*(P[6] + P[8]*(X[3] + X[1])) + P[9]*(X[3]*P[8] + P[6])*X[2])*(P[6] + P[7] + P[8]*(X[3] + X[1]) + P[9]*(X[4] + X[2])));
  y[1,2] <- (P[8]*P[9]*(P[12]*X[1]*(P[6] + P[8]*(X[3] + X[1])) + P[13]*X[2]*(P[7] + P[9]*(X[4] + X[2]))))/
    (2.*P[12]*P[13]*((X[4]*P[9] + P[7])*(P[6] + P[8]*(X[3] + X[1])) + P[9]*(X[3]*P[8] + P[6])*X[2])*(P[6] + P[7] + P[8]*(X[3] + X[1]) + P[9]*(X[4] + X[2])));
  y[2,1] <- y[1,2];
  y
}

snout_var <- mean(c(4.6,4.4,5.9,5.8,1.6,6.8,6.5,7.4)^2)
tubes_var <- mean(c(5.9,4.9,4.3,4.7,2.5,7.5,4.6,5.5)^2)
G1 <- snout_var
G2 <- tubes_var
th1 <- th_hat_longir+25
th2 <- th_hat_anceps+15

N <- nrow(data)
Lrat <- function(X){
  
  N <- nrow(data)
  
  mH <- colMeans(data)
  sH <- N*var(data)/(N-1)
  
  MOO <- muM(X,PARS)
  SEEG <- SIGMAM(X,PARS)
  
  c <- N*log(det(SEEG)/det(sH))
  
  for(i in 1:N){
    c <- c + as.numeric(data[i,] - MOO)%*%solve(SEEG)%*%as.numeric(data[i,] - MOO)
  }
  for(i in 1:N){
    c <- c - as.numeric(data[i,] - mH)%*%solve(sH)%*%as.numeric(data[i,] - mH)
  }
  return(c)
  
}

m1 <- 10^-4
m2 <- 10^-4
n1 <- 100
n2 <- 100
PARS <- c(colMeans(data),var(snout),var(tubes),var(data)[1,2],m1,m2,G1,G2,th1,th2,n1,n2);

# why not try minimizing wrt to m's and n's as well???

# N <- 1000
# fdata <- mvrnorm(N,colMeans(data),var(data))
sol <- optim(f=Lrat, par = runif(4,0,1e-4))
estimates <- sol$par
trials <- 0;
while( trials<1 ){
  trial <- runif(4,0,1e-4);
  s <- optim(f=Lrat, par = trial)
  if( s$value<sol$value ){
    estimates <- s$par;  
    sol <- s;
  } 
  trials <- trials+1;
}

print("********* USING RATIO AS METRIC *********")

print("Estimated Selection Strengths")
print(estimates)

print(muM(estimates,PARS))
print(SIGMAM(estimates,PARS))
print(colMeans(data))
print(var(data))

print("Likelihood Ratio")
print(Lrat(estimates))
