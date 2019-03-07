OR <- function(thing){
  thing <- as.vector(thing)
  ln <- length(thing)
  bool <- thing[1]||thing[2]
  if(ln>2){
    for(i in 3:ln) bool <- bool||thing[i]
  }
  return(bool)
}

AND <- function(thing){
  thing <- as.vector(thing)
  ln <- length(thing)
  bool <- thing[1]&&thing[2]
  if(ln>2){
    for(i in 3:ln) bool <- bool&&thing[i]
  }
  return(bool)
}

dist_fct <- function(value,stuff){
  p <- 0
  for(s in stuff){
    if(value > s) p <- p+1
  }
  return(p/length(stuff))
}

q_fct <- function(p,stuff){
  ordered_stuff <- stuff[order(stuff)]
  h <- floor(p*length(stuff))
  return(ordered_stuff[h])
}

likelihood <- function(d,m,S){
  N <- dim(d)[1]
  val <- 0
  for(i in 1:N) val <- val + dmvnorm(d[i,],m,S,log=T)
  return(val)
}

ml_sol <- function(mu1,mu2,V1,V2,C,n1,n2,th1,th2,G1,G2){
  
  A1 <- (-(G2*(n1)^2*(C - V1)*V1*((C)^2*(th1 - mu1) + (V1)^2*(-th2 + mu2) + C*V1*(-th1 + th2 - mu1 + mu2))) + 
           G1*n1*n2*V1*(V1*(V2)^2*(th1 - 2*th2 + 3*mu1 - 2*mu2) + (C)^3*(-th2 + mu2) + (C)^2*(V1*(-th2 + mu2) + V2*(th1 + th2 - 3*mu1 + mu2)) - C*V2*(V2*(th1 - mu1) + V1*(th1 - 3*th2 + mu1 + mu2))) + 
           C*sqrt((n1)^2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)^2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                            (G1*n2*V1*(V2)^2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)^3*(-th2 + mu2) + 
                               (C)^2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))^2)) - 
           V1*sqrt((n1)^2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)^2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                             (G1*n2*V1*(V2)^2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)^3*(-th2 + mu2) + 
                                (C)^2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))^2)))/
    (4.*G1*(n1)^2*n2*V1*(-(C)^2 + V1*V2)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)))
  
  A2 <- (G1*n1*n2*V1*(C - V2)*((V2)^2*(-th1 + mu1) + (C)^2*(th2 - mu2) + C*V2*(th1 - th2 + mu1 - mu2)) + 
           G2*(n1)^2*V1*((C)^3*(th1 - mu1) - (C)^2*(V2*(-th1 + mu1) + V1*(th1 + th2 + mu1 - 3*mu2)) + (V1)^2*V2*(2*th1 - th2 + 2*mu1 - 3*mu2) + C*V1*(V1*(th2 - mu2) + V2*(-3*th1 + th2 + mu1 + mu2))) - 
           C*sqrt((n1)^2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)^2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                            (G1*n2*V1*(V2)^2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)^3*(-th2 + mu2) + 
                               (C)^2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))^2)) + 
           V2*sqrt((n1)^2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)^2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                             (G1*n2*V1*(V2)^2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)^3*(-th2 + mu2) + 
                                (C)^2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))^2)))/
    (4.*G2*(n1)^2*n2*V1*(-(C)^2 + V1*V2)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)))
  
  B1 <- (G2*(n1)^2*V1*((C)^2*(th1 - mu1) + (V1)^2*(-th2 + mu2) + C*V1*(-th1 + th2 - mu1 + mu2)) + 
           G1*n1*n2*(V1*(V2)^2*(-th1 + mu1) + C*V1*V2*(-th1 + th2 - mu1 + mu2) + (C)^2*(2*V2*(th1 - mu1) + V1*(-th2 + mu2))) - 
           sqrt((n1)^2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)^2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                          (G1*n2*V1*(V2)^2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)^3*(-th2 + mu2) + 
                             (C)^2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))^2)))/
    (4.*G1*(n1)^2*n2*((C)^2 - V1*V2)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)))
  
  B2 <- (G1*n1*n2*V1*V2*((V2)^2*(-th1 + mu1) + (C)^2*(th2 - mu2) + C*V2*(th1 - th2 + mu1 - mu2)) + 
           G2*(n1)^2*V1*((C)^2*(V2*(-th1 + mu1) + 2*V1*(th2 - mu2)) + C*V1*V2*(th1 - th2 + mu1 - mu2) + (V1)^2*V2*(-th2 + mu2)) - 
           V2*sqrt((n1)^2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)^2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                             (G1*n2*V1*(V2)^2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)^3*(-th2 + mu2) + 
                                (C)^2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))^2)))/
    (4.*G2*(n1)^2*n2*V1*(-(C)^2 + V1*V2)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)))
  
  k <- (-(G1*n1*n2*V1*(C - V2)*(V2*(-th1 + mu1) + C*(th2 - mu2))) - G2*(n1)^2*(C - V1)*V1*(C*(th1 - mu1) + V1*(-th2 + mu2)) + 
          sqrt((n1)^2*(4*C*G1*n2*(G2*n1*V1 + G1*n2*V2)*((C)^2 - V1*V2)*(th1 - mu1)*(C*(V2*(-th1 + mu1) + V1*(th2 - mu2)) + V1*V2*(th1 - th2 + mu1 - mu2)) + 
                         (G1*n2*V1*(V2)^2*(-th1 + mu1) - C*V1*(G2*n1*V1 + G1*n2*V2)*(th1 - th2 + mu1 - mu2) + G2*n1*(V1)^3*(-th2 + mu2) + 
                            (C)^2*(G2*n1*V1*(th1 - mu1) + G1*n2*(2*V2*th1 - V1*th2 - 2*V2*mu1 + V1*mu2)))^2)))/(2.*C*n1*V1*(G2*n1*V1 + G1*n2*V2))
  
  names(A1) <- NULL
  names(A2) <- NULL
  names(B1) <- NULL
  names(B2) <- NULL
  names(k) <- NULL
  ML_sol <- list(A1,A2,B1,B2,k)
  names(ML_sol) <- c("A1","A2","B1","B2","k")
  return(ML_sol)
  
}


ml_sol_k0 <- function(mu1,mu2,V1,V2,C,n1,n2,th1,th2,G1,G2){
  
  A1 <- (mu2-mu1)/(2*n1*(V1*(mu2-th1)-C*(mu1-th1)))
  
  A2 <- (mu1-mu2)/(2*n2*(V2*(mu1-th2)-C*(mu2-th2)))
  
  B1 <- (mu1-th1)/(2*n1*(V1*(mu2-th1)-C*(mu1-th1)))
  
  B2 <- (mu2-th2)/(2*n2*(V2*(mu1-th2)-C*(mu2-th2)))
  
  names(A1) <- NULL
  names(A2) <- NULL
  names(B1) <- NULL
  names(B2) <- NULL
  ML_sol <- list(A1,A2,B1,B2,0)
  names(ML_sol) <- c("A1","A2","B1","B2","k")
  return(ML_sol)
  
}


mu <- function(A1,A2,B1,B2,k,n1,n2,th1,th2){
  
  mu1 <- (A1*(A2 + B2)*th1 + B1*(2*B2*k + A2*(th2 + k)))/(A2*B1 + A1*(A2 + B2))
  mu2 <- (A1*B2*th1 + A1*A2*th2 + A2*B1*th2 + (A1 + 2*B1)*B2*k)/(A2*B1 + A1*(A2 + B2))
  return(c(mu1,mu2))
  
}

SIGMA <- function(A1,A2,B1,B2,k,n1,n2,th1,th2,G1,G2){
  
  COV <- (B1*(A1 + B1)*G1*n1 + B2*(A2 + B2)*G2*n2)/
    (2.*(A2*B1 + A1*(A2 + B2))*((A1 + B1)*G1 + (A2 + B2)*G2)*n1*n2)
  VAR1 <- ((B1)^2*G1*n1 + A2*B1*G1*n2 + (A2 + B2)*(A1*G1 + (A2 + B2)*G2)*n2)/
    (2.*(A2*B1 + A1*(A2 + B2))*((A1 + B1)*G1 + (A2 + B2)*G2)*n1*n2)
  VAR2 <- (((A1 + B1)^2*G1 + (A2*B1 + A1*(A2 + B2))*G2)*n1 + (B2)^2*G2*n2)/
    (2.*(A2*B1 + A1*(A2 + B2))*((A1 + B1)*G1 + (A2 + B2)*G2)*n1*n2)
  
  y <- matrix(c(VAR1,COV,COV,VAR2),2,2)
  return(y)
  
}

Fhat <- function(x,L){
  N <- length(L)
  L0 <- NULL
  for(i in 1:N){
    if(L[i]<=x) L0 <- c(L0,L[i])
  }
  p <- length(L0)/N
  return(p)
}

SIGMA_gf <- function(A1,A2,B1,B2,k,n1,n2,th1,th2,G1,G2,m1,m2){
  
  COV <- (G1*G2*(B1*((A1 + B1)*G1 + m1)*n1 + B2*((A2 + B2)*G2 + m2)*n2))/(2.*((A1 + B1)*G1 + (A2 + B2)*G2 + m1 + m2)*(G2*(A2*B1*G1 + A1*(A2 + B2)*G1 + (A2 + B2)*m1) + ((A1 + B1)*G1 + m1)*m2)*n1*n2)
  VAR1 <- (G1*((B1)^2*G1*G2*n1 + B1*G1*(A2*G2 + m2)*n2 + ((A2 + B2)*G2 + m2)*(A1*G1 + (A2 + B2)*G2 + m1 + m2)*n2))/
    (2.*((A1 + B1)*G1 + (A2 + B2)*G2 + m1 + m2)*(G2*(A2*B1*G1 + A1*(A2 + B2)*G1 + (A2 + B2)*m1) + ((A1 + B1)*G1 + m1)*m2)*n1*n2)
  VAR2 <- (G2*(((A1)^2*(G1)^2 + (B1)^2*(G1)^2 + m1*((A2 + B2)*G2 + m1 + m2) + B1*G1*(A2*G2 + 2*m1 + m2) + A1*G1*(2*B1*G1 + (A2 + B2)*G2 + 2*m1 + m2))*n1 + (B2)^2*G1*G2*n2))/
    (2.*((A1 + B1)*G1 + (A2 + B2)*G2 + m1 + m2)*(G2*(A2*B1*G1 + A1*(A2 + B2)*G1 + (A2 + B2)*m1) + ((A1 + B1)*G1 + m1)*m2)*n1*n2)
  
  y <- matrix(c(VAR1,COV,COV,VAR2),2,2)
  return(y)
  
}
