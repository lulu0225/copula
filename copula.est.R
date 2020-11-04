library(quadprog)
library(Matrix)
library(graphics)
library(copula)  

# X: a n*2 matrix of data
# m1, m2: degrees of Bernstein polynomial estimated by grid search if is m1 = NULL and m2 = NULL 
# is.pqd: is.pqd = F for unconstrained copula; is.pqd = T for PQD-constrained copula

copula.est <- function(X, m1 = NULL, m2 = NULL, is.pqd = F){

# Cm = Gm^T theta + XX + YY
Gm.pqd <- function (u1, u2, m1, m2) {
  G0 <- matrix(0, m1 ,m2)
  for (k1 in 1:m1){
    for (k2 in 1:m2){
      G0[k1, k2] <- choose((m1+1),k1) * choose((m2+1),k2) * u1^k1 * (1 - u1)^(m1 - k1 + 1) * u2^k2 * (1 - u2)^(m2 - k2 + 1)
    }
  }
  as.vector(t(G0))
}


Cm.pqd <- function(u1, u2, m1, m2, theta){
  
  XX <- 0
  for (k1 in 1 : (m1+1)){
    XX <- XX + k1/(m1+1) * choose((m1+1),k1) * u1^k1 * (1 - u1)^(m1 - k1 + 1) * u2^(m2+1)
  }
  
  YY <- 0
  for (k2 in 1 : m2){
    YY <- YY + k2/(m2+1) * choose((m2+1),k2) * u2^k2 * (1 - u2)^(m2 - k2 + 1) * u1^(m1+1)
  }

  pro <- crossprod(Gm.pqd(u1, u2, m1, m2), theta) + XX + YY
  pro
}


opt.pqd <- function (m1, m2, U1, U2){
  
  # constraints A^T theta >= b0
  
  # PQD
  W2 <- matrix(0, m1, m2)
  
  for (k1 in 1 : m1){
    for (k2 in 1 : m2){
      
      W2[k1, k2] <- k1/(m1 + 1) * k2/(m2 + 1)
      
    }
  }
  
  w2 <- as.vector(t(W2))
  
  cons.theta <- function(x, y, m1, m2){
    Q <- matrix(0, m1, m2)
    Q[x, y] <- 1
    Q[x + 1, y] <- -1
    Q[x, y + 1] <- -1
    Q[x + 1, y + 1] <- 1
    Q
  }
  
  cons.theta.1 <- function(x,y,m1,m2){
    Q <- matrix(0, m1, m2)
    Q[x, y] <- 1
    Q[x + 1, y] <- -1
    Q
  }
  
  cons.theta.2 <- function(x,y,m1,m2){
    Q <- matrix(0, m1, m2)
    Q[x, y] <- 1
    Q[x, y + 1] <- -1
    Q
  }

  cons.theta.3 <- function(x,y,m1,m2){
    Q <- matrix(0, m1, m2)
    Q[x, y] <- -1
    Q[x + 1, y] <- 1
    Q
  }
  
  cons.theta.4 <- function(x,y,m1,m2){
    Q <- matrix(0, m1, m2)
    Q[x, y] <- -1
    Q[x, y + 1] <- 1
    Q
  }
  
  R <- c()
  for (k1 in 1: (m1 - 1)){
    for (k2 in 1 : (m2 - 1)){
      r <- as.vector(t(cons.theta(k1, k2, m1, m2)))
      R <- rbind(R, r)
    }
  }
  
  R1 <- c()
  for (k1 in 1: (m1 - 1)){
    r1 <- as.vector(t(cons.theta.1(k1, m2, m1, m2)))
    R1 <- rbind(R1, r1)
  }
  
  R2 <- c()
  for (k2 in 1: (m2 - 1)){
    r2 <- as.vector(t(cons.theta.2(m1, k2, m1, m2)))
    R2 <- rbind(R2, r2)
  }
  
  R3 <- c()
  for (k1 in 1: (m1 - 1)){
    r3 <- as.vector(t(cons.theta.3(k1, 1, m1, m2)))
    R3 <- rbind(R3, r3)
  }
  
  R4 <- c()
  for (k2 in 1: (m2 - 1)){
    r4 <- as.vector(t(cons.theta.4(1, k2, m1, m2)))
    R4 <- rbind(R4, r4)
  }
  
  r1 <- c(1, rep(0, m1*m2 - 1))
  r2 <- - c(rep(0,m2 - 1), -1, rep(0,(m1 - 1)*m2))
  r3 <- - c(rep(0,(m1 - 1)*m2), -1, rep(0,m2 - 1))
  r4 <- c(rep(0, m1*m2 - 1), 1)
  
  bb0 <- rep(0, (m1-1)*(m2-1))
  bb1 <- rep(-1/(m1+1), m1-1)
  bb2 <- rep(-1/(m2+1), m2-1)
  bb3 <- rep(0, m1-1)
  bb4 <- rep(0, m2-1)
  
  b1 <- 0
  b2 <- - 1/(m1 + 1)
  b3 <- - 1/(m2 + 1)
  b4 <- m1/(m1+1) + m2/(m2+1) -1
  
  if (is.pqd == TRUE){
  #PQD
  b0 <- c(w2, bb0, bb1, bb2, bb3, bb4, b1, b2, b3, b4)
  A <- cbind(diag((m1)*(m2)), t(R), t(R1), t(R2), t(R3), t(R4), r1, r2, r3, r4)
  }
  else{
  #non-PQD  
  b0 <- c(bb0, bb1, bb2, bb3, bb4, b1, b2, b3, b4)
  A <- cbind(t(R), t(R1), t(R2), t(R3), t(R4), r1, r2, r3, r4)
  }

  # get estimate of theta
  n <- length(U1)
  
  # empirical copula 
  Cn <- c()
  
  # theta minimize (1/2 theta^T D theta - d^T theta )
  D <- matrix (0, (m1 * m2), (m1 * m2))
  d <- c(rep(0, (m1 * m2)))
  
  XXX <- function(u1, u2, m1, m2){
    XX <- 0
    for (k1 in 1 : (m1+1)){
      XX <- XX + k1/(m1+1) * choose((m1+1),k1) * u1^k1 * (1 - u1)^(m1 - k1 + 1) * u2^(m2+1)
    }
    XX
  }
  
  YYY <- function(u1, u2, m1, m2){
    YY <- 0
    for (k2 in 1 : m2){
      YY <- YY + k2/(m2+1) * choose((m2+1),k2) * u2^k2 * (1 - u2)^(m2 - k2 + 1) * u1^(m1+1)
    }
    YY
  }
  
  for ( i in 1:n ){
    Cn[i] <- 1 / (n + 1) * sum(U1 <= U1[i] & U2 <= U2[i])
    D <- D + 2 * n^2 / (n+1)^2 / (n/(n+1)^2 * Cn[i] * (1 - Cn[i]) + .001) * tcrossprod(Gm.pqd(U1[i], U2[i], m1, m2))
    d <- d +  2 * ( n/(n + 1)*Cn[i] - n^2 / (n+1)^2 *(XXX(U1[i], U2[i], m1, m2) + YYY(U1[i], U2[i], m1, m2)) )/ (n/(n+1)^2 * Cn[i] * (1 - Cn[i]) + .001) * Gm.pqd(U1[i], U2[i], m1, m2)
    
  }
  
  D <- D/n
  d <- d/n
  
  # Compute the nearest positive definite matrix
  Dn <- nearPD(D)$mat
  
  sc <- norm(Dn, "2")
  
  D.test <- Dn / sc
  d.test <- d / sc
  
  theta <- solve.QP(Dmat = D.test,  dvec = d.test, Amat = A, bvec = b0, meq = 0, factorized=FALSE)$solution
  
  theta
}


# m1 and m2 "minimize" Dm
Dm <- function (m1, m2, U1, U2, x.points, y.points){
  theta <- opt.pqd(m1, m2, U1, U2)
  
  N <-length(x.points)
  Rn <- 0
  for (j in 1:N){
    Cn <- c()
    Cest <- c()
    for (i in 1:N){
      Cn[i] <- 1 / (n+1) * sum(U1 <= x.points[i] & U2 <= y.points[j])
      Cest[i] <- Cm.pqd(x.points[i], y.points[j], m1, m2, theta)
       Rn <- Rn + (Cn[i]-Cest[i])^2/(Cest[i]*(1-Cest[i]) + .001)
    }
  }
  Rn / N^2
  
}


# find optimal m1 and m2
find <- function(U1, U2, tol =  10^(-2)){
  
  x.points <- seq(0.001, 1, length.out = 100)
  y.points <- seq(0.001, 1, length.out = 100)
  
  n <- length(U1)
  a = seq(0.3,0.7,by=0.05)
  m = ceiling(n^a) - 1
  
  minval = 1000
  minline = 1000
  for (k in 1:9){
    i = k
    j = 1
    value <- c()
    while( i >= 1 ){
      value[j] <- Dm(m[i], m[j], U1, U2, x.points, y.points)
      
      if(value[j] < minline){
        minline = value[j]
        locline <- c(i, j)
      }
      
      i = i - 1
      j = j + 1
    }
    
    
    if ( k > 1 & (minval - minline)/minval < tol ) {
      break
    }
    
    minval <- minline
    minloc <- locline
  }
  
  
  if(sum(minloc) == 10){
    
    for (k in 2: 9){
      i = 9
      j = k
      value <- c()
      while(j <= 9){
        value[10 - i] <- Dm(m[i], m[j], U1, U2, x.points, y.points)
        
        if(value[10 - i] < minline){
          minline = value[10 - i]
          locline <- c(i, j)
        }
        
        i = i - 1
        j = j + 1
      }
      
      if ( k > 1 & (minval - minline)/minval < tol ) {
        break
      }
      
      minval <- minline
      minloc <- locline
    }
    
  }
  
  m1 <- m[minloc[1]]
  m2 <- m[minloc[2]]
  
  return(list('minval' = minval, 'm1' = m1, 'm2' = m2))
}

#########################

U1 <- X[,1]
U2 <- X[,2]

# transform to pseudo-observations
ec1 <- ecdf(U1)
ec2 <- ecdf(U2)
U1 <- n / (n+1) * ec1(U1)
U2 <- n / (n+1) * ec2(U2)

# grid seach for m1 and m2 if not imput
if (is.null(m1) == T & is.null(m2) == T){
  
m <- find(U1, U2, tol =  10^(-2))
m1 <- m$m1 
m2 <- m$m2 

}
else{
m1 <- m1-1
m2 <- m2-1 
}

theta <- opt.pqd(m1, m2, U1, U2)

# convert theta from vector to matrix
theta.matrix <- matrix(theta, m1, m2, byrow = T)

z <- matrix(0, 100, 100)
for (x in 1:100){
  for (y in 1:100){
    z[x, y] <- Cm.pqd(x.points[x], y.points[y], m1, m2, theta)
  }
}

x.points <- seq(0.001, 1, length.out = 100)
y.points <- seq(0.001, 1, length.out = 100)

# comtour plot of estimated copula
contour(x.points, y.points, z)

return(list('theta.matrix' = theta.matrix, 'm1' = m1+1, 'm2' = m2+1))

}
#######################################################################

# Example
# sample n=100 iid data from bivarate Frank copula with theta.true not equal to 0
n=100
theta.true <- 3
fr <- frankCopula(theta.true, 2)
X <- rCopula(n, fr)

copula.est(X, m1=NULL, m2=NULL, is.pqd = T)
