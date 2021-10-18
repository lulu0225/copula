#######################################################################
# function for testing PQD

# X: a n*2 matrix of data
# m1, m2: degrees of Bernstein polynomial estimated by grid search if m1 = NULL and m2 = NULL

#######################################################################

# check if the estimation function "copula.est" exists
exists("copula.est") 

copula.test <- function(X, m1 = NULL, m2 = NULL){

# Cm = Gm^T theta + XXX + YYY
Gm.pqd <- function (u1, u2, m1, m2) {
  
  G1 <- c()
  G2 <- c()
  
  for (k1 in 1:m1){
    G1[k1] <- choose((m1+1),k1) * u1^k1 * (1 - u1)^(m1 - k1 + 1)
  }
  
  for (k2 in 1:m2){
    G2[k2] <- choose((m2+1),k2) * u2^k2 * (1 - u2)^(m2 - k2 + 1)
  }
  
  G0 <- G1 %o% G2
  
  as.vector(t(G0))
}

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

# Cm = Gm^T theta + XXX + YYY
Cm.pqd <- function(u1, u2, m1, m2, theta){
  
  pro <- crossprod(Gm.pqd(u1, u2, m1, m2), theta) + XXX(u1, u2, m1, m2) + YYY(u1, u2, m1, m2)
  pro
  
}

  ###################################################
  ###################################################
  U1 <- X[,1]
  U2 <- X[,2]

  n <- length(U1)

  # transform to pseudo-observations
  ec1 <- ecdf(U1)
  ec2 <- ecdf(U2)
  U1 <- n / (n+1) * ec1(U1)
  U2 <- n / (n+1) * ec2(U2)
  
  # Estimation without PQD constraint
  EstWithoutPQD <- copula.est(X, m1, m2, is.pqd = F, print.contour = F)
  
  m1.nonpqd <- EstWithoutPQD$m1
  m2.nonpqd <- EstWithoutPQD$m2
  
  theta.nonpqd <- as.vector(t(EstWithoutPQD$theta.matrix))
  
  ###################################################
  ###################################################
  # KS distance
  
  # test statistic, C2 is unconstrained
  KS <- function (U1, U2, m1.nonpqd, m2.nonpqd, theta.nonpqd){
    n <- length(U1)
    
    x <- seq(0.05, 0.95, by = 0.05)
    y <- seq(0.05, 0.95, by = 0.05)
    
    sup <- -Inf
    for (u in x){
      for (v in y){
        
        C2 <- Cm.pqd(u, v, m1.nonpqd, m2.nonpqd, theta.nonpqd)
        ma <- sqrt(n) * (u*v - C2)
        sup <- max(ma,sup)
      }
    }
    sup
  }
  
  # test statistic
  ks.test <- KS(U1, U2, m1.nonpqd-1, m2.nonpqd-1, theta.nonpqd)

  ###################################################
  ###################################################
  #CVM distance
  
  CVM.2 <- function (U1, U2, m1.nonpqd, m2.nonpqd, theta.nonpqd){
    n <- length(U1)
 
    find.integral <- function(u){
      u1 <- u[1]
      u2 <- u[2]
      
      theta.matrix <- matrix(theta.nonpqd, m1.nonpqd, m2.nonpqd, byrow = T)
      theta.left <- rep(0, m1.nonpqd)
      theta.right <- seq(1/ (m1.nonpqd+1), 1-1/(m1.nonpqd+1), 1/(m1.nonpqd+1))
      
      matrix.1 <- cbind(theta.left, theta.matrix, theta.right)
      
      theta.upper <- rep(0, m2.nonpqd+2)
      theta.lower <- seq(0, 1, 1/(m2.nonpqd+1))
      
      theta.m <- rbind(theta.upper, matrix.1, theta.lower)
      ##############################
      
      g <- matrix(0, m1.nonpqd+1 ,m2.nonpqd+1)
      w <- matrix(0, m1.nonpqd+1 ,m2.nonpqd+1)
      
      for (k1 in 0:m1.nonpqd){
        for (k2 in 0: (m2.nonpqd)){
          
          g[k1+1, k2+1] <- (m1.nonpqd+1) * (m2.nonpqd+1) * choose(m1.nonpqd,k1) * choose(m2.nonpqd,k2) * u1^k1 * (1 - u1)^(m1.nonpqd - k1) * u2^k2 * (1 - u2)^(m2.nonpqd  - k2)
          w[k1+1, k2+1] <- theta.m[k1+2, k2+2] - theta.m[k1+2, k2+1] - theta.m[k1+1, k2+2] + theta.m[k1+1, k2+1]
          
        }
      }
      g.vector <- as.vector(t(g))
      w.vector <- as.vector(t(w))
      gr <- crossprod(g.vector, w.vector)
      
      n * (max(0, u1*u2 - Cm.pqd(u1, u2, m1.nonpqd, m2.nonpqd, theta.nonpqd)))^2 * gr
    }
    
    x <- seq(0.05, 1, by = 0.05)
    y <- seq(0.05, 1, by = 0.05)
    
    sum <- 0
    for (u in x){
      for (v in y){
        sum = sum + find.integral(c(u, v))
      }
    }
    sum/(20*20)
  }
  
  # test statistic
  cvm.test <- CVM.2(U1, U2, m1.nonpqd-1, m2.nonpqd-1, theta.nonpqd)

  ###################################################
  ###################################################
  #AD distance 
  
  AD.2 <- function (U1, U2, m1.nonpqd, m2.nonpqd, theta.nonpqd){
    n <- length(U1)
  
    find.integral <- function(u){
      u1 <- u[1]
      u2 <- u[2]
      
      theta.matrix <- matrix(theta.nonpqd, m1.nonpqd, m2.nonpqd, byrow = T)
      theta.left <- rep(0, m1.nonpqd)
      theta.right <- seq(1/ (m1.nonpqd+1), 1-1/(m1.nonpqd+1), 1/(m1.nonpqd+1))
      
      matrix.1 <- cbind(theta.left, theta.matrix, theta.right)
      
      theta.upper <- rep(0, m2.nonpqd+2)
      theta.lower <- seq(0, 1, 1/(m2.nonpqd+1))
      
      theta.m <- rbind(theta.upper, matrix.1, theta.lower)
      ##############################
      
      g <- matrix(0, m1.nonpqd+1 ,m2.nonpqd+1)
      w <- matrix(0, m1.nonpqd+1 ,m2.nonpqd+1)
      
      for (k1 in 0:m1.nonpqd){
        for (k2 in 0: (m2.nonpqd)){
          
          g[k1+1, k2+1] <- (m1.nonpqd+1) * (m2.nonpqd+1) * choose(m1.nonpqd,k1) * choose(m2.nonpqd,k2) * u1^k1 * (1 - u1)^(m1.nonpqd - k1) * u2^k2 * (1 - u2)^(m2.nonpqd  - k2)
          w[k1+1, k2+1] <- theta.m[k1+2, k2+2] - theta.m[k1+2, k2+1] - theta.m[k1+1, k2+2] + theta.m[k1+1, k2+1]
          
        }
      }
      g.vector <- as.vector(t(g))
      w.vector <- as.vector(t(w))
      gr <- crossprod(g.vector, w.vector)
      
      n * (max(0, u1*u2 - Cm.pqd(u1, u2, m1.nonpqd, m2.nonpqd, theta.nonpqd)))^2 /(u1 * u2 * (1-u1) * (1-u2)) * gr
    }
    
    x <- seq(0.05, .95, by = 0.05)
    y <- seq(0.05, .95, by = 0.05)
    
    sum <- 0
    for (u in x){
      for (v in y){
        
        sum = sum + find.integral(c(u, v))
      }
    }
    sum/(19*19)
  }
  
  # test statistic
  ad.test <- AD.2(U1, U2, m1.nonpqd-1, m2.nonpqd-1, theta.nonpqd)

  ###################################################
  ###################################################
  
  #draw samples from C1 which is PQD constrained to approximate 
  #the distribution of test statistic under the null of PQD
  
  # Estimation with PQD constraint
  EstWithPQD <- copula.est(X, m1, m2, is.pqd = T, print.contour = F)
  
  m1.pqd <- EstWithPQD$m1
  m2.pqd <- EstWithPQD$m2
  
  theta.pqd <- as.vector(t(EstWithPQD$theta.matrix))

  #convert theta to w 
  theta.to.w <- function(theta.pqd, m1.pqd, m2.pqd){
    theta.matrix <- matrix(theta.pqd, m1.pqd, m2.pqd, byrow = T)
    
    theta.left <- rep(0, m1.pqd)
    theta.right <- seq(1/ (m1.pqd+1), 1-1/(m1.pqd+1), 1/(m1.pqd+1))
    
    matrix.1 <- cbind(theta.left, theta.matrix, theta.right)
    
    theta.upper <- rep(0, m2.pqd+2)
    theta.lower <- seq(0, 1, 1/(m2.pqd+1))
    
    theta.m <- rbind(theta.upper, matrix.1, theta.lower)
    
    #calculate w which is the weight in the pdf
    der <- function(k1, k2, theta.m){
      theta.m[k1+1, k2+1] - theta.m[k1, k2+1] - theta.m[k1+1, k2] + theta.m[k1, k2]
    }
    
    w <- matrix(0, m1.pqd+1, m2.pqd+1)
    
    for (k1 in 1:(m1.pqd+1)){
      for (k2 in 1:(m2.pqd+1)){
        w[k1, k2] <- max(0, der(k1, k2, theta.m)) #positive probability
      }
    }
    as.vector(t(w))
  }
  
  w.vector <- theta.to.w(theta.pqd, m1.pqd-1, m2.pqd-1)
  
  ##########################
  
  #draw M=1000 samples from C1
  M = 1000
  ks.sim <- c()
  cvm.sim <- c()
  ad.sim <- c()
  
  for (j in 1 : M){
    
    MM <- sample(1:length(w.vector), n, prob = w.vector, replace = T)
    k1.pqd <- floor((MM-1)/(m2.pqd+1))
    k2.pqd <- MM - k1.pqd * (m2.pqd+1) - 1
    
    U1.pqd <- c()
    U2.pqd <- c()
    
    for (k in 1:n){
      U1.pqd[k] <- rbeta(1, k1.pqd[k]+1, m1.pqd + 1 - k1.pqd[k])
      U2.pqd[k] <- rbeta(1, k2.pqd[k]+1, m2.pqd + 1 - k2.pqd[k])
    }
    
    # Estimation without PQD constraint using simulated PQD data
    EstWithoutPQD.pqd <- copula.est(cbind(U1.pqd,U2.pqd), m1 = m1.nonpqd, m2 = m2.nonpqd, is.pqd = F, print.contour = F)
    
    theta.nonpqd.pqd <- as.vector(t(EstWithoutPQD.pqd$theta.matrix))
    
    #the distribution of test statistics under the null of PQD
   ks.sim[j] <- KS(U1.pqd, U2.pqd, m1.nonpqd-1, m2.nonpqd-1, theta.nonpqd.pqd)
   cvm.sim[j] <- CVM.2(U1.pqd, U2.pqd, m1.nonpqd-1, m2.nonpqd-1, theta.nonpqd.pqd)
   ad.sim[j] <- AD.2(U1.pqd, U2.pqd, m1.nonpqd-1, m2.nonpqd-1, theta.nonpqd.pqd)

  }
  
  return(list('ks.test' = ks.test, 
              'ks.quantiles'= quantile(ks.sim, c(0.025, 0.05, 0.1, 0.5, 0.90, 0.95, 0.975)),
              'cvm.test' = cvm.test, 
              'cvm.quantiles'= quantile(cvm.sim, c(0.025, 0.05, 0.1, 0.5, 0.90, 0.95, 0.975)),
              'ad.test' = ad.test, 
              'ad.quantiles'= quantile(ad.sim, c(0.025, 0.05, 0.1, 0.5, 0.90, 0.95, 0.975))
              ))
  
}


###################################################
###################################################

# Examples
# install.packages("copula")
library(copula)

#Example 1: NQD
#sample n=200 iid data from bivarate Frank copula with theta = -1
set.seed(123)
n = 200
fr <- frankCopula(-1, 2)
X <- rCopula(n, fr)

copula.test(X, m1 = NULL, m2 = NULL)


#Example 2: PQD
#sample n=200 iid data from bivarate Frank copula with theta = 1
set.seed(123)
n = 200
fr <- frankCopula(1, 2)
X <- rCopula(n, fr)

copula.test(X, m1 = NULL, m2 = NULL)

