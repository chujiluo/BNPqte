#---------------------------------------------------
# bart: mixed-type predictors
#---------------------------------------------------
MixData = function(n, p, sigma, binary) {
  k = p/2
  X = matrix(NA, nrow = n, ncol = p)
  X[, 1:k] = matrix(rbinom(n*k, size = 1, prob = 0.5), nrow = n, ncol = k)
  X[, (k+1):p] = matrix(runif(n*k), nrow = n, ncol = k)
  
  #f0 = 10 * sin(pi * X[, k+1] * X[, k+2]) + 20 * (X[, k+3] - 0.5)^2 + 10 * X[, 1] + 5 * X[, 2]
  f0 = 10 * sin(pi * X[, 1] * X[, k+1]) + 20 * (X[, k+3] - 0.5)^2 + 10 * X[, 2] + 5 * X[, k+2]
  X = as.data.frame(X)
  for (j in 1:k) {
    X[, j] = factor(X[, j], levels = c("0", "1"), labels = c("0", "1"))
  }
  dimnames(X)[[2]] = sapply(1:p, function(s) paste("X", s, sep=""))
  
  if (binary) {
    f0 = scale(f0, center = T, scale = F)
    p = pnorm(f0, mean = 0, sd = 1)
    Y = c()
    for (i in 1:n) {
      Y[i] = rbinom(1, 1, prob = p[i])
    }
    return(list(X=X, Y=Y, f0=f0, p=p))
    
  } else {
    Y = rnorm(n, mean = f0, sd = sigma)
    return(list(X=X, Y=Y, f0=f0, sigma=sigma))
  }
  
}

#---------------------------------------------------
# dpmm: mixture of three bivariate Normals
#---------------------------------------------------
ThreeNormals = function(n) {

  zeta1 = c(2, -1)
  zeta2 = c(1, 0)
  zeta3 = c(-1, -1)
  
  y = matrix(NA, nrow = n, ncol = 2)
  
  for (i in 1:n) {
    tmp = sample(1:3, 1)
    if (tmp==1){
      y[i, ] = rmvnorm(1, mean = zeta1, sigma = 0.5*diag(2))
    }
    if (tmp==2){
      y[i, ] = rmvnorm(1, mean = zeta2, sigma = 0.5*diag(2))
    }
    if (tmp==3){
      y[i, ] = rmvnorm(1, mean = zeta3, sigma = 0.5*diag(2))
    }
  }
  
  ## joint density of (y1, y2)
  dtrue = function(y1, y2) {
    z = c(y1, y2)
    out = (1/3)*dmvnorm(z, zeta1, 0.5*diag(2)) + (1/3)*dmvnorm(z, zeta2, 0.5*diag(2)) + (1/3)*dmvnorm(z, zeta3, 0.5*diag(2))
    
    return(out)
  }
  
  res = list(y = y, dtrue = dtrue)
  return(res)
}

#---------------------------------------------------
# dpmm: example from Dunson et al. (2007)
#---------------------------------------------------
DunsonExample = function(n) {
  x = runif(n)
  y1 = x + rnorm(n, 0, sqrt(0.01))
  y2 = x^4 + rnorm(n, 0, sqrt(0.04))
  u = runif(n)
  prob = exp(-2 * x)
  y = ifelse(u < prob, y1, y2)
  
  ## conditional density of grid given x
  dtrue = function(grid, x) {
    exp(-2 * x) * dnorm(grid, mean = x, sd = sqrt(0.01)) + (1 - exp(-2 * x)) * dnorm(grid, mean = x^4, sd = sqrt(0.04))
  }
  ## conditional mean of grid given x
  mtrue <- function(x) exp(-2 * x) * x + (1 - exp(-2 * x)) * x^4
  ## conditional CDF of grid given x
  ptrue = function(grid, x) {
    exp(-2 * x) * pnorm(grid, mean = x, sd = sqrt(0.01)) + (1 - exp(-2 * x)) * pnorm(grid, mean = x^4, sd = sqrt(0.04))
  }
  
  res = list(y = y, x = x, mtrue = mtrue, dtrue = dtrue, ptrue = ptrue)
  return(res)
}

#---------------------------------------------------
# qte: simulation 2 from Xu et al. (2018)
#---------------------------------------------------
QteExample2 = function(n) {
  ## x
  x1 = matrix(runif(n*5, 0, 1), ncol = 5)
  x2 = matrix(runif(n*5, 1, 2), ncol = 5)
  x3 = matrix(rbinom(n*10, 1, 0.5), ncol = 10)
  x = cbind(x1, x2, x3)
  
  treatment = c()
  y0 = c()
  y1 = c()
  y = c()
  
  for (i in 1:n) {
    ## treatment
    p = 1 / (1 + exp(2.125 - 0.5*x[i, 1]*x[i, 6] - x[i, 2]*x[i, 7] - 0.2*sum(x[i, 1:10]*x[i, 11:20])))
    treatment[i] = rbinom(1, 1, p)
    
    ## y0
    w0 = 1 / (1 + exp(-2*x[i, 2]*x[i, 7] + 2))
    u0 = runif(1)
    if(u0 < w0) {
      y0[i] = rnorm(1, mean = sum(x[i, 1:5]*x[i, 11:15]), sd = 2)
    } else {
      y0[i] = rnorm(1, mean = (8 + 0.6*sum(x[i, 6:10]*x[i, 16:20])), sd = 2)
    }
    
    ## y1
    w1 = 1 / (1 + exp(-2 + 2*x[i, 1]*x[i, 6]))
    u1 = runif(1)
    if(u1 < w1) {
      y1[i] = rnorm(1, sum(x[i, c(1, 3, 5, 7, 9)]*x[i, c(11, 13, 15, 17, 19)]), 2)
    } else {
      y1[i] = rnorm(1, 2 + 1.7*sum(x[i, c(2, 4, 6, 8, 10)]*x[i, c(12, 14, 16, 18, 20)]), 2)
    }
    
    ## y
    y[i] = treatment[i]*y1[i] + (1-treatment[i])*y0[i]
    
  }
  
  ## true density of potential outcome Y1: f(Y1)
  fy1 = function(y, k) {
    outcome = 0
    for(i in 1:k) {
      x1 = runif(5, 0, 1)
      x2 = runif(5, 1, 2)
      x3 = rbinom(10, 1, 0.5)
      x = c(x1, x2, x3)
      
      w1 = 1 / (1 + exp(-2 + 2*x[1]*x[6]))
      outcome = outcome + w1*dnorm(y, sum(x[c(1, 3, 5, 7, 9)]*x[c(11, 13, 15, 17, 19)]), 2) + 
        (1-w1)*dnorm(y, 2 + 1.7*sum(x[c(2, 4, 6, 8, 10)]*x[c(12, 14, 16, 18, 20)]), 2)
    }
    return(outcome/k)
  }
  
  ## true density of potential outcome Y0: f(Y0)
  fy0 = function(y, k) {
    outcome = 0
    for(i in 1:k) {
      x1 = runif(5, 0, 1)
      x2 = runif(5, 1, 2)
      x3 = rbinom(10, 1, 0.5)
      x = c(x1, x2, x3)
      
      w0 = 1 / (1 + exp(-2*x[2]*x[7] + 2))
      outcome = outcome + w0*dnorm(y, mean = sum(x[1:5]*x[11:15]), sd = 2) + 
        (1-w0)*dnorm(y, mean = (8 + 0.6*sum(x[6:10]*x[16:20])), sd = 2)
    }
    return(outcome/k)
  }
  
  res = list(x=x, treatment=treatment, y0=y0, y1=y1, y=y, fy1 = fy1, fy0 = fy0)
  return(res)
}

#---------------------------------------------------
# qte: simulation 3 from Xu et al. (2018)
#---------------------------------------------------
QteExample = function(n) {
  ## x
  x = matrix(runif(n*10, -2, 2), ncol = 10)
  
  treatment = c()
  y0 = c()
  y1 = c()
  y = c()
  
  for (i in 1:n) {
    ## treatment
    p = 1 / (1 + exp(-0.3 * sum(x[i, ])))
    treatment[i] = rbinom(1, 1, p)
    
    ## y0
    w0 = exp(-abs(x[i, 5]))
    u0 = runif(1)
    if(u0 < w0) {
      y0[i] = rnorm(1, mean = sum(0.2*x[i, 1:5])^4, sd = 1)
    } else {
      y0[i] = rnorm(1, mean = (2 + sum(0.2*x[i, 1:5]*x[i, 1:5])), sd = 1)
    }
    
    ## y1
    w1 = 1 / (1 + exp(-0.5*x[i, 3]*x[i, 4]))
    u1 = runif(1)
    if(u1 < w1) {
      y1[i] = rnorm(1, (3 + 0.5*x[i, 2]*x[i, 5] + 0.5*x[i, 1]^2), 0.5)
    } else {
      y1[i] = rnorm(1, (-0.5 + 0.5*x[i, 2]^2 - 0.5*x[i, 1]*x[i, 3]), 0.8)
    }
    
    ## y
    y[i] = treatment[i]*y1[i] + (1-treatment[i])*y0[i]
    
  }
  
  ## true density of potential outcome Y1: f(Y1)
  fy1 = function(y, k) {
    outcome = 0
    for(i in 1:k) {
      x = runif(10, -2, 2)
      
      w1 = 1 / (1 + exp(-0.5*x[3]*x[4]))
      outcome = outcome + w1*dnorm(y, (3 + 0.5*x[2]*x[5] + 0.5*x[1]^2), 0.5) + 
        (1-w1)*dnorm(y, (-0.5 + 0.5*x[2]^2 - 0.5*x[1]*x[3]), 0.8)
    }
    return(outcome/k)
  }
  
  ## true density of potential outcome Y0: f(Y0)
  fy0 = function(y, k) {
    outcome = 0
    for(i in 1:k) {
      x = runif(10, -2, 2)
      
      w0 = exp(-abs(x[5]))
      outcome = outcome + w0*dnorm(y, mean = sum(0.2*x[1:5])^4, sd = 1) + 
        (1-w0)*dnorm(y, mean = (2 + sum(0.2*x[1:5]*x[1:5])), sd = 1)
    }
    return(outcome/k)
  }
  
  res = list(x=x, treatment=treatment, y0=y0, y1=y1, y=y, fy1 = fy1, fy0 = fy0)
  return(res)
}