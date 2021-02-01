DPMdensity = function(y, 
                      ngrid=1000L, grid=NULL,
                      updateAlpha=TRUE, 
                      useHyperpriors=TRUE,
                      nclusters=50L, 
                      nskip=1000L, ndpost=1000L, keepevery=1L, 
                      state=NULL, status=TRUE, 
                      alpha=2.0, a0=2.0, b0=1.0, 
                      m=colMeans(y), m0=colMeans(y), S0=NA, 
                      lambda=0.5, gamma1=3.0, gamma2=2.0, 
                      nu=ncol(y)+2, Psi=NA , nu0=ncol(y)+2, Psi0=NA,
                      seed = 123
                      ) {
  
  #----------------------------------------------
  # set random seed
  #----------------------------------------------
  set.seed(seed = seed)
  
  
  #---------------------------------------------- 
  # check and process arguments
  #---------------------------------------------- 
  
  ##-----------------------------
  ## y
  ##-----------------------------
  if (is.matrix(y) & (ncol(y) > 1)) {
    n = dim(y)[1]
    d = dim(y)[2]
  } else {
    stop("y is required to be a matrix with more than 1 column.")
  }
  
  ##-----------------------------
  ## ngrid, grid: only evaluate grid points when d=2
  ##-----------------------------
  prediction = FALSE
  if((ngrid <= 0) & is.null(grid)) {
    ## both of ngrid and grid are invalid -> no density prediction
    grid1 = grid2 = rep(0, 2)
  } else {
    if (d == 2) {
      prediction = TRUE
      if (is.null(grid)) {
        left = right = rep(0, 2)
        for (j in 1:2) {
          left[j] = min(y[, j]) - 0.5 * sd(y[, j])
          right[j] = max(y[, j]) + 0.5 * sd(y[, j])
        }
        ngrid = as.integer(sqrt(ngrid))
        grid1 = seq(left[1], right[1], length.out = ngrid)
        grid2 = seq(left[2], right[2], length.out = ngrid)
      } else {
        if (is.matrix(grid)) {
          ngrid = nrow(grid)
          grid1 = grid[, 1]
          grid2 = grid[, 2]
        } else {
          stop("grid is required to be a matrix or NULL.")
        }
      }
    } else {
      ## d!=2 -> no density prediction
      grid1 = grid2 = rep(0, 2)
    }
  }
  
  
  ##-----------------------------
  ## state, status
  ##-----------------------------
  if (status == FALSE) {
    ## use previous analysis
    updateAlpha = state$updateAlpha
    useHyperpriors = state$useHyperpriors
    nclusters = state$nclusters
    a0 = state$a0
    b0 = state$b0
    m0 = state$m0
    S0 = state$S0
    gamma1 = state$gamma1
    gamma2 = state$gamma2
    nu = state$nu
    nu0 = state$nu0
    Psi0 = state$Psi0
    alpha = state$alpha
    m = state$m
    lambda = state$lambda
    Psi = state$Psi
    Omega = state$Omega
    Zeta = t(state$Zeta)
    kappa = state$kappa
    lw = state$lw
    a_gd = state$a_gd
    b_gd = state$b_gd
  } else {
    ## start new analysis
    
    ##-----------------------------
    ## alpha ~ Gamma(a0, b0) or fixed
    ##-----------------------------
    if (updateAlpha) {
      if ((a0 > 0) & (b0 > 0)) 
        alpha = rgamma(n = 1, shape = a0, rate = b0)
      else
        stop("a0 and b0 are required to be positive scalars.")
    } else {
      if (alpha > 0) 
        a0 = b0 = -1
      else
        stop("alpha is required to be a positive scalar.")
    }
    
    ##-----------------------------
    ## Hyperpriors for the base distribution (Normal-Inverse-Wishart: N(zeta|m, Omega/lambda)xIW(Omega|nu, Psi))
    ##-----------------------------
    if (nu < d) 
      stop("nu is required to be a scalar greater than ncol(y)-1.")
    
    if (useHyperpriors) {
      ### m ~ Normal(m0, S0)
      if (!(is.vector(m0) & (length(m0) == d)))
        stop("m0 is required to be a vector of length equal to ncol(y).")
      if (is.na(S0)) {
        S0 = diag(apply(y, 2, function(s) (range(s)[2]-range(s)[1])^2/16))
      }
      m = rmvnorm(n = 1, mean = m0, sigma = S0)
      
      ### lambda ~ Gamma(gamma1, gamma2)
      if ((gamma1 > 0) & (gamma2 > 0))
        lambda = rgamma(n = 1, shape = gamma1, rate = gamma2)
      else
        stop("gamma1 and gamma2 are required to be positive scalars.")
      
      ### Psi ~ Wishart(nu0, Psi0)
      if (nu0 < d)
        stop("nu0 is required to be a scalar greater than ncol(y)-1.")
      if (is.na(Psi0))
        Psi0 = S0 / nu0
      Psi = rWishart(n = 1, df = nu0, Sigma = Psi0)[ , , 1]
      
    } else {
      ### m, lambda and Psi are fixed
      if (is.vector(m) & (length(m) == d)) {
        m0 = rep(-1, d)
        S0 = diag(-1, d)
      } else {
        stop("m is required to be a vector of length equal to ncol(y).")
      }
      
      if (lambda > 0) 
        gamma1 = gamma2 = -1
      else 
        stop("lambda is required to be a positive scalar.")
      
      if (is.na(Psi)) {
        nu0 = -1
        Psi0 = diag(-1, d)
        Psi = diag(apply(y, 2, function(s) (range(s)[2]-range(s)[1])^2/16))
      } else if (!is.positive.definite(Psi)) {
        stop("Psi is required to be a positive definite matrix.")
      }
      
    }
    
    Omega = Zeta = kappa = lw = a_gd = b_gd = NULL   # will initialize in cpp function
  }
  
  
  #---------------------------------------------- 
  ## print information
  #---------------------------------------------- 
  cat("Fitting a DPM of multivariate Normals using blocked Gibbs sampling...", "\n")
  cat("Number of clusters: ", nclusters, "\n", "updateAlpha = ", updateAlpha, "\n", "useHyperpriors = ", useHyperpriors, "\n", sep = "")
  
  
  #---------------------------------------------- 
  ## call Cpp function
  #---------------------------------------------- 
  ptm <- proc.time()
  
  res = .Call("_NonParamQTE_cDPMdensity",
              n,
              d,
              y,
              prediction,
              ngrid,
              grid1,
              grid2,
              updateAlpha,
              useHyperpriors,
              a0,
              b0,
              alpha,
              m0,
              S0,
              m,
              gamma1,
              gamma2,
              lambda,
              nu0,
              Psi0,
              nu,
              Psi,
              nclusters,
              nskip,
              ndpost,
              keepevery,
              Zeta,
              Omega,
              a_gd,
              b_gd,
              lw,
              kappa
  )
  
  cat("Finished!", "\n")
  
  
  #---------------------------------------------- 
  # returns
  #---------------------------------------------- 
  res$proc.time = proc.time() - ptm
  res$prediction = prediction
  if(prediction) {
    res$grid1 = grid1
    res$grid2 = grid2
    res$predict.densities.mean = colMeans(res$predict.densities)
  }
  
  return(res)
}




