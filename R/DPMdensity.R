DPMdensity = function(y, 
                      ngrid=1000L, grid=NULL,
                      updateAlpha=TRUE, 
                      useHyperpriors=TRUE,
                      nclusters=50L, 
                      nskip=1000L, ndpost=1000L, keepevery=1L, printevery=1000L,
                      state=NULL, status=TRUE, 
                      alpha=10.0, a0=10.0, b0=1.0, 
                      m=colMeans(y), m0=colMeans(y), S0=NULL, 
                      lambda=0.5, gamma1=3.0, gamma2=2.0, 
                      nu=ncol(y)+2, Psi=NULL , nu0=ncol(y)+2, Psi0=NULL,
                      seed = 123
                      ) {
  
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
  if ((d == 2) & ((ngrid > 0) | !is.null(grid))) {
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
    prediction = FALSE
    ngrid = 0
    grid1 = grid2 = NULL
  }
  
  
  ##-----------------------------
  ## state, status
  ##-----------------------------
  if (status == FALSE) {
    ## use previous analysis
    nclusters = state$nclusters
    updateAlpha = state$updateAlpha
    a0 = state$a0
    b0 = state$b0
    alpha = state$alpha
    useHyperpriors = state$useHyperpriors
    m0 = state$m0
    S0 = state$S0
    m = state$m
    gamma1 = state$gamma1
    gamma2 = state$gamma2
    lambda = state$lambda
    nu0 = state$nu0
    Psi0 = state$Psi0
    nu = state$nu
    Psi = state$Psi
    Zeta = t(state$Zeta)
    Omega = state$Omega
    lw = state$lw
    a_gd = state$a_gd
    b_gd = state$b_gd
    kappa = state$kappa
  } else {
    ## start new analysis
    
    ##-----------------------------
    ## alpha ~ Gamma(a0, b0) or fixed
    ##-----------------------------
    if (updateAlpha) {
      if ((a0 > 0) & (b0 > 0)) 
        alpha = -1
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
      if (is.null(S0)) 
        S0 = diag(apply(y, 2, function(s) (range(s)[2]-range(s)[1])^2/16))
      m = NULL
      
      ### lambda ~ Gamma(gamma1, gamma2)
      if ((gamma1 > 0) & (gamma2 > 0))
        lambda = -1
      else
        stop("gamma1 and gamma2 are required to be positive scalars.")
      
      ### Psi ~ Wishart(nu0, Psi0)
      if (nu0 < d)
        stop("nu0 is required to be a scalar greater than ncol(y)-1.")
      if (is.null(Psi0))
        Psi0 = S0 / nu0
      Psi = NULL
      
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
      
      if (is.null(Psi)) {
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
  cat("Fitting a DPM of Multivariate Normals using Blocked Gibbs Sampling...", "\n")
  cat(" - Number of observations: ", n, "; Dimension: ", d, ".\n", sep = "")
  if(prediction)
    cat(" - Prediction = ", prediction, "; ngrid1 = ", ngrid, ", ngrid2 = ",  ngrid, ".\n", sep = "")
  else
    cat(" - Prediction = ", prediction, ".\n", sep = "")
  if(status)
    cat("Start a new analysis...", "\n", sep = "")
  else
    cat("Use previous analysis...", "\n", sep = "")
  cat(" - Number of clusters: ", nclusters, "; updateAlpha = ", updateAlpha, "; useHyperpriors = ", useHyperpriors, ".\n", sep = "")
  cat(" - Number of MCMC: ", nskip+ndpost*keepevery, ".\n", sep = "")
  
  
  #----------------------------------------------
  # set random seed
  #----------------------------------------------
  set.seed(seed = seed)
  
  
  #---------------------------------------------- 
  ## call Cpp function
  #---------------------------------------------- 
  ptm <- proc.time()
  
  res = .Call("_BNPqte_cDPMdensity",
              n,
              d,
              y,
              status,
              prediction,
              ngrid,
              updateAlpha,
              useHyperpriors,
              a0,
              b0,
              m0,
              S0,
              gamma1,
              gamma2,
              nu0,
              Psi0,
              nu,
              nclusters,
              nskip,
              ndpost,
              keepevery,
              printevery,
              alpha,
              lambda,
              m,
              Psi,
              Zeta,
              Omega,
              a_gd,
              b_gd,
              lw,
              kappa,
              grid1,
              grid2
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
  }
  
  return(res)
}




