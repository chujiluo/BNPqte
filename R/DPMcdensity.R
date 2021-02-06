DPMcdensity = function(y, x, 
                       xpred=NULL, ngrid=1000L, grid=NULL, 
                       type="pdf", compute.band=TRUE,
                       updateAlpha=TRUE, 
                       useHyperpriors=TRUE,
                       nclusters=50L, 
                       nskip=1000L, ndpost=1000L, keepevery=1L, 
                       state=NULL, status=TRUE, 
                       alpha=10.0, a0=10.0, b0=1.0, 
                       m=colMeans(cbind(y, x)), m0=colMeans(cbind(y, x)), S0=NULL, 
                       lambda=0.5, gamma1=3.0, gamma2=2.0, 
                       nu=ncol(as.matrix(x))+3, Psi=NULL , nu0=ncol(as.matrix(x))+3, Psi0=NULL,
                       seed = 123
) {
  
  #---------------------------------------------- 
  # check and process arguments
  #---------------------------------------------- 
  
  ##-----------------------------
  ## y
  ##-----------------------------
  if(is.vector(y)) {
    n = length(y)
  } else {
    stop("y is required to be a vector.")
  }
  
  ##-----------------------------
  ## x
  ##-----------------------------
  if(is.vector(x)) {
    if(length(x) != n)
      stop("The length of x is required to be the same as length(y), when x is a vector.")
    else
      d = 2
  } else {
    if(is.matrix(x)) {
      if(nrow(x) != n)
        stop("nrow(x) is required to be the same as length(y), when x is a matrix")
      else
        d = ncol(x) + 1
    } else {
      stop("x is required to be a vector or matrix.")
    }
  }
  data = cbind(y, x)  # n x d matrix (d >= 2)
  x = as.matrix(x)
  
  ##-----------------------------
  ## ngrid, grid, type, xpred
  ##-----------------------------
  if(is.null(xpred)) {
    pdf = cdf = FALSE
    ngrid = npred = 0
    grid = xpred = NULL
  } else {
    if((ngrid <= 0) & is.null(grid)) {
      pdf = cdf = FALSE
      ngrid = npred = 0
      grid = xpred = NULL
    } else {
      xpred = as.matrix(xpred)
      npred = nrow(xpred)
      if(ncol(xpred) != (d-1))
        stop("xpred is required to have the same number of columns as x.")
      
      if(!(type %in% c("pdf", "cdf", "both"))) {
        stop("type is required to be one of pdf, cdf and both.")
      } else {
        if(type == "pdf") {
          pdf = TRUE
          cdf = FALSE
        }
        if(type == "cdf") {
          pdf = FALSE
          cdf = TRUE
        }
        if(type == "both")
          pdf = cdf = TRUE
      }
      
      if(is.null(grid)) {
        ### ngrid > 0
        grid <- seq(from = (min(y) - 0.25 * sd(y)), to = (max(y) + 0.25 * sd(y)), length.out = ngrid)
      } else {
        ### grid is provided
        grid = as.vector(grid)
        ngrid = length(grid)
      }
    }
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
      stop("nu is required to be a scalar greater than ncol(cbind(y, x))-1.")
    
    if (useHyperpriors) {
      ### m ~ Normal(m0, S0)
      if (!(is.vector(m0) & (length(m0) == d)))
        stop("m0 is required to be a vector of length equal to ncol(cbind(y, x)).")
      if (is.null(S0)) 
        S0 = diag(apply(data, 2, function(s) (range(s)[2]-range(s)[1])^2/16))
      m = NULL
      
      ### lambda ~ Gamma(gamma1, gamma2)
      if ((gamma1 > 0) & (gamma2 > 0))
        lambda = -1
      else
        stop("gamma1 and gamma2 are required to be positive scalars.")
      
      ### Psi ~ Wishart(nu0, Psi0)
      if (nu0 < d)
        stop("nu0 is required to be a scalar greater than ncol(cbind(y, x))-1.")
      if (is.null(Psi0))
        Psi0 = S0 / nu0
      Psi = NULL
      
    } else {
      ### m, lambda and Psi are fixed
      if (is.vector(m) & (length(m) == d)) {
        m0 = rep(-1, d)
        S0 = diag(-1, d)
      } else {
        stop("m is required to be a vector of length equal to ncol(cbind(y, x)).")
      }
      
      if (lambda > 0) 
        gamma1 = gamma2 = -1
      else 
        stop("lambda is required to be a positive scalar.")
      
      if (is.null(Psi)) {
        nu0 = -1
        Psi0 = diag(-1, d)
        Psi = diag(apply(data, 2, function(s) (range(s)[2]-range(s)[1])^2/16))
      } else if (!is.positive.definite(Psi)) {
        stop("Psi is required to be a positive definite matrix.")
      }
      
    }
    
    Omega = Zeta = kappa = lw = a_gd = b_gd = NULL
  }
  
  
  #---------------------------------------------- 
  ## print information
  #---------------------------------------------- 
  cat("Fitting a Weighted Dependent DPM of Multivariate Normals using Blocked Gibbs Sampling...", "\n")
  cat("- Number of observations: ", n, "; Number of covariates: ", d-1, ".\n", sep = "")
  if(pdf | cdf)
    cat("- Prediction = TRUE; Prediction type = ", type, "; ngrid(y) = ", ngrid, "; ngrid(x) = ", npred, ".\n", sep = "")
  else
    cat("- Prediction = FALSE.\n", sep = "")
  if(status)
    cat("- Start a new analysis.", "\n", sep = "")
  else
    cat("- Use previous analysis.", "\n", sep = "")
  cat("- Number of clusters: ", nclusters, "; updateAlpha = ", updateAlpha, "; useHyperpriors = ", useHyperpriors, ".\n", sep = "")
  
  
  #----------------------------------------------
  # set random seed
  #----------------------------------------------
  set.seed(seed = seed)
  
  
  #---------------------------------------------- 
  ## call Cpp function
  #---------------------------------------------- 
  ptm <- proc.time()
  
  res = .Call("_BNPqte_cDPMcdensity",
              n,
              d,
              data,
              y,
              x,
              status,
              pdf,
              cdf,
              ngrid,
              npred,
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
              grid,
              xpred
  )
  
  cat("Finished!", "\n")
  
  
  #---------------------------------------------- 
  # returns
  #---------------------------------------------- 
  res$proc.time = proc.time() - ptm
  res$pdf = pdf
  res$cdf = cdf
  if(pdf | cdf) {
    res$xpred = xpred;
    res$grid = grid;
    if(pdf) {
      res$predict.pdf.mean = apply(simplify2array(res$predict.pdf), c(1, 2), mean)
      if(compute.band) {
        #
      }
    }
    if(cdf) {
      res$predict.cdf.mean = apply(simplify2array(res$predict.cdf), c(1, 2), mean)
      if(compute.band) {
        #
      }
    }
  }
  
  return(res)
}