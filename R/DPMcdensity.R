DPMcdensity = function(y, x, 
                       method="truncated", nclusters=50L,
                       xpred=NULL, ngrid=1000L, grid=NULL, 
                       type.pred=c("pdf"), compute.band=TRUE, type.band="HPD",
                       updateAlpha=TRUE, useHyperpriors=TRUE,
                       status=TRUE, state=NULL, 
                       nskip=1000L, ndpost=1000L, keepevery=1L, printevery=1000L,
                       alpha=10.0, a0=10.0, b0=1.0, 
                       m=NULL, m0=NULL, S0=NULL, 
                       lambda=0.5, gamma1=3.0, gamma2=2.0, 
                       nu=NULL, Psi=NULL, nu0=NULL, Psi0=NULL,
                       diag=FALSE,
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
  ## ngrid, grid, type.pred, xpred, compute.band, type.band
  ##-----------------------------
  pdf = cdf = meanReg = FALSE
  if((!is.null(xpred)) & ((ngrid > 0) | (!is.null(grid)))) {
    xpred = as.matrix(xpred)
    npred = nrow(xpred)
    
    if(ncol(xpred) != (d-1))
      stop("xpred is required to have the same number of columns as x.")
    
    type.pred = as.vector(type.pred)
    if(!all(type.pred %in% c("pdf", "cdf", "meanReg"))) {
      stop("Available type.pred are: pdf, cdf and meanReg.")
    } else {
      if("pdf" %in% type.pred)
        pdf = TRUE
      if("cdf" %in% type.pred)
        cdf = TRUE
      if("meanReg" %in% type.pred)
        meanReg = TRUE
    }
    
    if(is.null(grid)) {
      ### ngrid > 0
      grid = seq(from = (min(y) - 0.25 * sd(y)), to = (max(y) + 0.25 * sd(y)), length.out = ngrid)
    } else {
      ### grid is provided
      grid = as.vector(grid)
      ngrid = length(grid)
    }
  } else {
    ngrid = npred = 0
    grid = xpred = NULL
  }
  
  hpd = bci = FALSE
  if(compute.band) {
    if(!(type.band %in% c("HPD", "BCI")))
      stop("Only two available bands are provided: HPD (Highest posterior interval ) and BCI (Bayesian credible interval).")
    else
      if(type.band == "HPD")
        hpd = TRUE
      else
        bci = TRUE
  }

  
  ##-----------------------------
  ## method
  ##-----------------------------
  if(!(method %in% c("truncated", "neal")))
    stop("Only two available sampling methods: truncated or neal.")
  
  
  ##-----------------------------
  ## state, status
  ##-----------------------------
  if (status == FALSE) {
    ## use previous analysis
    method = state$method
    
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
    kappa = state$kappa
    
    if(method == "truncated") {
      lw = state$lw
      a_gd = state$a_gd
      b_gd = state$b_gd
    }
    
  } else {
    ## start new analysis
    
    ##-----------------------------
    ## alpha ~ Gamma(a0, b0) or fixed
    ##-----------------------------
    if (updateAlpha) {
      if ((a0 > 0) & (b0 > 0)) 
        alpha = 1.0   # initialize
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
    if(is.null(nu)) {
      nu = ncol(data) + 2
    } else {
      if (nu < d) 
        stop("nu is required to be a scalar greater than ncol(cbind(y, x))-1.")
    }
    
    if (useHyperpriors) {
      ### m ~ Normal(m0, S0)
      if(is.null(m0)) {
        m0 = colMeans(data)
      } else {
        if (!(is.vector(m0) & (length(m0) == d)))
          stop("m0 is required to be a vector of length equal to ncol(cbind(y, x)).")
      }
      if (is.null(S0)) 
        S0 = diag(apply(data, 2, function(s) (range(s)[2]-range(s)[1])^2/16))
      m = m0 + rnorm(d, 0, 100)   # initialize
      
      ### lambda ~ Gamma(gamma1, gamma2)
      if ((gamma1 > 0) & (gamma2 > 0))
        lambda = rgamma(1, shape = gamma1, rate = gamma2)  # initialize
      else
        stop("gamma1 and gamma2 are required to be positive scalars.")
      
      ### Psi ~ Wishart(nu0, Psi0)
      if(is.null(nu0)) {
        nu0 = ncol(data) + 2
      } else {
        if (nu0 < d) 
          stop("nu0 is required to be a scalar greater than ncol(cbind(y, x))-1.")
      }
      if (is.null(Psi0))
        Psi0 = S0 / nu0
      Psi = nu0 * Psi0   # initialize
      
    } else {
      ### m, lambda and Psi are fixed
      if(is.null(m)) {
        m = colMeans(data)
      } else {
        if (is.vector(m) & (length(m) == d)) {
          m0 = rep(-1, d)
          S0 = diag(-1, d)
        } else {
          stop("m is required to be a vector of length equal to ncol(cbind(y, x)).")
        }
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
    
    if(method == "truncated") {
      a_gd = rep(1.0, (nclusters-1))
      b_gd = rep(alpha, (nclusters-1))
      lw = NULL
    }
    
    Omega = Zeta = kappa = NULL   # will initialize in cpp function
  }
  
  
  #---------------------------------------------- 
  ## print information
  #---------------------------------------------- 
  cat("*****Into main of Weight-Dependent DPMM\n")
  cat("*****Data: n, d: ", n, ", ", d, "\n", sep = "")
  if(any(pdf, cdf, meanReg))
    cat("*****Prediction: type, ngrid, nxpred: ", paste(type.pred, collapse = ", "), ", ", ngrid, ", ", npred, "\n", sep = "")
  else
    cat("*****Prediction: FALSE\n")
  if(method == "truncated") {
    cat("*****Posterior sampling method: Blocked Gibbs Sampling with", nclusters, "clusters\n")
  } else {
    cat("*****Posterior sampling method: Algorithm 8 with m = 1 in Neal (2000)\n")
  }
  cat("*****Prior: updateAlpha, useHyperpriors: ", updateAlpha, ", ", useHyperpriors, "\n", sep="")
  cat("*****MCMC: nskip, ndpost, keepevery, printevery: ", nskip, ", ", ndpost, ", ", keepevery, ", ", printevery, "\n", sep = "")
  if(status)
    cat("*****Start a new MCMC...", "\n", sep = "")
  else
    cat("*****Continue previous MCMC...", "\n", sep = "")
  
  
  #----------------------------------------------
  # set random seed
  #----------------------------------------------
  set.seed(seed = seed)
  
  
  #---------------------------------------------- 
  ## call Cpp function
  #---------------------------------------------- 
  ptm <- proc.time()
  
  if(method == "truncated") {
    res = .Call("_BNPqte_cDPMcdensity",
                n,
                d,
                data,
                y,
                x,
                status,
                diag,
                pdf,
                cdf,
                meanReg,
                ngrid,
                npred,
                hpd,
                bci,
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
                a_gd,
                b_gd,
                Zeta,
                Omega,
                lw,
                kappa,
                grid,
                xpred
    )
  } else {
    res = .Call("_BNPqte_cDPMcdensityNeal",
                n,
                d,
                data,
                y,
                x,
                status,
                pdf,
                cdf,
                meanReg,
                ngrid,
                npred,
                hpd,
                bci,
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
                kappa,
                grid,
                xpred
    )
  }
  
  cat("Finished!", "\n")
  
  
  #---------------------------------------------- 
  # returns
  #---------------------------------------------- 
  res$proc.time = proc.time() - ptm
  
  if(any(pdf, cdf, meanReg)) {
    res$prediction = TRUE
    res$type.pred = type.pred
    res$compute.band = compute.band
    res$type.band = type.band
    res$xpred = xpred
    res$grid = grid
  } else {
    res$prediction = FALSE
  }
  
  attr(res, 'class') <- 'DPMcdensity'
  
  return(res)
}