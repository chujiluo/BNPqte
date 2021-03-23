predict.DPMcdensity = function(object, xpred, grid, 
                               type.pred=c("pdf"), compute.band=TRUE, type.band="HPD",
                               printevery=100L) {
  
  #---------------------------------------------- 
  # process arguments
  #---------------------------------------------- 
  d = length(object$state$m)
  
  ## xpred
  xpred = as.matrix(xpred)
  npred = nrow(xpred)
  if(ncol(xpred) != (d-1))
    stop(paste("xpred is required to have the same number of columns as x, i.e. ncol(xpred) = ", d-1, ".", sep = ""))
  
  ## grid
  grid = as.vector(grid)
  ngrid = length(grid)
  
  ## type.pred
  pdf = cdf = meanReg = FALSE
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
  
  ## compute.band, type.band
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
  
  ## extract posterior samples from object
  method = object$method
  updateAlpha = object$updateAlpha
  useHyperpriors = object$useHyperpriors
  
  ZetaList = object$posterior$Zeta # ndpost, nclusterxd, each row for cluster
  OmegaList = object$posterior$Omega # ndpost - nclusters - dxd
  ndpost = length(ZetaList)
  
  if(method == "truncated") {
    lwList = object$posterior$lw
    nclusters = ncol(lwList)
  } else {
    kappaList = object$posterior$kappa # ndpost x n
    nclusterList = object$posterior$nclusters  # ndpost
    
    alphaList = as.vector(object$posterior$alpha)  # ndpost or 1
    lambdaList = as.vector(object$posterior$lambda)  # ndpost or 1
    
    if(useHyperpriors) {
      mList = object$posterior$m  # ndpost x d
      PsiList = object$posterior$Psi   # ndpost, d x d
    } else {
      mList = t(as.matrix(object$posterior$m))  # d x 1
      PsiList = list(object$posterior$Ps)  # 1, d x d
    }
    
    nu = object$state$nu
    n = ncol(kappaList)
  }
  
  #---------------------------------------------- 
  # call cpp function
  #---------------------------------------------- 
  cat("Predicting ", sep = "")
  if(method == "truncated") {
    res = .Call("_BNPqte_cpDPMcdensity",
                ngrid,
                npred,
                grid,
                xpred,
                d,
                nclusters,
                ndpost,
                ZetaList,
                OmegaList,
                lwList,
                pdf,
                cdf,
                meanReg,
                hpd,
                bci,
                printevery
    )
  } else {
    res = .Call("_BNPqte_cpDPMcdensityNeal",
                ngrid,
                npred,
                grid,
                xpred,
                n,
                d,
                updateAlpha,
                useHyperpriors,
                ZetaList,
                OmegaList,
                kappaList,
                nclusterList,
                alphaList,
                lambdaList,
                mList,
                PsiList,
                nu,
                ndpost,
                printevery,
                pdf,
                cdf,
                meanReg,
                hpd,
                bci
    )
  }
  
  cat("Finished!\n")
  
  
  #---------------------------------------------- 
  # returns
  #---------------------------------------------- 
  res$prediction = TRUE
  res$type.pred = type.pred
  res$compute.band = compute.band
  res$xpred = xpred
  res$grid = grid
  
  attr(res, 'class') <- 'DPMcdensity'
  
  return(res)
}
