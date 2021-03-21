predict.DPMdensity = function(object, grid, printevery=100L) {
  
  #---------------------------------------------- 
  # process arguments
  #---------------------------------------------- 
  if (is.matrix(grid) & (ncol(grid)==2)) {
    ngrid = nrow(grid)
    grid1 = grid[, 1]
    grid2 = grid[, 2]
  } else {
    stop("grid is required to be a matrix with 2 columns.")
  }
  
  
  d = 2
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
    res = .Call("_BNPqte_cpDPMdensity",
                ngrid,
                grid1,
                grid2,
                d,
                nclusters,
                ndpost,
                ZetaList,
                OmegaList,
                lwList,
                printevery
    )
  } else {
    res = .Call("_BNPqte_cpDPMdensityNeal",
                ngrid,
                grid1,
                grid2,
                d,
                n,
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
                printevery)
  }
  
  cat("Finished!\n")
  
  #---------------------------------------------- 
  # returns
  #---------------------------------------------- 
  
  attr(res, 'class') <- 'DPMdensity'
  
  return(res)
}