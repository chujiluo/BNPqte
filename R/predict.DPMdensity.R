predict.DPMdensity = function(object, grid) {

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
  ZetaList = object$posterior$Zeta # ndpost, nclusterxd, each row for cluster
  OmegaList = object$posterior$Omega # ndpost - nclusters - dxd
  lwList = object$posterior$lw
  nclusters = ncol(lwList)
  ndpost = nrow(lwList)

  
  #---------------------------------------------- 
  # call cpp function
  #---------------------------------------------- 
  res = .Call("_BNPqte_cpDPMdensity",
              ngrid,
              grid1,
              grid2,
              d,
              nclusters,
              ndpost,
              ZetaList,
              OmegaList,
              lwList
  )
  
  
  #---------------------------------------------- 
  # returns
  #---------------------------------------------- 
  
  attr(res, 'class') <- 'DPMdensity'
  
  return(res)
}
