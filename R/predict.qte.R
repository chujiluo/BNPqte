predict.qte = function(object, probs, 
                       compute.band=TRUE, type.band="HPD", alphas=c(0.05)) {
  
  # process args
  if(!("cdf" %in% object$type.pred))
    stop("cdf prediction is not implemented in object, 
         so no quantiles or quantiles treatment effects can be predicted.")
  
  
  probs = sort(as.vector(probs))
  if(!all((probs>0) & (probs<1)))
    stop("probs are required to be within (0, 1).")
  nprobs = length(probs)

  
  if(compute.band) {
    if(!all(type.band %in% c("HPD", "BCI"))) {
      stop("Only two available bands are provided: HPD (Highest posterior interval ) 
           and BCI (Bayesian credible interval).")
    }
    else {
      if(!all((alphas > 0) & (alphas < 1)))
        stop("alphas are required to be a vector of elements between 0 and 1.")
      alphas = as.vector(alphas)
      nalphas = length(alphas)
    }
  }
  
  
  # read in cdfs
  group0.cdfs = object$control.cdfs     # (bart.ndpost*dpm.ndpost) x dpm.ngrid
  group1.cdfs = object$treatment.cdfs
  
  grids = object$grid
  ngrid = length(object$grid)
  ndpost = nrow(group0.cdfs)
  
  
  # calculate quantiles
  quantiles0 = .Call(`_BNPqte_quantile_fun`, ngrid, nprobs, ndpost, grids, probs, group0.cdfs)
  quantiles1 = .Call(`_BNPqte_quantile_fun`, ngrid, nprobs, ndpost, grids, probs, group1.cdfs)
  
  
  # quantile treatment effects
  qtes = quantiles1 - quantiles0
  
  if(compute.band) {
    lower.band = matrix(NA, nrow = nprobs, ncol = nalphas)
    upper.band = matrix(NA, nrow = nprobs, ncol = nalphas)
    
    for (i in 1:nprobs) {
      band = .Call(`_BNPqte_credible_interval`, nrow(qtes), qtes[, i], alphas, type.band)
      lower.band[i, ] = band[1, ]
      upper.band[i, ] = band[2, ]
    }
  }
  
  
  # returns
  res = list()
  res$probs = probs
  res$compute.band = compute.band
  if(compute.band) {
    res$type.band = type.band
    res$alphas = alphas
  }
  
  res$control.quantiles = quantiles0
  res$treatment.quantiles = quantiles1
  res$quantiles.avg = cbind(colMeans(quantiles0), colMeans(quantiles1))
  colnames(res$control.quantiles) = paste(probs*100, "%", sep = "")
  colnames(res$treatment.quantiles) = paste(probs*100, "%", sep = "")
  colnames(res$quantiles.avg) = c("control", "treatment")
  rownames(res$quantiles.avg) = paste(probs*100, "%", sep = "")
  
  res$qtes = qtes
  res$qtes.avg = colMeans(qtes)
  colnames(res$qtes) = paste(probs*100, "%", sep = "")
  names(res$qtes.avg) = paste(probs*100, "%", sep = "")
  
  if(compute.band) {
    res$lower.band = lower.band
    res$upper.band = upper.band
    rownames(res$lower.band) = paste(probs*100, "%", sep = "")
    rownames(res$upper.band) = paste(probs*100, "%", sep = "")
    colnames(res$lower.band) = alphas
    colnames(res$upper.band) = alphas
  }
  
  attr(res, 'class') <- 'qte'
  
  return(res)
}