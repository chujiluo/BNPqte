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
    if(nalphas == 1) {
      control.quantiles.ci = matrix(NA, nrow = nprobs, ncol = 2)
      for (i in 1:nprobs) {
        band = .Call(`_BNPqte_credible_interval`, nrow(quantiles0), quantiles0[, i], alphas, type.band)
        control.quantiles.ci[i, 1] = band[1, ]
        control.quantiles.ci[i, 2] = band[2, ]
      }
      rownames(control.quantiles.ci) = paste(probs*100, "%", sep = "")
      colnames(control.quantiles.ci) = c(alphas, 1-alphas)
      
      treatment.quantiles.ci = matrix(NA, nrow = nprobs, ncol = 2)
      for (i in 1:nprobs) {
        band = .Call(`_BNPqte_credible_interval`, nrow(quantiles1), quantiles1[, i], alphas, type.band)
        treatment.quantiles.ci[i, 1] = band[1, ]
        treatment.quantiles.ci[i, 2] = band[2, ]
      }
      rownames(treatment.quantiles.ci) = paste(probs*100, "%", sep = "")
      colnames(treatment.quantiles.ci) = c(alphas, 1-alphas)
      
      qtes.ci = matrix(NA, nrow = nprobs, ncol = 2)
      for (i in 1:nprobs) {
        band = .Call(`_BNPqte_credible_interval`, nrow(qtes), qtes[, i], alphas, type.band)
        qtes.ci[i, 1] = band[1, ]
        qtes.ci[i, 2] = band[2, ]
      }
      rownames(qtes.ci) = paste(probs*100, "%", sep = "")
      colnames(qtes.ci) = c(alphas, 1-alphas)
    } else {
      qtes.ci = list()
      control.quantiles.ci = list()
      treatment.quantiles.ci = list()
      
      for (i in 1:nalphas) {
        qtes.ci[[i]] = matrix(NA, nrow = nprobs, ncol = 2)
        control.quantiles.ci[[i]] = matrix(NA, nrow = nprobs, ncol = 2)
        treatment.quantiles.ci[[i]] = matrix(NA, nrow = nprobs, ncol = 2)
        
        rownames(qtes.ci[[i]]) = paste(probs*100, "%", sep = "")
        rownames(control.quantiles.ci[[i]]) = paste(probs*100, "%", sep = "")
        rownames(treatment.quantiles.ci[[i]]) = paste(probs*100, "%", sep = "")
        
        colnames(qtes.ci[[i]]) = c(alphas[i], 1-alphas[i])
        colnames(control.quantiles.ci[[i]]) = c(alphas[i], 1-alphas[i])
        colnames(treatment.quantiles.ci[[i]]) = c(alphas[i], 1-alphas[i])
      }
      
      for (i in 1:nprobs) {
        band1 = .Call(`_BNPqte_credible_interval`, nrow(qtes), qtes[, i], alphas, type.band)
        band2 = .Call(`_BNPqte_credible_interval`, nrow(quantiles0), quantiles0[, i], alphas, type.band)
        band3 = .Call(`_BNPqte_credible_interval`, nrow(quantiles1), quantiles1[, i], alphas, type.band)
        
        for (j in 1:nalphas) {
          qtes.ci[[j]][i, 1] = band1[1, j]
          qtes.ci[[j]][i, 2] = band1[2, j]
          
          control.quantiles.ci[[j]][i, 1] = band2[1, j]
          control.quantiles.ci[[j]][i, 2] = band2[2, j]
          
          treatment.quantiles.ci[[j]][i, 1] = band3[1, j]
          treatment.quantiles.ci[[j]][i, 2] = band3[2, j]
        }
      }
    }
  }
  
  
  # returns
  res = list()
  res$type.pred = "cdf"
  res$probs = probs
  res$compute.band = compute.band
  if(compute.band) {
    res$type.band = type.band
    res$alphas = alphas
  }
  
  res$control.quantiles = quantiles0
  res$treatment.quantiles = quantiles1
  colnames(res$control.quantiles) = paste(probs*100, "%", sep = "")
  colnames(res$treatment.quantiles) = paste(probs*100, "%", sep = "")
  
  res$control.quantiles.avg = colMeans(quantiles0)
  res$treatment.quantiles.avg = colMeans(quantiles1)
  names(res$control.quantiles.avg) = paste(probs*100, "%", sep = "")
  names(res$treatment.quantiles.avg) = paste(probs*100, "%", sep = "")
  
  res$qtes = qtes
  res$qtes.avg = colMeans(qtes)
  colnames(res$qtes) = paste(probs*100, "%", sep = "")
  names(res$qtes.avg) = paste(probs*100, "%", sep = "")
  
  if(compute.band) {
    res$qtes.ci = qtes.ci
    res$control.quantiles.ci = control.quantiles.ci
    res$treatment.quantiles.ci = treatment.quantiles.ci
  }
  
  attr(res, 'class') <- 'qte'
  
  return(res)
}