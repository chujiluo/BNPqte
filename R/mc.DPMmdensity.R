# Estimate marginal density, cdf and its quantiles by fitting a DPMM on the joint data and 
# using Bayesian boostrap to approximate the distribution of propensity scores;

# Do the above task for multiple data sets in parallel
## The posterior of DPMM is sampled using blocked Gibbs sampling
mc.DPMmdensity = function(n0, n1, d, y0, y1, propensity0, propensity1, 
                          diag, pdf, cdf, ngrid, grid, npred, xpred,
                          updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu, 
                          nclusters, nskip, ndpost, keepevery, 
                          diri, probs, nprobs) {
  
  res = list()
  
  if(is.matrix(propensity0)) {
    for(i in 1:nrow(propensity0)) {
      tmp0 = .Call("_BNPqte_cDPMmdensity",
                   n0, d, cbind(y0, propensity0[i, ]), y0, as.matrix(propensity0[i, ]),
                   diag, pdf, cdf, ngrid, grid, npred, as.matrix(xpred[i, ]),
                   updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu, 
                   nclusters, nskip, ndpost, keepevery,
                   diri[i, ], probs, nprobs)
      
      tmp1 = .Call("_BNPqte_cDPMmdensity",
                   n1, d, cbind(y1, propensity1[i, ]), y1, as.matrix(propensity1[i, ]),
                   diag, pdf, cdf, ngrid, grid, npred, as.matrix(xpred[i, ]),
                   updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu,
                   nclusters, nskip, ndpost, keepevery,
                   diri[i, ], probs, nprobs)
      
      if(i == 1) {
        if(pdf) {
          group0.pdfs = tmp0$predict.pdfs
          group1.pdfs = tmp1$predict.pdfs
        }
        if(cdf) {
          group0.cdfs = tmp0$predict.cdfs
          group1.cdfs = tmp1$predict.cdfs
          quantiles0 = tmp0$predict.quantiles
          quantiles1 = tmp1$predict.quantiles
          qtes = tmp1$predict.quantiles - tmp0$predict.quantiles
        }
        if(diag) {
          group0.lmpps = tmp0$logMPPs
          group1.lmpps = tmp1$logMPPs
          
          group0.ylogliks = tmp0$ylogliks
          group1.ylogliks = tmp1$ylogliks
        }
      } else {
        if(pdf) {
          group0.pdfs = rbind(group0.pdfs, tmp0$predict.pdfs)
          group1.pdfs = rbind(group1.pdfs, tmp1$predict.pdfs)
        }
        if(cdf) {
          group0.cdfs = rbind(group0.cdfs, tmp0$predict.cdfs)
          group1.cdfs = rbind(group1.cdfs, tmp1$predict.cdfs)
          quantiles0 = rbind(quantiles0, tmp0$predict.quantiles)
          quantiles1 = rbind(quantiles1, tmp1$predict.quantiles)
          qtes = rbind(qtes, tmp1$predict.quantiles - tmp0$predict.quantiles)
        }
        if(diag) {
          group0.lmpps = rbind(group0.lmpps, tmp0$logMPPs)
          group1.lmpps = rbind(group1.lmpps, tmp1$logMPPs)
          
          group0.ylogliks = rbind(group0.ylogliks, tmp0$ylogliks)
          group1.ylogliks = rbind(group1.ylogliks, tmp1$ylogliks)
        }
      }
      
      
      rm(tmp0)
      rm(tmp1)
      gc()
    }
    
    if(pdf) {
      res$group0.pdfs = group0.pdfs
      res$group1.pdfs = group1.pdfs
    }
    if(cdf) {
      res$group0.cdfs = group0.cdfs
      res$group1.cdfs = group1.cdfs
      res$quantiles0 = quantiles0
      res$quantiles1 = quantiles1
      res$qtes = qtes
    }
    if(diag) {
      res$group0.lmpps = group0.lmpps
      res$group1.lmpps = group1.lmpps
      
      res$group0.ylogliks = group0.ylogliks
      res$group1.ylogliks = group1.ylogliks
    }
  } else {
    ## only one row in propensity0, propensity1, xpred and diri
    tmp0 = .Call("_BNPqte_cDPMmdensity",
                 n0, d, cbind(y0, propensity0), y0, as.matrix(propensity0),
                 diag, pdf, cdf, ngrid, grid, npred, as.matrix(xpred),
                 updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu, 
                 nclusters, nskip, ndpost, keepevery,
                 diri, probs, nprobs)
    
    tmp1 = .Call("_BNPqte_cDPMmdensity",
                 n1, d, cbind(y1, propensity1), y1, as.matrix(propensity1),
                 diag, pdf, cdf, ngrid, grid, npred, as.matrix(xpred),
                 updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu,
                 nclusters, nskip, ndpost, keepevery,
                 diri, probs, nprobs)
    
    if(pdf) {
      res$group0.pdfs = tmp0$predict.pdfs
      res$group1.pdfs = tmp1$predict.pdfs
    }
    if(cdf) {
      res$group0.cdfs = tmp0$predict.cdfs
      res$group1.cdfs = tmp1$predict.cdfs
      res$quantiles0 = tmp0$predict.quantiles
      res$quantiles1 = tmp1$predict.quantiles
      res$qtes = tmp1$predict.quantiles - tmp0$predict.quantiles
    }
    if(diag) {
      res$group0.lmpps = tmp0$logMPPs
      res$group1.lmpps = tmp1$logMPPs
      
      res$group0.ylogliks = tmp0$ylogliks
      res$group1.ylogliks = tmp1$ylogliks
    }
    
    rm(tmp0)
    rm(tmp1)
    gc()
  }
  
  return(res)
}


# Do the above task for multiple data sets in parallel
## The posterior of DPMM is sampled using Algorithm 8 with m = 1 in Neal, 2000
mc.DPMmdensityNeal = function(n0, n1, d, y0, y1, propensity0, propensity1, 
                              diag, pdf, cdf, ngrid, grid, npred, xpred,
                              updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu, 
                              nclusters, nskip, ndpost, keepevery, 
                              diri, probs, nprobs) {
  
  res = list()
  
  if(is.matrix(propensity0)) {
    for(i in 1:nrow(propensity0)) {
      tmp0 = .Call("_BNPqte_cDPMmdensityNeal",
                   n0, d, cbind(y0, propensity0[i, ]), y0, as.matrix(propensity0[i, ]),
                   diag, pdf, cdf, ngrid, grid, npred, as.matrix(xpred[i, ]),
                   updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu, 
                   nclusters, nskip, ndpost, keepevery, diri[i, ], probs, nprobs)
      
      tmp1 = .Call("_BNPqte_cDPMmdensityNeal",
                   n1, d, cbind(y1, propensity1[i, ]), y1, as.matrix(propensity1[i, ]),
                   diag, pdf, cdf, ngrid, grid, npred, as.matrix(xpred[i, ]),
                   updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu, 
                   nclusters, nskip, ndpost, keepevery, diri[i, ], probs, nprobs)
      
      if(i == 1) {
        if(pdf) {
          group0.pdfs = tmp0$predict.pdfs
          group1.pdfs = tmp1$predict.pdfs
        }
        if(cdf) {
          group0.cdfs = tmp0$predict.cdfs
          group1.cdfs = tmp1$predict.cdfs
          quantiles0 = tmp0$predict.quantiles
          quantiles1 = tmp1$predict.quantiles
          qtes = tmp1$predict.quantiles - tmp0$predict.quantiles
        }
        if(diag) {
          group0.ylogliks = tmp0$ylogliks
          group1.ylogliks = tmp1$ylogliks
        }
      } else {
        if(pdf) {
          group0.pdfs = rbind(group0.pdfs, tmp0$predict.pdfs)
          group1.pdfs = rbind(group1.pdfs, tmp1$predict.pdfs)
        }
        if(cdf) {
          group0.cdfs = rbind(group0.cdfs, tmp0$predict.cdfs)
          group1.cdfs = rbind(group1.cdfs, tmp1$predict.cdfs)
          quantiles0 = rbind(quantiles0, tmp0$predict.quantiles)
          quantiles1 = rbind(quantiles1, tmp1$predict.quantiles)
          qtes = rbind(qtes, tmp1$predict.quantiles - tmp0$predict.quantiles)
        }
        if(diag) {
          group0.ylogliks = rbind(group0.ylogliks, tmp0$ylogliks)
          group1.ylogliks = rbind(group1.ylogliks, tmp1$ylogliks)
        }
      }
      
      
      rm(tmp0)
      rm(tmp1)
      gc()
    }
    
    if(pdf) {
      res$group0.pdfs = group0.pdfs
      res$group1.pdfs = group1.pdfs
    }
    if(cdf) {
      res$group0.cdfs = group0.cdfs
      res$group1.cdfs = group1.cdfs
      res$quantiles0 = quantiles0
      res$quantiles1 = quantiles1
      res$qtes = qtes
    }
    if(diag) {
      res$group0.ylogliks = group0.ylogliks
      res$group1.ylogliks = group1.ylogliks
    }
    
  } else {
    ## only one row in propensity0, propensity1, xpred and diri
    tmp0 = .Call("_BNPqte_cDPMmdensityNeal",
                 n0, d, cbind(y0, propensity0), y0, as.matrix(propensity0),
                 diag, pdf, cdf, ngrid, grid, npred, as.matrix(xpred),
                 updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu, 
                 nclusters, nskip, ndpost, keepevery,
                 diri, probs, nprobs)
    
    tmp1 = .Call("_BNPqte_cDPMmdensityNeal",
                 n1, d, cbind(y1, propensity1), y1, as.matrix(propensity1),
                 diag, pdf, cdf, ngrid, grid, npred, as.matrix(xpred),
                 updateAlpha, useHyperpriors, alpha, a0, b0, lambda, gamma1, gamma2, nu0, nu,
                 nclusters, nskip, ndpost, keepevery,
                 diri, probs, nprobs)
    
    if(pdf) {
      res$group0.pdfs = tmp0$predict.pdfs
      res$group1.pdfs = tmp1$predict.pdfs
    }
    if(cdf) {
      res$group0.cdfs = tmp0$predict.cdfs
      res$group1.cdfs = tmp1$predict.cdfs
      res$quantiles0 = tmp0$predict.quantiles
      res$quantiles1 = tmp1$predict.quantiles
      res$qtes = tmp1$predict.quantiles - tmp0$predict.quantiles
    }
    if(diag) {
      res$group0.ylogliks = tmp0$ylogliks
      res$group1.ylogliks = tmp1$ylogliks
    }
    
    rm(tmp0)
    rm(tmp1)
    gc()
  }
  
  return(res)
}
