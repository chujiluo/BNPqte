qte = function(y, x, treatment,
               probs=c(0.1, 0.25, 0.50, 0.75, 0.90), 
               compute.band=TRUE, type.band="HPD", alphas=c(0.05),
               bart.link="probit", bart.params=list(split.prob="polynomial"), 
               dpm.params=list(),
               mc.cores=1L, seed = 123) {
  
  res = list()
  
  #---------------------------------------------------------------------------------------------------- 
  # check and process arguments
  #---------------------------------------------------------------------------------------------------- 
  ##---------------------------------------------------------------------------
  ## y: continuous outcome
  ##---------------------------------------------------------------------------
  y = as.vector(y)
  if(!is.numeric(y))
    stop("y is required to be a numeric vector.")
  n = length(y)
  
  ##---------------------------------------------------------------------------
  ## x: mixed-type covariates (NOTE: categorical covariates need to be factor.)
  ##---------------------------------------------------------------------------
  if(nrow(x) != n)
    stop("x is a matrix or data frame with number of rows equal to length(y).")
  
  ##---------------------------------------------------------------------------
  ## treatment: binary
  ##---------------------------------------------------------------------------
  treatment = as.vector(treatment)
  if(length(treatment) != n)
    stop("treatment is a vector with length equal to length(y).")
  if(!all(unique(treatment) %in% c(0, 1)))
    stop("treatment is required to be a vector consists of 0 and 1.")
  
  d = 2
  group0 = which(treatment == 0)
  group1 = which(treatment == 1)
  n0 = length(group0)
  n1 = length(group1)
  y0 = y[group0]
  y1 = y[group1]
  
  ##---------------------------------------------------------------------------
  ## probs
  ##---------------------------------------------------------------------------
  probs = sort(as.vector(probs))
  if(!all((probs>0) & (probs<1)))
    stop("probs are required to be within (0, 1).")
  nprobs = length(probs)
  
  ##---------------------------------------------------------------------------
  ## compute.band, type.band, alphas
  ##---------------------------------------------------------------------------
  if(compute.band) {
    if(!all(type.band %in% c("HPD", "BCI"))) {
      stop("Only two available bands are provided: HPD (Highest posterior interval ) and BCI (Bayesian credible interval).")
    }
    else {
      if(!all((alphas > 0) & (alphas < 1)))
        stop("alphas are required to be a vector of elements between 0 and 1.")
      alphas = as.vector(alphas)
      nalphas = length(alphas)
    }
  }
  
  ##---------------------------------------------------------------------------
  ## bart.link
  ##---------------------------------------------------------------------------
  if(!(bart.link %in% c("probit", "logit")))
    stop("Only two available links: probit or logit, are provided for fitting a BART model for binary responses.")
  
  ##---------------------------------------------------------------------------
  ## bart.params (currently for BART::pbart and BART::lbart)
  ## if nothing is specified, then default values will be used
  ##---------------------------------------------------------------------------
  if(is.null(bart.params$split.prob)) {
    bart.split.prob = "polynomial"
  } else {
    bart.split.prob = bart.params$split.prob
  }
  if(!(bart.split.prob %in% c("polynomial", "exponential")))
    stop("Only two types of available splitting probability for BART: polynomial or exponential.")
  if(is.null(bart.params$sparse)) {
    bart.sparse = FALSE
  } else {
    bart.sparse = bart.params$sparse
  }
  if(is.null(bart.params$theta)) {
    bart.theta = 0
  } else {
    bart.theta = bart.params$theta
  }
  if(is.null(bart.params$omega)) {
    bart.omega = 1
  } else {
    bart.omega = bart.params$omega
  }
  if(is.null(bart.params$a)) {
    bart.a = 0.5
  } else {
    bart.a = bart.params$a
  }
  if(is.null(bart.params$b)) {
    bart.b = 1
  } else {
    bart.b = bart.params$b
  }
  if(is.null(bart.params$augment)) {
    bart.augment = FALSE
  } else {
    bart.augment = bart.params$augment
  }
  if(is.null(bart.params$rho)) {
    bart.rho = NULL
  } else {
    bart.rho = bart.params$rho
  }
  if(is.null(bart.params$xinfo)) {
    bart.xinfo = matrix(0.0,0,0)
  } else {
    bart.xinfo = bart.params$xinfo
  }
  if(is.null(bart.params$usequants)) {
    bart.usequants = FALSE
  } else {
    bart.usequants = bart.params$usequants
  }
  if(is.null(bart.params$cont)) {
    bart.cont = FALSE
  } else {
    bart.cont = bart.params$cont
  }
  if(is.null(bart.params$rm.const)) {
    bart.rm.const = TRUE
  } else {
    bart.rm.const = bart.params$rm.const
  }
  if(is.null(bart.params$k)) {
    bart.k = 2.0
  } else {
    bart.k = bart.params$k
  }
  if(is.null(bart.params$power)) {
    bart.power = 2.0
  } else {
    bart.power = bart.params$power
  }
  if(is.null(bart.params$base)) {
    bart.base = 0.95
  } else {
    bart.base = bart.params$base
  }
  if(is.null(bart.params$binaryOffset)) {
    bart.binaryOffset = NULL
  } else {
    bart.binaryOffset = bart.params$binaryOffset
  }
  if(is.null(bart.params$ntree)) {
    bart.ntree = 50L
  } else {
    bart.ntree = bart.params$ntree
  }
  if(is.null(bart.params$numcut)) {
    bart.numcut = 100L
  } else {
    bart.numcut = bart.params$numcut
  }
  if(is.null(bart.params$ndpost)) {
    bart.ndpost = 5L
  } else {
    bart.ndpost = bart.params$ndpost
  }
  if(is.null(bart.params$nskip)) {
    bart.nskip = 500L
  } else {
    bart.nskip = bart.params$nskip
  }
  if(is.null(bart.params$keepevery)) {
    bart.keepevery = 100L
  } else {
    bart.keepevery = bart.params$keepevery
  }
  if(is.null(bart.params$nkeeptrain)) {
    bart.nkeeptrain = bart.ndpost
  } else {
    bart.nkeeptrain = bart.params$nkeeptrain
  }
  if(is.null(bart.params$nkeeptest)) {
    bart.nkeeptest = bart.ndpost
  } else {
    bart.nkeeptest = bart.params$nkeeptest
  }
  if(is.null(bart.params$nkeeptreedraws)) {
    bart.nkeeptreedraws = bart.ndpost
  } else {
    bart.nkeeptreedraws = bart.params$nkeeptreedraws
  }
  if(is.null(bart.params$printevery)) {
    bart.printevery = 100L
  } else {
    bart.printevery = bart.params$printevery
  }
  if(is.null(bart.params$transposed)) {
    bart.transposed = FALSE
  } else {
    bart.transposed = bart.params$transposed
  }
  if(is.null(bart.params$tau.interval)) {
    bart.tau.interval = 0.95
  } else {
    bart.tau.interval = bart.params$tau.interval
  }
  
  
  #---------------------------------------------------------------------------------------------------- 
  # set seed
  #---------------------------------------------------------------------------------------------------- 
  set.seed(seed)

  
  #---------------------------------------------------------------------------------------------------- 
  # fit probit or logit BART (Will implement two splitting probs in the future)
  #---------------------------------------------------------------------------------------------------- 
  cat("---------------------------------\n")
  cat(" Modeling Treatment ~ Covariates\n")
  cat("---------------------------------\n")
  
  if(bart.link == "probit") {
    ## fit probit bart
    bart_obj = pbart(x.train=x, y.train=treatment,
                     sparse=bart.sparse, theta=bart.theta, omega=bart.omega,
                     a=bart.a, b=bart.b, augment=bart.augment, rho=bart.rho,
                     xinfo=bart.xinfo, usequants=bart.usequants,
                     cont=bart.cont, rm.const=bart.rm.const,
                     k=bart.k, power=bart.power, base=bart.base,
                     binaryOffset=bart.binaryOffset,
                     ntree=bart.ntree, numcut=bart.numcut,
                     ndpost=bart.ndpost, nskip=bart.nskip, keepevery=bart.keepevery,
                     nkeeptrain=bart.nkeeptrain, nkeeptest=bart.nkeeptest,
                     nkeeptreedraws=bart.nkeeptreedraws,
                     printevery=bart.printevery, transposed=bart.transposed)
  } else {
    ## fit logit bart
    bart_obj = lbart(x.train=x, y.train=treatment,
                     sparse=bart.sparse, a=bart.a, b=bart.b, augment=bart.augment, rho=bart.rho,
                     xinfo=bart.xinfo, usequants=bart.usequants,
                     cont=bart.cont, rm.const=bart.rm.const, tau.interval=bart.tau.interval,
                     k=bart.k, power=bart.power, base=bart.base,
                     binaryOffset=bart.binaryOffset,
                     ntree=bart.ntree, numcut=bart.numcut,
                     ndpost=bart.ndpost, nskip=bart.nskip, keepevery=bart.keepevery,
                     nkeeptrain=bart.nkeeptrain, nkeeptest=bart.nkeeptest,
                     nkeeptreedraws=bart.nkeeptreedraws,
                     printevery=bart.printevery, transposed=bart.transposed)
  }
  
  ## estimated propensity scores: bart.ndpost x n
  propensity = bart_obj$prob.train
  propensity0 = propensity[, group0, drop = FALSE]
  propensity1 = propensity[, group1, drop = FALSE]
  
  res$propensity = propensity
  res$bart.parmas = list(link = bart.link, ntree = bart.ntree, 
                         nskip = bart.nskip, ndpost = bart.ndpost, keepevery = bart.keepevery,
                         split.prob = bart.split.prob, sparse = bart.sparse,
                         binaryOffset = bart_obj$binaryOffset, 
                         vip = bart_obj$varcount.mean, 
                         varprob.mean = bart_obj$varprob.mean)
  
  rm(bart_obj)
  gc()
  
  #---------------------------------------------------------------------------------------------------- 
  # fit DPMM
  #---------------------------------------------------------------------------------------------------- 
   
  ##---------------------------------------------------------------------------
  ## process dpm.params: if not specified by users, then use default values
  ##---------------------------------------------------------------------------
  if(is.null(dpm.params$nclusters)) {
    dpm.nclusters = 50L
  } else {
    dpm.nclusters = dpm.params$nclusters
  }
  if(is.null(dpm.params$ngrid)) {
    dpm.ngrid = 100L
  } else {
    dpm.ngrid = dpm.params$ngrid
  }
  if(is.null(dpm.params$grid)) {
    dpm.grid = NULL
  } else {
    dpm.grid = dpm.params$grid
  }
  if(is.null(dpm.params$type.pred)) {
    dpm.type.pred = c("cdf")
  } else {
    dpm.type.pred = dpm.params$type.pred
  }
  dpm.pdf = dpm.cdf = FALSE
  if((dpm.ngrid <= 0) & is.null(dpm.grid)) {
    stop("Either positive ngrid or non-empty grid is required to be provided in dpm.params.")
  } else {
    dpm.type.pred = as.vector(dpm.type.pred)
    if(!all(dpm.type.pred %in% c("pdf", "cdf"))) {
      stop("Only two types of available type.pred: pdf and cdf.")
    } else {
      if("pdf" %in% dpm.type.pred)
        dpm.pdf = TRUE
      if("cdf" %in% dpm.type.pred)
        dpm.cdf = TRUE
    }
    
    if(is.null(dpm.grid)) {
      # common grids for two groups
      dpm.grid = seq(from = (min(y) - 0.25 * sd(y)), to = (max(y) + 0.25 * sd(y)), length.out = dpm.ngrid)
    } else {
      dpm.grid = sort(as.vector(dpm.grid))
      dpm.ngrid = length(dpm.grid)
    }
  }
  
  if(is.null(dpm.params$xpred)) {
    dpm.xpred = NULL
  } else {
    dpm.xpred = dpm.params$xpred
  }
  if(is.null(dpm.params$npred)) {
    dpm.npred = 0
  } else {
    dpm.npred = dpm.params$npred
  }
  # as.matrix(propensity[i, ])
  if(is.null(dpm.xpred) & (dpm.npred == 0)) {
    dpm.xpred = propensity
    dpm.npred = ncol(propensity)
  } else {
    if(is.null(dpm.xpred)) {
      dpm.xpred = seq(from = (min(propensity) - 0.25 * sd(propensity)), 
                      to = (max(propensity) + 0.25 * sd(propensity)), 
                      length.out = dpm.npred)
      dpm.xpred = t(replicate(bart.ndpost, dpm.xpred))
    } else {
      if((!is.vector(dpm.xpred)) | (!is.matrix(dpm.xpred))) {
        stop("xpred in dpm.params is required to be a vector, matrix or NULL.")
      } else {
        if(is.vector(dpm.xpred)) {
          dpm.npred = length(dpm.xpred)
          dpm.xpred = t(replicate(bart.ndpost, dpm.xpred))
        } else {
          if(nrow(dpm.xpred) != bart.ndpost) {
            stop("When xpred in dpm.params is provided as a matrix, 
                 it is required to have the same number of rows as ndpost in bart.params.")
          } else {
            dpm.npred = ncol(dpm.xpred)
          }
        }
      }
    }
  }
  
  if(is.null(dpm.params$diag)) {
    dpm.diag = FALSE
  } else {
    dpm.diag = dpm.params$diag
  }
  if(is.null(dpm.params$nskip)) {
    dpm.nskip = 500L
  } else {
    dpm.nskip = dpm.params$nskip
  }
  if(is.null(dpm.params$ndpost)) {
    dpm.ndpost = 200L
  } else {
    dpm.ndpost = dpm.params$ndpost
  }
  if(is.null(dpm.params$keepevery)) {
    dpm.keepevery = 2L
  } else {
    dpm.keepevery = dpm.params$keepevery
  }
  
  if(is.null(dpm.params$updateAlpha)) {
    dpm.updateAlpha = TRUE
  } else {
    dpm.updateAlpha = dpm.params$updateAlpha
  }
  if(is.null(dpm.params$useHyperpriors)) {
    dpm.useHyperpriors = TRUE
  } else {
    dpm.useHyperpriors = dpm.params$useHyperpriors
  }
  if(is.null(dpm.params$alpha)) {
    dpm.alpha = 10.0
  } else {
    dpm.alpha = dpm.params$alpha
  }
  if(is.null(dpm.params$a0)) {
    dpm.a0 = 10.0
  } else {
    dpm.a0 = dpm.params$a0
  }
  if(is.null(dpm.params$b0)) {
    dpm.b0 = 1.0
  } else {
    dpm.b0 = dpm.params$b0
  }
  if(is.null(dpm.params$lambda)) {
    dpm.lambda = 0.5
  } else {
    dpm.lambda = dpm.params$lambda
  }
  if(is.null(dpm.params$gamma1)) {
    dpm.gamma1 = 3.0
  } else {
    dpm.gamma1 = dpm.params$gamma1
  }
  if(is.null(dpm.params$gamma2)) {
    dpm.gamma2 = 2.0
  } else {
    dpm.gamma2 = dpm.params$gamma2
  }
  if(is.null(dpm.params$nu)) {
    dpm.nu = 4
  } else {
    dpm.nu = dpm.params$nu
  }
  if(is.null(dpm.params$nu0)) {
    dpm.nu0 = 4
  } else {
    dpm.nu0 = dpm.params$nu0
  }
  if (dpm.updateAlpha) {
    if ((dpm.a0 > 0) & (dpm.b0 > 0)) 
      dpm.alpha = -1
    else
      stop("a0 and b0 in dpm.params are required to be positive scalars.")
  } else {
    if (dpm.alpha > 0)
      dpm.a0 = dpm.b0 = -1
    else
      stop("alpha in dpm.params is required to be a positive scalar.")
  }
  if (dpm.nu < d) 
    stop("nu in dpm.params is required to be a scalar greater than 1.")
  if (dpm.useHyperpriors) {
    ### lambda ~ Gamma(gamma1, gamma2)
    if ((dpm.gamma1 > 0) & (dpm.gamma2 > 0))
      dpm.lambda = -1
    else
      stop("gamma1 and gamma2 in dpm.params are required to be positive scalars.")
    ### Psi ~ Wishart(nu0, Psi0)
    if (dpm.nu0 < d)
      stop("nu0 in dpm.params is required to be a scalar greater than 1.")
  } else {
    if (dpm.lambda > 0) 
      dpm.gamma1 = dpm.gamma2 = -1
    else 
      stop("lambda in dpm.params is required to be a positive scalar.")
    
    dpm.nu0 = -1
  }
  
  
  ##---------------------------------------------------------------------------
  ## define returns from DPMM
  ##---------------------------------------------------------------------------
  if(dpm.pdf) {
    group0.pdfs = matrix(NA, nrow = bart.ndpost, ncol = dpm.ngrid)
    group1.pdfs = matrix(NA, nrow = bart.ndpost, ncol = dpm.ngrid)
  }
  if(dpm.cdf) {
    group0.cdfs = matrix(NA, nrow = (bart.ndpost*dpm.ndpost), ncol = dpm.ngrid)
    group1.cdfs = matrix(NA, nrow = (bart.ndpost*dpm.ndpost), ncol = dpm.ngrid)
    quantiles0 = matrix(NA, nrow = (bart.ndpost*dpm.ndpost), ncol = nprobs)
    quantiles1 = matrix(NA, nrow = (bart.ndpost*dpm.ndpost), ncol = nprobs)
    qtes = matrix(NA, nrow = (bart.ndpost*dpm.ndpost), ncol = nprobs)
  }
  if(dpm.diag) {
    group0.lmpps = matrix(NA, nrow = bart.ndpost, ncol = dpm.ndpost)
    group1.lmpps = matrix(NA, nrow = bart.ndpost, ncol = dpm.ndpost)
  }
  
  ##---------------------------------------------------------------------------
  ## Bayesian Boostrap
  ##---------------------------------------------------------------------------
  diri = .Call("_BNPqte_rdirichlet", bart.ndpost, rep(1.0, dpm.npred))
  
  ##---------------------------------------------------------------------------
  ## print info. for DPMMs
  ##---------------------------------------------------------------------------
  cat("-------------------------------\n")
  cat(" Modeling Outcome ~ Propensity\n")
  cat("-------------------------------\n")
  cat("*****Into main of Weight-Dependent DPMMs\n")
  cat("*****Data: n0(control), n1(treatment): ", n0, ", ", n1, "\n", sep = "")
  if(dpm.cdf | dpm.pdf)
    cat("*****Prediction: type, ngrid, nxpred: ", 
        paste(dpm.type.pred, collapse = ", "), ", ", dpm.ngrid, ", ", dpm.npred, "\n", sep = "")
  else
    cat("*****Prediction: FALSE\n")
  cat("*****Number of clusters:", dpm.nclusters, "\n")
  cat("*****Prior: updateAlpha, useHyperpriors: ", dpm.updateAlpha, ", ", dpm.useHyperpriors, "\n", sep="")
  cat("*****MCMC: nskip, ndpost, keepevery: ", 
      dpm.nskip, ", ", dpm.ndpost, ", ", dpm.keepevery, "\n", sep = "")
  cat("*****Number of DPMMs:", bart.ndpost*2, "\n")
  
  ##---------------------------------------------------------------------------
  ## fitting DPMMs
  ##---------------------------------------------------------------------------
  if(mc.cores == 1) {
    for (i in 1:bart.ndpost) {
      tmp0 = .Call("_BNPqte_cDPMmdensity",
                   n0,
                   d,
                   cbind(y0, propensity0[i, ]),
                   y0,
                   as.matrix(propensity0[i, ]),
                   dpm.diag,
                   dpm.pdf,
                   dpm.cdf,
                   dpm.ngrid,
                   dpm.grid,
                   dpm.npred,
                   as.matrix(dpm.xpred[i, ]),
                   dpm.updateAlpha,
                   dpm.useHyperpriors,
                   dpm.alpha,
                   dpm.a0,
                   dpm.b0,
                   dpm.lambda,
                   dpm.gamma1,
                   dpm.gamma2,
                   dpm.nu0,
                   dpm.nu,
                   dpm.nclusters,
                   dpm.nskip,
                   dpm.ndpost,
                   dpm.keepevery,
                   diri[i, ],
                   probs,
                   nprobs)
      
      tmp1 = .Call("_BNPqte_cDPMmdensity",
                   n1,
                   d,
                   cbind(y1, propensity1[i, ]),
                   y1,
                   as.matrix(propensity1[i, ]),
                   dpm.diag,
                   dpm.pdf,
                   dpm.cdf,
                   dpm.ngrid,
                   dpm.grid,
                   dpm.npred,
                   as.matrix(dpm.xpred[i, ]),
                   dpm.updateAlpha,
                   dpm.useHyperpriors,
                   dpm.alpha,
                   dpm.a0,
                   dpm.b0,
                   dpm.lambda,
                   dpm.gamma1,
                   dpm.gamma2,
                   dpm.nu0,
                   dpm.nu,
                   dpm.nclusters,
                   dpm.nskip,
                   dpm.ndpost,
                   dpm.keepevery,
                   diri[i, ],
                   probs,
                   nprobs)
      
      if(dpm.pdf) {
        group0.pdfs[i, ] = tmp0$predict.pdf.avg
        group1.pdfs[i, ] = tmp1$predict.pdf.avg
      }
      if(dpm.cdf) {
        group0.cdfs[((i-1)*dpm.ndpost+1):(i*dpm.ndpost), ] = tmp0$predict.cdfs
        group1.cdfs[((i-1)*dpm.ndpost+1):(i*dpm.ndpost), ] = tmp1$predict.cdfs
        quantiles0[((i-1)*dpm.ndpost+1):(i*dpm.ndpost), ] = tmp0$predict.quantiles
        quantiles1[((i-1)*dpm.ndpost+1):(i*dpm.ndpost), ] = tmp1$predict.quantiles
        qtes[((i-1)*dpm.ndpost+1):(i*dpm.ndpost), ] = tmp1$predict.quantiles - tmp0$predict.quantiles
      }
      if(dpm.diag) {
        group0.lmpps[i, ] = tmp0$logMPPs
        group1.lmpps[i, ] = tmp1$logMPPs
      }
      
      rm(tmp0)
      rm(tmp1)
      gc()
      cat("-------DPMM fit", 2*i, "out of", 2*bart.ndpost, "\n")
    }
  } else {
    ## parallel computing
    print("Will be added soon.")
  }
  
   
  #---------------------------------------------------------------------------------------------------- 
  # returns
  #---------------------------------------------------------------------------------------------------- 
  res$grid = dpm.grid
  res$xpred = dpm.xpred
  res$type.pred = dpm.type.pred
  
  if(dpm.cdf) {
    res$probs = probs
    res$alphas = alphas
    res$compute.band = compute.band
    if(compute.band)
      res$type.band = type.band
     
    res$control.cdfs = group0.cdfs  # (bart.ndpost*dpm.ndpost) x dpm.ngrid
    res$treatment.cdfs = group1.cdfs
    colnames(res$control.cdfs) = dpm.grid
    colnames(res$treatment.cdfs) = dpm.grid
     
    res$quantiles.avg = cbind(colMeans(quantiles0), colMeans(quantiles1))  # nprobs x 2
    colnames(res$quantiles.avg) = c("control", "treatment")
    rownames(res$quantiles.avg) = paste(probs*100, "%", sep = "")
    
    res$qtes.avg = colMeans(qtes)  # nprobs
    names(res$qtes.avg) = paste(probs*100, "%", sep = "")
     
    if(compute.band) {
      lower.band = matrix(NA, nrow = nprobs, ncol = nalphas)
      upper.band = matrix(NA, nrow = nprobs, ncol = nalphas)
      
      for (i in 1:nprobs) {
        band = .Call(`_BNPqte_credible_interval`, nrow(qtes), qtes[, i], alphas, type.band)
        lower.band[i, ] = band[1, ]
        upper.band[i, ] = band[2, ]
      }
       
      res$lower.band = lower.band
      res$upper.band = upper.band
      rownames(res$lower.band) = paste(probs*100, "%", sep = "")
      rownames(res$upper.band) = paste(probs*100, "%", sep = "")
      colnames(res$lower.band) = alphas
      colnames(res$upper.band) = alphas
    }
  }
  
  if(dpm.pdf) {
    res$pdf.avg = rbind(colMeans(group0.pdfs), colMeans(group1.pdfs))  # 2 x dpm.ngrid
    rownames(res$pdf.avg) = c("control", "treatment")
    colnames(res$pdf.avg) = dpm.grid
  }
   
  if(dpm.diag) {
    res$control.lmpps = group0.lmpps
    res$treatment.lmpps = group1.lmpps
  }
   
  res$n0 = n0
  res$n1 = n1
  res$p = ncol(x)
  
  res$dpm.params = list(nclusters = dpm.nclusters, updateAlpha = dpm.updateAlpha, useHyperpriors = dpm.useHyperpriors,
                        nskip = dpm.nskip, ndpost = dpm.ndpost, keepevery = dpm.keepevery)
  
  attr(res, 'class') <- 'qte'
   
  print("Finished!")
   
  return(res)
  
}