qte = function(y, x, treatment,
               probs=c(0.1, 0.25, 0.50, 0.75, 0.90), 
               compute.band=TRUE, type.band="HPD", alphas=c(0.05),
               bart.link="probit", bart.params=list(split.prob="polynomial"), 
               dpm.params=list(method="truncated", nclusters=50L),
               mc.cores=1L, nice = 19L, seed = 123) {
  
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
  if(is.null(bart.params$split.prob)) {
    bart.split.prob = "polynomial"
  } else {
    bart.split.prob = bart.params$split.prob
  }
  if(!(bart.split.prob %in% c("polynomial", "exponential"))) {
    stop("Only two types of available splitting probability for BART: polynomial or exponential.")
  } else {
    if(bart.split.prob == "polynomial") {
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
    }
    if(bart.split.prob == "exponential") {
      bart.power = -1.0
      if(is.null(bart.params$base)) {
        bart.base = 0.5
      } else {
        bart.base = bart.params$base
      }
    }
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
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  
  #---------------------------------------------------------------------------------------------------- 
  # fit probit or logit BART
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
                     k=bart.k, power=bart.power, base=bart.base, split.prob=bart.split.prob,
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
                     k=bart.k, power=bart.power, base=bart.base, split.prob=bart.split.prob,
                     binaryOffset=bart.binaryOffset,
                     ntree=bart.ntree, numcut=bart.numcut,
                     ndpost=bart.ndpost, nskip=bart.nskip, keepevery=bart.keepevery,
                     nkeeptrain=bart.nkeeptrain, nkeeptest=bart.nkeeptest,
                     nkeeptreedraws=bart.nkeeptreedraws,
                     printevery=bart.printevery, transposed=bart.transposed)
  }
  
  ## estimated propensity scores: bart.ndpost x n
  propensity = bart_obj$prob.train
  res$propensity = propensity
  
  propensity = log(propensity / (1 - propensity))  ## logit transformation on PS
  propensity0 = propensity[, group0, drop = FALSE]
  propensity1 = propensity[, group1, drop = FALSE]
  
  res$bart.parmas = list(link = bart.link, ntree = bart.ntree, 
                         nskip = bart.nskip, ndpost = bart.ndpost, keepevery = bart.keepevery,
                         split.prob = bart.split.prob, sparse = bart.sparse,
                         binaryOffset = bart_obj$binaryOffset, 
                         vip = bart_obj$vip, 
                         within_type_vip = bart_obj$within_type_vip,
                         pvip = bart_obj$pvip,
                         varprob.mean = bart_obj$varprob.mean,
                         mi = bart_obj$mi)
  
  rm(bart_obj)
  gc()
  
  #---------------------------------------------------------------------------------------------------- 
  # fit DPMM
  #---------------------------------------------------------------------------------------------------- 
   
  ##---------------------------------------------------------------------------
  ## process dpm.params: if not specified by users, then use default values
  ##---------------------------------------------------------------------------
  if(is.null(dpm.params$method)) {
    dpm.method = "truncated"
  } else {
    dpm.method = dpm.params$method
  }
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
    dpm.type.pred = c("cdf", "pdf")
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
      dpm.alpha = 1.0
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
      dpm.lambda = rgamma(1, shape = dpm.gamma1, rate = dpm.gamma2)
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
  if(dpm.method == "truncated") {
    cat("*****Posterior sampling method: Blocked Gibbs Sampling with", dpm.nclusters, "clusters\n")
  } else {
    cat("*****Posterior sampling method: Algorithm 8 with m = 1 in Neal (2000)\n")
  }
  cat("*****Prior: updateAlpha, useHyperpriors: ", dpm.updateAlpha, ", ", dpm.useHyperpriors, "\n", sep="")
  cat("*****MCMC: nskip, ndpost, keepevery: ", 
      dpm.nskip, ", ", dpm.ndpost, ", ", dpm.keepevery, "\n", sep = "")
  cat("*****Number of DPMMs:", bart.ndpost*2, "\n")
  
  ##---------------------------------------------------------------------------
  ## fitting DPMMs
  ##---------------------------------------------------------------------------
  if(bart.ndpost == 1)
    mc.cores = 1
  
  if(mc.cores == 1) {
    for (i in 1:bart.ndpost) {
      if(dpm.method == "truncated") {
        tmp0 = .Call("_BNPqte_cDPMmdensity",
                     n0, d, cbind(y0, propensity0[i, ]), y0, as.matrix(propensity0[i, ]),
                     dpm.diag, dpm.pdf, dpm.cdf, dpm.ngrid, dpm.grid, dpm.npred, as.matrix(dpm.xpred[i, ]),
                     dpm.updateAlpha, dpm.useHyperpriors, dpm.alpha, dpm.a0, dpm.b0, dpm.lambda, dpm.gamma1, dpm.gamma2,
                     dpm.nu0, dpm.nu, dpm.nclusters, dpm.nskip, dpm.ndpost, dpm.keepevery, diri[i, ], probs, nprobs)
        
        tmp1 = .Call("_BNPqte_cDPMmdensity",
                     n1, d, cbind(y1, propensity1[i, ]), y1, as.matrix(propensity1[i, ]),
                     dpm.diag, dpm.pdf, dpm.cdf, dpm.ngrid, dpm.grid, dpm.npred, as.matrix(dpm.xpred[i, ]),
                     dpm.updateAlpha, dpm.useHyperpriors, dpm.alpha, dpm.a0, dpm.b0, dpm.lambda, dpm.gamma1, dpm.gamma2,
                     dpm.nu0, dpm.nu, dpm.nclusters, dpm.nskip, dpm.ndpost, dpm.keepevery, diri[i, ], probs, nprobs)
      } else {
        tmp0 = .Call("_BNPqte_cDPMmdensityNeal",
                     n0, d, cbind(y0, propensity0[i, ]), y0, as.matrix(propensity0[i, ]),
                     dpm.diag, dpm.pdf, dpm.cdf, dpm.ngrid, dpm.grid, dpm.npred, as.matrix(dpm.xpred[i, ]),
                     dpm.updateAlpha, dpm.useHyperpriors, dpm.alpha, dpm.a0, dpm.b0, dpm.lambda, dpm.gamma1, dpm.gamma2,
                     dpm.nu0, dpm.nu, dpm.nclusters, dpm.nskip, dpm.ndpost, dpm.keepevery, diri[i, ], probs, nprobs)
        
        tmp1 = .Call("_BNPqte_cDPMmdensityNeal",
                     n1, d, cbind(y1, propensity1[i, ]), y1, as.matrix(propensity1[i, ]),
                     dpm.diag, dpm.pdf, dpm.cdf, dpm.ngrid, dpm.grid, dpm.npred, as.matrix(dpm.xpred[i, ]),
                     dpm.updateAlpha, dpm.useHyperpriors, dpm.alpha, dpm.a0, dpm.b0, dpm.lambda, dpm.gamma1, dpm.gamma2,
                     dpm.nu0, dpm.nu, dpm.nclusters, dpm.nskip, dpm.ndpost, dpm.keepevery, diri[i, ], probs, nprobs)
      }
      
      if(dpm.pdf) {
        if(i == 1) {
          group0.pdfs = tmp0$predict.pdfs
          group1.pdfs = tmp1$predict.pdfs
        } else {
          group0.pdfs = rbind(group0.pdfs, tmp0$predict.pdfs)
          group1.pdfs = rbind(group1.pdfs, tmp1$predict.pdfs)
        }
      }
      if(dpm.cdf) {
        if(i == 1) {
          group0.cdfs = tmp0$predict.cdfs
          group1.cdfs = tmp1$predict.cdfs
          quantiles0 = tmp0$predict.quantiles
          quantiles1 = tmp1$predict.quantiles
          qtes = tmp1$predict.quantiles - tmp0$predict.quantiles
        } else {
          group0.cdfs = rbind(group0.cdfs, tmp0$predict.cdfs)
          group1.cdfs = rbind(group1.cdfs, tmp1$predict.cdfs)
          quantiles0 = rbind(quantiles0, tmp0$predict.quantiles)
          quantiles1 = rbind(quantiles1, tmp1$predict.quantiles)
          qtes = rbind(qtes, tmp1$predict.quantiles - tmp0$predict.quantiles)
        }
      }
      if(dpm.diag) {
        if(i == 1) {
          if(dpm.method == "truncated") {
            group0.lmpps = tmp0$logMPPs
            group1.lmpps = tmp1$logMPPs
          }
          group0.ylogliks = tmp0$ylogliks
          group1.ylogliks = tmp1$ylogliks
        } else {
          if(dpm.method == "truncated") {
            group0.lmpps = rbind(group0.lmpps, tmp0$logMPPs)
            group1.lmpps = rbind(group1.lmpps, tmp1$logMPPs)
          }
          group0.ylogliks = rbind(group0.ylogliks, tmp0$ylogliks)
          group1.ylogliks = rbind(group1.ylogliks, tmp1$ylogliks)
        }
      }
      
      rm(tmp0)
      rm(tmp1)
      gc()
      cat("-------DPMM fit", 2*i, "out of", 2*bart.ndpost, "\n")
    }
  } else {
    if(.Platform$OS.type!='unix')
      stop('parallel::mcparallel/mccollect do not exist on windows')
    
    parallel::mc.reset.stream()
    
    mc.cores.detected = detectCores()
    if(mc.cores > mc.cores.detected) 
      mc.cores = mc.cores.detected
    
    if(mc.cores > bart.ndpost)
      mc.cores = bart.ndpost
    
    mc.bart.ndpost = floor(bart.ndpost / mc.cores)
    
    cat("-------Fit DPMMs in parallel with", mc.cores, "cores...\n")
    
    ### parallel computing
    if(dpm.method == "truncated") {
      for (i in 1:mc.cores) {
        if (i <= bart.ndpost %% mc.cores){
          start_row = (i - 1) * (mc.bart.ndpost + 1) + 1
          end_row = i * (mc.bart.ndpost + 1)
        }
        else{
          start_row = (i - 1) * mc.bart.ndpost + 1 + bart.ndpost %% mc.cores
          end_row = i * mc.bart.ndpost + bart.ndpost %% mc.cores
        }
        
        parallel::mcparallel({psnice(value=nice);
          mc.DPMmdensity(n0, n1, d, y0, y1, propensity0[start_row:end_row, ], propensity1[start_row:end_row, ], 
                         dpm.diag, dpm.pdf, dpm.cdf, dpm.ngrid, dpm.grid, dpm.npred, dpm.xpred[start_row:end_row, ],
                         dpm.updateAlpha, dpm.useHyperpriors, dpm.alpha, dpm.a0, dpm.b0, dpm.lambda, dpm.gamma1, dpm.gamma2, 
                         dpm.nu0, dpm.nu, dpm.nclusters, dpm.nskip, dpm.ndpost, dpm.keepevery, 
                         diri[start_row:end_row, ], probs, nprobs)},
          silent=(i!=1))
      }
    } else {
      for (i in 1:mc.cores) {
        if (i <= bart.ndpost %% mc.cores){
          start_row = (i - 1) * (mc.bart.ndpost + 1) + 1
          end_row = i * (mc.bart.ndpost + 1)
        }
        else{
          start_row = (i - 1) * mc.bart.ndpost + 1 + bart.ndpost %% mc.cores
          end_row = i * mc.bart.ndpost + bart.ndpost %% mc.cores
        }
        
        parallel::mcparallel({psnice(value=nice);
          mc.DPMmdensityNeal(n0, n1, d, y0, y1, propensity0[start_row:end_row, ], propensity1[start_row:end_row, ], 
                             dpm.diag, dpm.pdf, dpm.cdf, dpm.ngrid, dpm.grid, dpm.npred, dpm.xpred[start_row:end_row, ],
                             dpm.updateAlpha, dpm.useHyperpriors, dpm.alpha, dpm.a0, dpm.b0, dpm.lambda, dpm.gamma1, dpm.gamma2, 
                             dpm.nu0, dpm.nu, dpm.nclusters, dpm.nskip, dpm.ndpost, dpm.keepevery, 
                             diri[start_row:end_row, ], probs, nprobs)},
          silent=(i!=1))
      }
    }
    
    ### collect parallel results
    parallel.list <- parallel::mccollect()
    
    if(dpm.pdf) {
      group0.pdfs = parallel.list[[1]]$group0.pdfs
      group1.pdfs = parallel.list[[1]]$group1.pdfs
    }
    
    if(dpm.cdf) {
      group0.cdfs = parallel.list[[1]]$group0.cdfs
      group1.cdfs = parallel.list[[1]]$group1.cdfs
      quantiles0 = parallel.list[[1]]$quantiles0
      quantiles1 = parallel.list[[1]]$quantiles1
      qtes = parallel.list[[1]]$qtes
    }
    
    if(dpm.diag) {
      if(dpm.method == "truncated") {
        group0.lmpps = parallel.list[[1]]$group0.lmpps
        group1.lmpps = parallel.list[[1]]$group1.lmpps
      }
      group0.ylogliks = parallel.list[[1]]$group0.ylogliks
      group1.ylogliks = parallel.list[[1]]$group1.ylogliks
    }
    
    for(i in 2:mc.cores) {
      if(dpm.pdf) {
        group0.pdfs = rbind(group0.pdfs, parallel.list[[i]]$group0.pdfs)
        group1.pdfs = rbind(group1.pdfs, parallel.list[[i]]$group1.pdfs)
      }
      if(dpm.cdf) {
        group0.cdfs = rbind(group0.cdfs, parallel.list[[i]]$group0.cdfs)
        group1.cdfs = rbind(group1.cdfs, parallel.list[[i]]$group1.cdfs)
        quantiles0 = rbind(quantiles0, parallel.list[[i]]$quantiles0)
        quantiles1 = rbind(quantiles1, parallel.list[[i]]$quantiles1)
        qtes = rbind(qtes, parallel.list[[i]]$qtes)
      }
      if(dpm.diag) {
        if(dpm.method == "truncated") {
          group0.lmpps = rbind(group0.lmpps, parallel.list[[i]]$group0.lmpps)
          group1.lmpps = rbind(group1.lmpps, parallel.list[[i]]$group1.lmpps)
        }
        group0.ylogliks = rbind(group0.ylogliks, parallel.list[[i]]$group0.ylogliks)
        group1.ylogliks = rbind(group1.ylogliks, parallel.list[[i]]$group1.ylogliks)
      }
    }
    
    rm(parallel.list)
    gc()
  }
  
  cat("DPMM Finished!\n")
   
  #---------------------------------------------------------------------------------------------------- 
  # returns
  #---------------------------------------------------------------------------------------------------- 
  cat("Collecting returns ...")
  
  res$grid = dpm.grid
  res$xpred = dpm.xpred
  res$type.pred = dpm.type.pred
  res$compute.band = compute.band
  if(compute.band) {
    res$type.band = type.band
    res$alphas = alphas
  }
  
  if(dpm.cdf) {
    res$probs = probs
     
    # cdfs
    res$control.cdfs = group0.cdfs  # (bart.ndpost*dpm.ndpost) x dpm.ngrid
    res$treatment.cdfs = group1.cdfs
    colnames(res$control.cdfs) = dpm.grid
    colnames(res$treatment.cdfs) = dpm.grid
     
    # quantiles
    res$control.quantiles.avg = colMeans(quantiles0)
    res$treatment.quantiles.avg = colMeans(quantiles1)
    names(res$control.quantiles.avg) = paste(probs*100, "%", sep = "")
    names(res$treatment.quantiles.avg) = paste(probs*100, "%", sep = "")
    
    # quantile treatment effects
    res$qtes.avg = colMeans(qtes)  # nprobs
    names(res$qtes.avg) = paste(probs*100, "%", sep = "")
     
    # credible intervals for quantiles and quantile treatment effects
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
        res$control.quantiles.ci = control.quantiles.ci
        
        treatment.quantiles.ci = matrix(NA, nrow = nprobs, ncol = 2)
        for (i in 1:nprobs) {
          band = .Call(`_BNPqte_credible_interval`, nrow(quantiles1), quantiles1[, i], alphas, type.band)
          treatment.quantiles.ci[i, 1] = band[1, ]
          treatment.quantiles.ci[i, 2] = band[2, ]
        }
        rownames(treatment.quantiles.ci) = paste(probs*100, "%", sep = "")
        colnames(treatment.quantiles.ci) = c(alphas, 1-alphas)
        res$treatment.quantiles.ci = treatment.quantiles.ci
        
        qtes.ci = matrix(NA, nrow = nprobs, ncol = 2)
        for (i in 1:nprobs) {
          band = .Call(`_BNPqte_credible_interval`, nrow(qtes), qtes[, i], alphas, type.band)
          qtes.ci[i, 1] = band[1, ]
          qtes.ci[i, 2] = band[2, ]
        }
        rownames(qtes.ci) = paste(probs*100, "%", sep = "")
        colnames(qtes.ci) = c(alphas, 1-alphas)
        res$qtes.ci = qtes.ci
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
        
        res$qtes.ci = qtes.ci
        res$control.quantiles.ci = control.quantiles.ci
        res$treatment.quantiles.ci = treatment.quantiles.ci
      }
    }
  }
  
  if(dpm.pdf) {
    res$control.pdfs.avg = colMeans(group0.pdfs)
    res$treatment.pdfs.avg = colMeans(group1.pdfs)
    names(res$control.pdfs.avg) = dpm.grid
    names(res$treatment.pdfs.avg) = dpm.grid
    
    if(compute.band) {
      if(nalphas == 1) {
        control.pdfs.ci = matrix(NA, nrow = dpm.ngrid, ncol = 2)
        treatment.pdfs.ci = matrix(NA, nrow = dpm.ngrid, ncol = 2)
        
        for (i in 1:dpm.ngrid) {
          band1 = .Call(`_BNPqte_credible_interval`, nrow(group0.pdfs), group0.pdfs[, i], alphas, type.band)
          band2 = .Call(`_BNPqte_credible_interval`, nrow(group1.pdfs), group1.pdfs[, i], alphas, type.band)
          
          control.pdfs.ci[i, 1] = band1[1, ]
          control.pdfs.ci[i, 2] = band1[2, ]
          treatment.pdfs.ci[i, 1] = band2[1, ]
          treatment.pdfs.ci[i, 2] = band2[2, ]
        }
        
        rownames(control.pdfs.ci) = dpm.grid
        colnames(control.pdfs.ci) = c(alphas, 1-alphas)
        res$control.pdfs.ci = control.pdfs.ci
        
        rownames(treatment.pdfs.ci) = dpm.grid
        colnames(treatment.pdfs.ci) = c(alphas, 1-alphas)
        res$treatment.pdfs.ci = treatment.pdfs.ci
      } else {
        control.pdfs.ci = list()
        treatment.pdfs.ci = list()
        
        for (i in 1:nalphas) {
          control.pdfs.ci[[i]] = matrix(NA, nrow = dpm.ngrid, ncol = 2)
          treatment.pdfs.ci[[i]] = matrix(NA, nrow = dpm.ngrid, ncol = 2)
          
          rownames(control.pdfs.ci[[i]]) = dpm.grid
          rownames(treatment.pdfs.ci[[i]]) = dpm.grid
          
          colnames(control.pdfs.ci[[i]]) = c(alphas[i], 1-alphas[i])
          colnames(treatment.pdfs.ci[[i]]) = c(alphas[i], 1-alphas[i])
        }
        
        for (i in 1:dpm.ngrid) {
          band1 = .Call(`_BNPqte_credible_interval`, nrow(group0.pdfs), group0.pdfs[, i], alphas, type.band)
          band2 = .Call(`_BNPqte_credible_interval`, nrow(group1.pdfs), group1.pdfs[, i], alphas, type.band)
          
          for (j in 1:nalphas) {
            control.pdfs.ci[[j]][i, 1] = band1[1, j]
            control.pdfs.ci[[j]][i, 2] = band1[2, j]
            
            treatment.pdfs.ci[[j]][i, 1] = band2[1, j]
            treatment.pdfs.ci[[j]][i, 2] = band2[2, j]
          }
        }
        
        res$control.pdfs.ci = control.pdfs.ci
        res$treatment.pdfs.ci = treatment.pdfs.ci
      }
    }
  }
  
  if(dpm.diag) {
    if(dpm.method == "truncated") {
      res$control.lmpps = group0.lmpps
      res$treatment.lmpps = group1.lmpps
    }
    res$control.ylogliks = group0.ylogliks
    res$treatment.ylogliks = group1.ylogliks
  }
  
  res$n0 = n0
  res$n1 = n1
  res$p = ncol(x)
  
  res$dpm.params = list(method = dpm.method, nclusters = dpm.nclusters, 
                        updateAlpha = dpm.updateAlpha, useHyperpriors = dpm.useHyperpriors,
                        nskip = dpm.nskip, ndpost = dpm.ndpost, keepevery = dpm.keepevery)
  
  attr(res, 'class') <- 'qte'
   
  cat("Finished!\n")
   
  return(res)
  
}