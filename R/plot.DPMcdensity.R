plot.DPMcdensity = function(object, xpred.idx=c(1), 
                            true_pdf=NULL, true_pdf_fun=NULL,
                            true_cdf=NULL, true_cdf_fun=NULL,
                            true_meanReg=NULL, true_meanReg_fun=NULL) {
  if(!object$prediction)
    stop("No pdf, cdf or mean regression curve is evaluated in the DPMcdensity() function.")
  
  
  if("pdf" %in% object$type.pred) {
    
    if(!is.null(true_pdf)) {
      if((nrow(true_pdf)!=length(xpred.idx)) | (ncol(true_pdf)!=length(object$grid)))
        stop("true_pdf is required to be a matrix with nrow=length(xpred.idx) and ncol=length(object$grid).")
    } else {
      if(!is.null(true_pdf_fun)) {
        true_pdf = matrix(NA, nrow = length(xpred.idx), ncol = length(object$grid))
        for (i in 1:nrow(true_pdf)) {
          for (j in 1:ncol(true_pdf)) {
            true_pdf[i, j] = true_pdf_fun(object$grid[j], object$xpred[xpred.idx[i], ])
          }
        }
      }
    }
    
    for (i in 1:length(xpred.idx)) {
      if(object$compute.band) {
        plot(object$grid, object$predict.pdf.upper[xpred.idx[i], ], type="l", lty=2,
             ylab = "f(y|x)", 
             xlab = ifelse(ncol(object$xpred)==1,
                           paste0("y|x=", round(object$xpred[xpred.idx[i], ], 3)),
                           paste0("y|x=c(", paste(round(object$xpred[xpred.idx[i], ], 3), collapse = ","), ")")),
             ylim = c(0, max(object$predict.pdf.upper[xpred.idx[i], ])+0.5))
        lines(object$grid, object$predict.pdf.lower[xpred.idx[i], ], lty=2)
        lines(object$grid, object$predict.pdf.avg[xpred.idx[i], ])
      } else {
        plot(object$grid, object$predict.pdf.avg[xpred.idx[i], ], type="l",
             ylim = c(0, max(object$predict.pdf.avg[xpred.idx[i], ])+1.0),
             ylab = "f(y|x)", 
             xlab = ifelse(ncol(object$xpred)==1,
                           paste0("y|x=", round(object$xpred[xpred.idx[i], ], 3)),
                           paste0("y|x=c(", paste(round(object$xpred[xpred.idx[i], ], 3), collapse = ","), ")")))
      }
      if(!is.null(true_pdf))
        lines(object$grid, true_pdf[i, ], col="red")
    }
  }
  
  
  if("cdf" %in% object$type.pred) {
    
    if(!is.null(true_cdf)) {
      if((nrow(true_cdf)!=length(xpred.idx)) | (ncol(true_cdf)!=length(object$grid)))
        stop("true_cdf is required to be a matrix with nrow=length(xpred.idx) and ncol=length(object$grid).")
    } else {
      if(!is.null(true_cdf_fun)) {
        true_cdf = matrix(NA, nrow = length(xpred.idx), ncol = length(object$grid))
        for (i in 1:nrow(true_cdf)) {
          for (j in 1:ncol(true_cdf)) {
            true_cdf[i, j] = true_cdf_fun(object$grid[j], object$xpred[xpred.idx[i], ])
          }
        }
      }
    }
    
    for (i in 1:length(xpred.idx)) {
      if(object$compute.band) {
        plot(object$grid, object$predict.cdf.upper[xpred.idx[i], ], type="l", lty=2,
             ylab = "F(y|x)", 
             xlab = ifelse(ncol(object$xpred)==1,
                           paste0("y|x=", round(object$xpred[xpred.idx[i], ], 3)),
                           paste0("y|x=c(", paste(round(object$xpred[xpred.idx[i], ], 3), collapse = ","), ")")),
             ylim = c(0, 1.0))
        lines(object$grid, object$predict.cdf.lower[xpred.idx[i], ], lty=2)
        lines(object$grid, object$predict.cdf.avg[xpred.idx[i], ])
      } else {
        plot(object$grid, object$predict.cdf.avg[xpred.idx[i], ], type="l",
             ylim = c(0, 1.0),
             ylab = "F(y|x)", 
             xlab = ifelse(ncol(object$xpred)==1,
                           paste0("y|x=", round(object$xpred[xpred.idx[i], ], 3)),
                           paste0("y|x=c(", paste(round(object$xpred[xpred.idx[i], ], 3), collapse = ","), ")")))
      }
      if(!is.null(true_cdf))
        lines(object$grid, true_cdf[i, ], col="red")
    }
  }
  
  
  if("meanReg" %in% object$type.pred) {
    
    if(!is.null(true_meanReg)) {
      if(length(true_meanReg) != nrow(object$xpred))
        stop("true_meanReg is required to be a vector with length=nrow(object$xpred).")
    } else {
      if(!is.null(true_meanReg_fun)) {
        true_meanReg = c()
        for (i in 1:nrow(object$xpred)) {
          true_meanReg[i] = true_meanReg_fun(object$xpred[i, ])
        }
      }
    }
    
    if(ncol(object$xpred) == 1) {
      if(object$compute.band) {
        plot(object$xpred, object$predict.meanReg.upper, type="l", lty=2,
             xlab="x", ylab="E(y|x)",
             ylim = c(min(object$predict.meanReg.lower)-0.5, max(object$predict.meanReg.upper)+0.5))
        lines(object$xpred, object$predict.meanReg.lower, lty=2)
        lines(object$xpred, object$predict.meanReg.avg)
      } else {
        plot(object$xpred, object$predict.meanReg.avg, type="l",
             xlab="x", ylab="E(y|x)",
             ylim = c(min(object$predict.meanReg.avg)-0.5, max(object$predict.meanReg.avg)+0.5))
      }
      if(!is.null(true_meanReg))
        lines(object$xpred, true_meanReg, col="red")
    } else {
      if(object$compute.band) {
        plot(object$predict.meanReg.upper, type="l", lty=2,
             xlab="x-index", ylab="E(y|x)",
             ylim = c(min(object$predict.meanReg.lower)-0.5, max(object$predict.meanReg.upper)+0.5))
        lines(object$predict.meanReg.lower, lty=2)
        lines(object$predict.meanReg.avg)
      } else {
        plot(object$predict.meanReg.avg, type="l", 
             xlab="x-index", ylab="E(y|x)",
             ylim = c(min(object$predict.meanReg.avg)-0.5, max(object$predict.meanReg.avg)+0.5))
      }
      if(!is.null(true_meanReg))
        lines(true_meanReg, col="red")
    }
  }
  
  
}
