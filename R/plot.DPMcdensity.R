label_both_w_equal <- function(labels) label_both(labels, sep = " = ")

plot.DPMcdensity = function(object, xpred.idx=c(1), 
                            true_pdf=NULL, true_pdf_fun=NULL,
                            true_cdf=NULL, true_cdf_fun=NULL,
                            true_meanReg=NULL, true_meanReg_fun=NULL) {
  
  if(!object$prediction)
    stop("No pdf, cdf or mean regression curve is evaluated in the DPMcdensity() function.")
  
  if(ncol(object$xpred) == 1)
    stringX = paste(object$xpred[xpred.idx, ])
  else
    stringX = paste("c(", apply(as.matrix(object$xpred[xpred.idx, ]), 1, toString), ")", sep="")
  
  
  if("pdf" %in% object$type.pred) {
    # prepare data.frame for plotting
    plot_df = data.frame(grid = rep(object$grid, length(xpred.idx)),
                         x = rep(stringX, each=length(object$grid)),
                         avg = c(t(object$predict.pdf.avg[xpred.idx, ])))
    
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
    if(!is.null(true_pdf))
      plot_df$true = c(t(true_pdf))
    
    if(object$compute.band) {
      plot_df$upper = c(t(object$predict.pdf.upper[xpred.idx, ]))
      plot_df$lower = c(t(object$predict.pdf.lower[xpred.idx, ]))
    }
    
    # plot
    fig = ggplot(plot_df, aes(x=grid, y=avg)) +
      geom_line(aes(y=avg, colour="Posterior Mean\nEstimates")) 
    if(!is.null(true_pdf))
      fig = fig + geom_line(aes(y=true, colour="True Density")) +
      scale_colour_manual("Conditional Density", values=c("black", "red"))
    else
      fig = fig + scale_colour_manual("Conditional Density", values=c("black"))
    if(object$compute.band)
      fig = fig + geom_ribbon(aes(ymin=lower, ymax=upper, fill=paste("95%", object$type.band)), alpha=0.5) +
      scale_fill_manual("Pointwise C.I.",values="grey70")
    
    fig = fig + facet_wrap(~x, labeller = label_both_w_equal) +
      ylab("f(y|x)") +
      xlab("y") +
      theme_bw()
    
    print(fig)
  }
  
  
  if("cdf" %in% object$type.pred) {
    # prepare data.frame for plotting
    plot_df = data.frame(grid = rep(object$grid, length(xpred.idx)),
                         x = rep(stringX, each=length(object$grid)),
                         avg = c(t(object$predict.cdf.avg[xpred.idx, ])))
    
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
    
    if(!is.null(true_cdf))
      plot_df$true = c(t(true_cdf))
    
    if(object$compute.band) {
      plot_df$upper = c(t(object$predict.cdf.upper[xpred.idx, ]))
      plot_df$lower = c(t(object$predict.cdf.lower[xpred.idx, ]))
    }
    
    # plot
    fig = ggplot(plot_df, aes(x=grid, y=avg)) +
      geom_line(aes(y=avg, colour="Posterior Mean\nEstimates")) 
    if(!is.null(true_cdf))
      fig = fig + geom_line(aes(y=true, colour="True CDF")) +
      scale_colour_manual("Conditional CDF", values=c("black", "red"))
    else
      fig = fig + scale_colour_manual("Conditional CDF", values=c("black"))
    if(object$compute.band)
      fig = fig + geom_ribbon(aes(ymin=lower, ymax=upper, fill=paste("95%", object$type.band)), alpha=0.5) +
      scale_fill_manual("Pointwise C.I.",values="grey70")
    
    fig = fig + facet_wrap(~x, labeller = label_both_w_equal) +
      ylab("F(y|x)") +
      xlab("y") +
      theme_bw()
    
    print(fig)
  }
  
  
  if("meanReg" %in% object$type.pred) {
    # prepare data.frame for plotting
    plot_df = data.frame(x = object$xpred,
                         avg = object$predict.meanReg.avg)
    
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
    if(!is.null(true_meanReg))
      plot_df$true = true_meanReg
    
    if(object$compute.band) {
      plot_df$upper = object$predict.meanReg.upper
      plot_df$lower = object$predict.meanReg.lower
    }
    
    plot_df$title = "Conditional Mean Regression Curves"
    
    # plot
    fig = ggplot(plot_df, aes(x=x, y=avg)) +
      geom_line(aes(y=avg, colour="Posterior Mean\nEstimates")) 
    if(!is.null(true_meanReg))
      fig = fig + geom_line(aes(y=true, colour="True Mean")) +
      scale_colour_manual("Conditional Mean", values=c("black", "red"))
    else
      fig = fig + scale_colour_manual("Conditional Mean", values=c("black"))
    if(object$compute.band)
      fig = fig + geom_ribbon(aes(ymin=lower, ymax=upper, fill=paste("95%", object$type.band)), alpha=0.5) +
      scale_fill_manual("Pointwise C.I.",values="grey70")
    
    fig = fig  + facet_grid(. ~ title) +
      ylab("E(Y|x)") +
      xlab("x") +
      theme_bw()
    
    print(fig)
  }
}
