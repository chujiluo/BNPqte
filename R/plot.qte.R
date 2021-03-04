plot.qte = function(object, 
                    true.fy1=NULL, true.fy0=NULL,
                    true.quantile1=NULL, true.quantile0=NULL, true.qte=NULL) {
  
  if(object$compute.band) {
    nalphas = length(object$alphas)
    if(nalphas > 1)
      warning("Only the credible interval corresponding to the first alpha in object$alphas is plotted.")
  }
  
  # pdf plot
  if("pdf" %in% object$type.pred) {
    ngrid = length(object$grid)
    
    ## prepare data frame
    plot_df = data.frame(grid = rep(object$grid, 2),
                         avg = c(object$control.pdfs.avg, object$treatment.pdfs.avg),
                         grp = c(rep("T=0", ngrid), rep("T=1", ngrid)))
    
    if((!is.null(true.fy1)) & (!is.null(true.fy0))) {
      true.fy1 = as.vector(true.fy1)
      true.fy0 = as.vector(true.fy0)
      
      if(length(true.fy1) != ngrid)
        stop("true.fy1 is required to be a vector of marginal densities evaluated at object$grid, 
             should have the same length as object$grid.")
      if(length(true.fy0) != ngrid)
        stop("true.fy0 is required to be a vector of marginal densities evaluated at object$grid, 
             should have the same length as object$grid.")
      
      plot_df$true = c(true.fy0, true.fy1)
    }
    
    if(object$compute.band) {
      if(nalphas == 1) {
        plot_df$lower = c(object$control.pdfs.ci[, 1], object$treatment.pdfs.ci[, 1])
        plot_df$upper = c(object$control.pdfs.ci[, 2], object$treatment.pdfs.ci[, 2])
      } else {
        plot_df$lower = c(object$control.pdfs.ci[[1]][, 1], object$treatment.pdfs.ci[[1]][, 1])
        plot_df$upper = c(object$control.pdfs.ci[[1]][, 2], object$treatment.pdfs.ci[[1]][, 2])
      }
    }
    
    ## plot
    fig = ggplot(plot_df, aes(x=grid, y=avg)) +
      geom_line(aes(y=avg, colour="Estimated")) 
    if((!is.null(true.fy1)) & (!is.null(true.fy0)))
      fig = fig + geom_line(aes(y=true, colour="True")) +
      scale_colour_manual("Density", values=c("black", "red"))
    else
      fig = fig + scale_colour_manual("Density", values=c("black"))
    if(object$compute.band)
      fig = fig + geom_ribbon(aes(ymin=lower, ymax=upper, 
                                  fill=paste(paste0((1-object$alphas[1])*100, "%"), 
                                             object$type.band)), alpha=0.5) +
      scale_fill_manual("Pointwise C.I.",values="grey70")
    
    fig = fig + facet_wrap(~grp, labeller = label_value) +
      ylab(expression(f[t](y))) +
      xlab("y(t)") +
      theme_bw()
    
    print(fig)
  }
  
  # quantile plot
  if(("cdf" %in% object$type.pred) & ((!is.null(true.quantile1)) & (!is.null(true.quantile0)))) {
    nprobs = length(object$probs)
    
    true.quantile1 = as.vector(true.quantile1)
    true.quantile0 = as.vector(true.quantile0)
    
    if(length(true.quantile1) != nprobs)
      stop("true.quantile1 is required to be a vector of quantiles corresponding to object$probs, 
             should have the same length as object$probs")
    if(length(true.quantile0) != nprobs)
      stop("true.quantile0 is required to be a vector of quantiles corresponding to object$probs, 
             should have the same length as object$probs")
    
    ## prepare data frame
    plot_df = data.frame(avg = c(object$control.quantiles.avg, object$treatment.quantiles.avg),
                         true = c(true.quantile0, true.quantile1),
                         grp = c(rep("Y(0)", nprobs), rep("Y(1)", nprobs)))
    
    if(object$compute.band) {
      if(nalphas == 1) {
        plot_df$lower = c(object$control.quantiles.ci[, 1], object$treatment.quantiles.ci[, 1])
        plot_df$upper = c(object$control.quantiles.ci[, 2], object$treatment.quantiles.ci[, 2])
      } else {
        plot_df$lower = c(object$control.quantiles.ci[[1]][, 1], object$treatment.quantiles.ci[[1]][, 1])
        plot_df$upper = c(object$control.quantiles.ci[[1]][, 2], object$treatment.quantiles.ci[[1]][, 2])
      }
    }
    
    
    ## plot
    fig = ggplot(plot_df, aes(x=true, y=avg)) +
      geom_line(aes(y=avg, colour = "Q-Q")) + 
      geom_abline(aes(slope = 1, intercept = 0, colour="Y=X"), show.legend = FALSE) + 
      scale_colour_manual("", values=c("black", "red"))
    if(object$compute.band)
      fig = fig + geom_ribbon(aes(ymin=lower, ymax=upper, 
                                  fill=paste(paste0((1-object$alphas[1])*100, "%"), 
                                             object$type.band)), alpha=0.5) +
      scale_fill_manual("",values="grey70")
    fig = fig + facet_wrap(~grp, labeller = label_value) +
      ylab("Estimated Quantiles") +
      xlab("True Quantiles") +
      theme_bw()
    
    print(fig)
  }
  
  # width of qte C.I.s
  if(("cdf" %in% object$type.pred) & object$compute.band) {
    ## prepare data frame
    width = c()
    if(nalphas == 1) {
      width = object$qtes.ci[, 2] - object$qtes.ci[, 1]
    } else {
      for (i in 1:nalphas) {
        width = c(width, (object$qtes.ci[[i]][, 2] - object$qtes.ci[[i]][, 1]))
      }
    }
    
    plot_df = data.frame(probs = rep(paste0(object$probs*100, "th"), nalphas),
                         alphas = as.factor(rep(object$alphas, each = length(object$probs))),
                         width = width)
    plot_df$title = "Widths of C.I. for QTEs"
    
    ## plot
    fig = ggplot(plot_df, aes(x=probs, y=width, group = alphas)) +
      geom_point(aes(shape = alphas, colour = alphas)) 
    
    fig = fig + facet_grid(. ~ title) +
      ylab("Width of CI") +
      xlab("Probability") +
      theme_bw()
    
    print(fig)
  }
  
  # bias of qte
  if(("cdf" %in% object$type.pred) & (!is.null(true.qte))) {
    true.qte = as.vector(true.qte)
    if(length(true.qte) != length(object$probs))
      stop("true.qte is required to be a vector of quantile treatment effects corrsponding to object$probs, 
             should have the same length as object$probs")
    
    plot_df = data.frame(probs = paste0(object$probs*100, "th"),
                         bias = (true.qte - object$qtes.avg))
    plot_df$title = "Bias of QTEs"
    
    ## plot
    fig = ggplot(plot_df, aes(x=probs, y=bias)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype='dotted', col = 'red', show.legend = FALSE)
    
    fig = fig + facet_grid(. ~ title) +
      ylab("Bias") +
      xlab("Probability") +
      theme_bw()
    
    print(fig)
  }
  
}
