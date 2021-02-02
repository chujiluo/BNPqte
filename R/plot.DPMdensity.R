plot.DPMdensity = function(object, diff=FALSE, true_density=NULL, true_density_fun=NULL) {
  if(!object$prediction)
    stop("No density is evaluated in the DPMdensity() function.")
  
  
  if(diff) {
    if(is.null(true_density) & is.null(true_density_fun)) {
      
      stop("Either true_density or true_density_fun is required to provide when diff=TRUE.")
      
    } else {
      
      if(!is.null(true_density)) {
        if((nrow(true_density)!=length(object$grid1)) | (ncol(true_density)!=length(object$grid2)))
          stop("true_density is required to be a matrix with nrow=length(object$grid1) and ncol=length(object$grid2).")
      } else {
        true_density = matrix(NA, nrow = length(object$grid1), ncol = length(object$grid2))
        for (i in 1:nrow(true_density)) {
          for (j in 1:ncol(true_density)) {
            true_density[i, j] = true_density_fun(object$grid1[i], object$grid2[j])
          }
        }
      }
      
      # plot the difference density (true-predicted) surface
      fig = plot_ly(x = object$grid1, 
                    y = object$grid2, 
                    z = (true_density - object$predict.densities.mean)) %>% 
        add_surface() %>% 
        layout(title = "(True Density - Predicted Density) Surface", scene = list(xaxis = list(title = "Y1"),
                                                                                  yaxis = list(title = "Y2"),
                                                                                  zaxis = list(title = "Density")))
    }
    
  } else {
    # plot predicted density surface only
    fig = plot_ly(x = object$grid1, 
                  y = object$grid2, 
                  z = object$predict.densities.mean) %>% 
      add_surface() %>% 
      layout(title = "Predicted Density Surface", scene = list(xaxis = list(title = "Y1"),
                                                               yaxis = list(title = "Y2"),
                                                               zaxis = list(title = "Density")))
  }
  
  fig
  
}
