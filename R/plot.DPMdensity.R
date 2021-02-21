plot.DPMdensity = function(object, diff=FALSE, true_density=NULL, true_density_fun=NULL, type.plot="contour") {
  if(!object$prediction)
    stop("No density is evaluated in the DPMdensity() function.")
  
  if(!(type.plot %in% c("surface", "contour")))
    stop("Only two types of plots are available for 2D density: surface or contour.")
  
  
  if(diff) {
    if(is.null(true_density) & is.null(true_density_fun)) {
      stop("Either true_density or true_density_fun is required to provide when diff=TRUE.")
    } else {
      # evaluate the true densities
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
      
      if(type.plot == "surface") {
        # plot the difference density (true-predicted) surface (shown in Viewer pane)
        # note: in plot_ly, x corresponds to columns of z and y corresponds to rows of z
        fig = plot_ly(x = object$grid2, 
                      y = object$grid1, 
                      z = (true_density - object$predict.pdf.avg)) %>% 
          add_surface() %>% 
          layout(title = "Difference btw. True and Posterior Mean Estimates",
                 scene = list(xaxis = list(title = "Y2"),
                              yaxis = list(title = "Y1"),
                              zaxis = list(title = "Density")))
      } else {
        # plot true density, predicted density and density difference contours
        grid_df = expand.grid(grid1 = object$grid1, grid2 = object$grid2)
        plot_df = data.frame(Y1 = rep(grid_df$grid1, 3),
                             Y2 = rep(grid_df$grid2, 3),
                             den = c(as.vector(object$predict.pdf.avg), 
                                     as.vector(true_density),
                                     as.vector(true_density-object$predict.pdf.avg)),
                             type = c(rep("Posterior Mean Estimates", length(as.vector(object$predict.pdf.avg))),
                                      rep("True Density", length(as.vector(object$predict.pdf.avg))),
                                      rep("Density Difference", length(as.vector(object$predict.pdf.avg)))))
        
        fig = ggplot(plot_df, aes(Y1, Y2, z = den))  + 
          geom_contour(aes(colour = ..level..), show.legend = TRUE) + 
          facet_grid(.~type) + theme_bw()
        fig$labels$colour = "Density"
      }
      
    }
    
  } else {
    
    if(type.plot == "surface") {
      # plot predicted density surface only
      fig = plot_ly(x = object$grid2, 
                    y = object$grid1, 
                    z = object$predict.pdf.avg) %>% 
        add_surface() %>% 
        layout(title = "Posterior Mean Estimates",
               scene = list(xaxis = list(title = "Y2"),
                            yaxis = list(title = "Y1"),
                            zaxis = list(title = "Density")))
    } else {
      # plot predicted density contours
      plot_df = expand.grid(Y1 = object$grid1, Y2 = object$grid2)
      plot_df$den = as.vector(object$predict.pdf.avg)
      plot_df$title = "Posterior Mean Estimates"
      
      fig = ggplot(plot_df, aes(Y1, Y2, z = den))  + 
        geom_contour(aes(colour = ..level..), show.legend = TRUE) + 
        theme_bw() +
        facet_grid(. ~ title)
      fig$labels$colour = "Density"
    }
  }
  
  fig
  
}
