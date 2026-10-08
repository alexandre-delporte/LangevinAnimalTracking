
#' Get estimate of the trajectory (conditional mean) from particle filter 
#' @param PF_results List as returned by \code{particle_filter_2D}
#' @return a matrix of estimated filtered states based on conditional mean
#' @export
get_PF_mean<-function(PF_results) {

  particles <- PF_results$particles
  weights <- PF_results$weights
  
  N<-dim(particles)[3]
  num_states<-dim(particles)[2]
  
  estimated_PF_trajectory <- matrix(0, nrow = N, ncol = num_states) 
  for (t in 1:nrow(estimated_PF_trajectory)) {
    estimated_PF_trajectory[t, ] <- colSums(particles[, , t] * weights[, t])
  }
  
  colnames(estimated_PF_trajectory)<-c("X1","X2","V1","V2")
  return (estimated_PF_trajectory)
}

#' Plot trajectory within polygon boundary
#' Visualizes and compares true and estimated trajectories within a specified domain
#'
#' @param U_true Matrix of true states (x1, x2) from simulation
#' @param U_est Matrix of estimated states (x1, x2) from particle filter
#' @param polygon SpatialPolygons object defining the boundary (from `sp` package)
#' @param true_line_opacity Opacity for true trajectory line (0-1, default: 0.5)
#' @param est_line_opacity Opacity for estimated trajectory line (0-1, default: 0.1)
#' @param zoom_box Optional list with elements `xmin`, `xmax`, `ymin`, `ymax` 
#' to zoom in on a specific area of the plot (default: NULL for no zoom)
#' @param labels Character vector of length 2 for 
#' legend labels of true and estimated paths
#'
#' @return A ggplot object showing:
#' \itemize{
#'   \item Polygon boundary 
#'   \item True trajectory 
#'   \item Estimated trajectory 
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom sp SpatialPolygons
compare_trajectories_in_polygon <- function(U_true, U_est,
                                            polygon, true_line_opacity = 0.5,
                                            est_line_opacity = 0.1,
                                            zoom_box = NULL,
                                            labels=c("True path","Estimated path")) {

  # Extract the coordinates of the polygon
  polygon_coords <- as.data.frame(polygon@coords)

  # Create data frames for the true and estimated trajectory points
  U_true_df <- as.data.frame(U_true)
  colnames(U_true_df) <- c("X1", "X2")
  U_true_df$Type <- labels[1]

  U_est_df <- as.data.frame(U_est)
  colnames(U_est_df) <- c("X1", "X2")
  U_est_df$Type <- labels[2]

  # Combine data for easier legend handling
  combined_df <- rbind(U_true_df, U_est_df)
  combined_df$Type <- factor(combined_df$Type, levels = labels)

  colors<-c("black","orange")
  names(colors)=labels

  p <- ggplot() +
    geom_polygon(data = polygon_coords, aes(x = V1, y = V2), fill = NA, color = "grey40") +
    geom_point(data = combined_df, aes(x = X1, y = X2, color = Type), size = 1, alpha = 0.2) +
    geom_path(data = U_true_df, aes(x = X1, y = X2, color = Type), size = 1, alpha = true_line_opacity) +
    geom_path(data = U_est_df, aes(x = X1, y = X2, color = Type), size = 1, alpha = est_line_opacity) +
    scale_color_manual(values = colors,
                       name = "Trajectory") +
    labs(title ="",
         x = "X1",
         y = "X2") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Apply zoom if zoom_box is provided
  if (!is.null(zoom_box)) {
    p <- p + 
      coord_cartesian(xlim = c(zoom_box$xmin, zoom_box$xmax),
                      ylim = c(zoom_box$ymin, zoom_box$ymax))
  }
  
  print(p)
  return(p)
}

#' Plot cloud of particles to see where particles are landing 
#' at a specific time step
#' @param t0 time 
#' @param data data frame of observations with columns Y1 and Y2
#' @param particles array of particles as returned by \code{particle_filter2D_cpp}
#' @param weights matrix of particle weights as returned by \code{particle_filter2D_cpp}
#' @param polygon SpatialPolygons object defining the boundary (from `sp` package`)
#' @export
#' @import ggplot2
#' @importFrom sp SpatialPolygons
#'
particle_cloud<-function(t0,data,particles,weights,polygon=NULL) {

  dfp <- data.frame(x = particles[,1,t0], y = particles[,2,t0], w=weights[,t0])
  p<-ggplot() + geom_point(data=dfp, aes(x,y,size=w),alpha=0.6) +
    geom_point(aes(x=data$Y1[t0], y=data$Y2[t0]), colour="red", size=3) +
    coord_equal() + ggtitle(paste("Particle cloud at t=", t0))
  if (!(is.null(polygon))) {
    polygon_coords <- as.data.frame(polygon@coords)
    p<-p+geom_polygon(data = polygon_coords, aes(x = V1, y = V2), fill = "lightblue", color = "blue",alpha=0.1)
  }
  p
}
