
#' Density for scaled t-distribution used in Fastloc-GPS errors
#' @param y point where we want to evaluate the density
#' @param mean mean parameter
#' @param scale positive scale parameter 
#' @param df degree of freedom
#' @param log logical; if TRUE, probabilities p are given as log(p)
#' @return numeric 
#' @importFrom stats dt
dscaledt <- function(y, mean, scale, df,log=FALSE) {
  
  z <- (y - mean)/scale

  dt(z, df = df, log = log)/scale
}

#' Generate random samples from scaled t-distribution
#' @param n number of samples
#' @param scale positive scale parameter
#' @param df degree of freedom
#' @export
rscaledt <- function(n, scale, df) {
  
   scale * rt(n, df = df)
}

#' Sample from mixture of multivariate student for error in ARGOS positions
#'@param n number of samples
#'@param params parameters a, df, sigma_obs, rho and p in the distribution (see Brost et al. 2015).
#'@return matrix of size n x 2 with samples
#' @importFrom mvtnorm rmvt
#' @importFrom stats runif
#' @export
rmvt_mixture <- function(n,params) {
  
  # Covariance matrices
  Sigma1 <- params$sigma_obs^2 * matrix(c(1, 
    params$rho * sqrt(params$a),
    params$rho * sqrt(params$a), 1), nrow = 2)
  
  # Covariance matrix NW-SE direction
  Sigma2 <- params$sigma_obs^2 * matrix(c(1,
    -params$rho * sqrt(params$a),
    -params$rho * sqrt(params$a),1), nrow = 2)
  
  #Sample from one component of the mixture
  u<-runif(1)
  
  if (u <params$p) {
    return (mvtnorm::rmvt(n,sigma = Sigma1, df = params$df))
  }
  else {
    return (mvtnorm::rmvt(n,sigma = Sigma2, df = params$df))
  }
}
 #' Density of mixture of multivariate student for error in ARGOS positions
 #'@param x point where to evaluate the density
#' @param mean mean of the student distributions
#' @param params parameters a, df, sigma_obs, rho and p in the distribution (see Brost et al. 2015).
#' @param log logical; if TRUE, probabilities p are given as log(p)
#' @return numeric density value
#' @importFrom mvtnorm dmvt
#' 
dmvt_mixture <- function(x, mean, params,log=FALSE) {
  # Covariance matrices
  Sigma1 <- params$sigma_obs^2 * matrix(c(
    1, params$rho * sqrt(params$a),
    params$rho * sqrt(params$a), 1), nrow = 2)
  
  # Covariance matrix NW-SE direction
  Sigma2 <- params$sigma_obs^2 * matrix(c(1,
    -params$rho * sqrt(params$a),
    -params$rho * sqrt(params$a),1), nrow = 2)
  
  # Compute densities
  dens1 <- mvtnorm::dmvt(x - mean, sigma = Sigma1, df = params$df, log = FALSE)
  dens2 <- mvtnorm::dmvt(x - mean, sigma = Sigma2, df = params$df, log = FALSE)
  
  # Return mixture density
  density = params$p * dens1 + (1 - params$p) * dens2
  if (log) { return (log(density))} else {return(density)}
}

#' Plot heatmap of mixture of multivariate Student-t density
#'
#' @param mean vector of length 2: mean of the Student-t distributions
#' @param params list with elements a, df, sigma_obs, rho, p
#' @param xlim,ylim numeric vectors of length 2 defining plotting limits
#' @param n number of grid points per dimension (default = 100)
#' @return ggplot heatmap of the mixture density
#' @import ggplot2
plot_dmvt_mixture <- function(mean, params, xlim, ylim, n = 100) {
  # Create grid
  x_seq <- seq(xlim[1], xlim[2], length.out = n)
  y_seq <- seq(ylim[1], ylim[2], length.out = n)
  grid <- expand.grid(x = x_seq, y = y_seq)
  
  # Evaluate density for each point
  grid$density <- apply(grid, 1, function(row) {
    dmvt_mixture(as.numeric(row), mean = mean, params = params)
  })
  
  ggplot(grid, aes(x = x, y = y, fill = density)) +
    geom_raster(interpolate = FALSE) +
    scale_fill_viridis_c(option = "viridis") +
    coord_fixed() +
    labs(
      x = expression(x[1]),
      y = expression(x[2]),
      fill = "Density"
    ) +
    theme_minimal(base_size = 12)
}

