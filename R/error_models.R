
#' Brost et al. (2015) ARGOS location class error parameters
#'
#' Returns a named list of error distribution parameters for each ARGOS
#' location class, based on the mixture of bivariate Student-t model of
#' Brost et al. (2015). The parameters are \code{sigma_obs}, \code{rho},
#' \code{a}, \code{df}, and \code{p} as used in \code{dmvt_mixture}.
#'
#' Location classes in increasing precision order: "B", "A", "0", "1", "2", "3".
#' Class "Z" (invalid/deleted) is treated as class "B" (least precise).
#'
#' @return A named list with one element per location class. Each element is a
#'   list with components \code{sigma_obs}, \code{rho}, \code{a}, \code{df},
#'   \code{p}.
#' @references Appendix D in Brost, B.M., Hooten, M.B., Hanks, E.M., & Small, R.J. (2015).
#'   Animal movement constraints improve resource selection inference in the
#'   presence of telemetry error. \emph{Ecology}, 96(10), 2590-2597.
#' @export
argos_class_params <- function() {
  list(
    "3" = list(sigma_obs = 2.259, rho = 0.85, a = 0.69, df = 16.5, p = 0.5),
    "2" = list(sigma_obs = 1.089, rho = 0.73, a = 0.42, df = 2.3, p = 0.5),
    "1" = list(sigma_obs = 1.365, rho = 0.40, a = 0.64, df = 3, p = 0.5),
    "0" = list(sigma_obs = 2.720, rho = 0.16, a = 0.50, df = 3, p = 0.5),
    "A" = list(sigma_obs = 2.702, rho = 0.21, a = 0.90, df = 3, p = 0.5),
    "B" = list(sigma_obs = 13.338, rho = 0.30, a = 0.74, df = 3, p = 0.5)
  )
}

#' Build a per-observation error parameter list from ARGOS location classes
#'
#' Constructs the \code{obs_error_params} argument for \code{particle_filter2D}
#' from a vector of ARGOS location class codes. Each element of the returned
#' list contains the error distribution parameters for the corresponding
#' observation.
#'
#' @param loc_classes Character vector of length \code{N} giving the ARGOS
#'   location class for each observation. Valid values are \code{"3"},
#'   \code{"2"}, \code{"1"}, \code{"0"}, \code{"A"}, \code{"B"}, \code{"Z"}.
#' @param class_params Named list of error parameters per class, as returned
#'   by \code{argos_class_params()}. Override this to use custom parameters.
#'
#' @return A list of length \code{N}, where each element is an error parameter
#'   list (with \code{sigma_obs}, \code{rho}, \code{a}, \code{df}, \code{p}).
#'
#' @examples
#' \dontrun{
#' classes <- c("2", "1", "3", "0", "B", "A", "2")
#' obs_ep  <- make_argos_obs_params(classes)
#' # Pass to particle filter:
#' particle_filter2D(..., error_dist = "argos",
#'                   error_params = argos_class_params()[["B"]],  # fallback
#'                   obs_error_params = obs_ep)
#' }
#' @export
make_argos_obs_params <- function(loc_classes,
                                  class_params = argos_class_params()) {
  unknown <- setdiff(unique(loc_classes), names(class_params))
  if (length(unknown) > 0) {
    warning("Unknown ARGOS location class(es): ",
            paste(unknown, collapse = ", "),
            ". Falling back to class 'B' parameters.")
    loc_classes[loc_classes %in% unknown] <- "B"
  }
  lapply(loc_classes, function(cls) class_params[[cls]])
}

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

