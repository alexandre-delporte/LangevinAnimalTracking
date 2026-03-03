
#' Stochastic Gradient Descent with Fisher Information Preconditioning for SDE
#' parameter estimation
#'
#' Implements the Fisher-SGD algorithm of Baey et al. (2023) for maximum
#' likelihood estimation of parameters in latent variable models. The algorithm
#' preconditions the stochastic gradient by a positive definite estimate of the
#' Fisher information matrix (FIM), which scales the parameter updates across
#' directions and ensures descent. The learning rate follows a three-phase
#' schedule: exponential growth (pre-heating), a constant phase at 1 (heating),
#' and a Robbins-Monro decreasing phase. 
#'
#'
#' @param data A data frame of observations with columns \code{time},
#'   \code{Y1}, and \code{Y2} (observed positions).
#' @param sde_params A named list of initial SDE parameter values. Must
#'   contain at least \code{tau} (persistence time scale), \code{nu} (speed
#'   scale), and \code{omega} (angular velocity). Names are used to label
#'   columns in the output matrix.
#' @param fixpar A character vector of parameter names to hold fixed during
#'   optimisation (their gradients and updates are zeroed out). Default
#'   \code{NULL} means all parameters are estimated.
#' @param SGD_iter A positive integer giving the total number of SGD
#'   iterations to run (across all three phases).
#' @param potential_params A named list of parameters for the mixture of
#'   Gaussian potential \eqn{H}, with elements:
#'   \describe{
#'     \item{\code{alpha}}{Numeric vector of mixture weights.}
#'     \item{\code{B}}{List of precision matrices for each Gaussian component.}
#'     \item{\code{x_star}}{Matrix of attraction centres (one row per
#'       component).}
#'   }
#' @param error_dist A character string specifying the observation error
#'   distribution (e.g. \code{"gaussian"}, \code{"argos"}).
#' @param error_params A named list of parameters for the observation error
#'   distribution (e.g. \code{list(sigma = 0.2)}).
#' @param scheme A character string specifying the numerical splitting scheme
#'   used to approximate the SDE transition density. Either
#'   \code{"Lie-Trotter"} (default) or \code{"Strang"}.
#' @param polygon A \code{SpatialPolygons} object defining the spatial domain
#'   \eqn{\mathcal{D}} within which animal movement is constrained.
#' @param U0 A numeric vector of length 4 giving the initial state
#'   \eqn{(x_1, x_2, v_1, v_2)} (position and velocity) at time 0.
#' @param lambda A positive scalar giving the domain penalisation parameter
#'   \eqn{\lambda} used in the reflected Langevin push.
#' @param num_particles A positive integer giving the number of particles
#'   used in the particle filter.
#' @param split_around_fixed_point Logical. If \code{TRUE}, the potential is
#'   split around its fixed point. Currently not implemented and will throw an
#'   error. Default \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, progress messages are printed at
#'   each iteration including the current learning rate, gradient, Fisher
#'   matrix eigenvalues, step direction, and updated parameters. Default
#'   \code{FALSE}.
#' @param gamma0 A small positive scalar giving the initial learning rate at
#'   the start of the pre-heating phase. The authors recommend \eqn{10^{-4}}.
#'   Default \code{1e-4}.
#' @param K_preheat A positive integer giving the number of pre-heating
#'   iterations. The authors recommend 1000; this can be increased (with a
#'   corresponding decrease in \code{gamma0}) if early iterations are
#'   unstable. Default \code{1000}.
#' @param alpha A scalar in \eqn{(1/2, 1]} controlling the decay rate of the
#'   learning rate in the decreasing phase: \eqn{\gamma_k = (k -
#'   K_{\text{end-heating}})^{-\alpha}}. Must satisfy \eqn{\alpha > 1/2} to
#'   ensure \eqn{\sum \gamma_k^2 < \infty}. Default \code{2/3}.
#' @param C_heating A small positive scalar giving the smoothing constant of
#'   the 3rd-order exponential mean filter applied to the gradient norm during
#'   the heating phase. The heating phase ends when the filtered norm stops
#'   decreasing. The authors recommend \eqn{1/1000}. Default \code{1/1000}.
#' @param smoothing Logical. If \code{TRUE}, the latent trajectory is obtained
#'   via forward-filtering backward-sampling (FFBS) rather than from the
#'   particle filter mean, giving lower-variance gradient estimates at higher
#'   computational cost. Default \code{FALSE}.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{theta}}{A matrix of dimension \code{(SGD_iter + 1)} \eqn{\times}
#'     \code{p} containing the parameter estimates at each iteration. Row 1
#'     corresponds to the initial value \code{sde_params}. Column names match
#'     \code{names(sde_params)}.}
#'   \item{\code{Delta}}{A matrix of dimension \code{(n - 1)} \eqn{\times}
#'     \code{p} containing the running stochastic approximation of the
#'     per-step score contributions at the final iteration. Can be used to
#'     reconstruct the final Fisher information matrix estimate as
#'     \eqn{I = \Delta^\top \Delta}.}
#'   \item{\code{K_end_heating}}{An integer giving the iteration at which the
#'     heating phase ended and the decreasing learning rate schedule began.
#'     \code{NA} if the heating phase had not yet ended when the algorithm
#'     terminated.}
#' }
#'
#' @references
#' Baey, C., Delattre, M., Kuhn, E., Leger, J.-B., and Lemler, S. (2023).
#' Efficient preconditioned stochastic gradient descent for estimation in
#' latent variable models. \emph{Proceedings of the 40th International
#' Conference on Machine Learning}, PMLR 202.
#'
#' @seealso
#' \code{\link{particle_filter2D_cpp}} for the underlying particle filter,
#' \code{\link{forward_filtering_backward_sampling}} for the smoother,
#' \code{\link{llk_gradient}} for the pseudo-likelihood gradient computation.
#'
#' @export
SGD_Fisher <- function(data, sde_params, fixpar = NULL, SGD_iter,
                       potential_params = NULL, error_dist = NULL, error_params = NULL,
                       scheme = "Lie-Trotter", polygon, U0, lambda,
                       num_particles, split_around_fixed_point = FALSE,
                       verbose = FALSE, gamma0 = 1e-4, K_preheat = 1000,
                       alpha = 2/3, C_heating = 1/1000, smoothing = FALSE) {
  
  if (split_around_fixed_point) {
    stop("Not implemented yet with splitting around fixed point")
  }
  
  p <- length(sde_params)
  param_names <- names(sde_params)
  if (is.null(param_names)) stop("sde_params must be a named list")
  
  fix_idx <- which(param_names %in% fixpar)
  
  theta <- matrix(0, nrow = SGD_iter + 1, ncol = p)
  colnames(theta) <- param_names
  theta[1, ] <- unlist(sde_params)
  
  n <- nrow(data)
  Delta <- matrix(0, ncol = p, nrow = n - 1)
  I_p   <- diag(1, p)
  
  alpha_pot <- potential_params$alpha
  B         <- potential_params$B
  x_star    <- potential_params$x_star
  
  # Heating phase tracking (3rd-order exponential mean filter on gradient norm)
  heating_finished       <- FALSE
  K_end_heating          <- NA
  grad_norm_history      <- c()
  grad_norm_filtered      <- NA
  grad_norm_filtered_prev <- NA
  
  for (k in 1:SGD_iter) {
    if (verbose) message("SGD iter ", k, "\n")
    
    # Run particle filter
    filter <- particle_filter2D_cpp(
      observations     = as.matrix(data[, c("time", "Y1", "Y2")]),
      sde_params       = as.list(theta[k, ]),
      potential_params = potential_params,
      error_params     = error_params,
      error_dist       = error_dist,
      polygon_coords   = polygon@coords,
      U0 = U0, lambda = lambda,
      num_particles    = num_particles,
      split_around_fixed_point = split_around_fixed_point,
      scheme = scheme, ESS_threshold = 0.6,
      proposal_weight = 0.5,
      verbose = FALSE, print_timing = FALSE
    )
    
    if (smoothing) {
      backward_samples <- forward_filtering_backward_sampling(
        data, 1, forward_filter = filter,
        as.list(theta[k, ]), potential_params,
        error_params, error_dist, polygon,
        U0, lambda, num_particles,
        scheme = scheme,
        split_around_fixed_point = FALSE,
        verbose = FALSE
      )
      z <- t(apply(backward_samples, c(2, 3), mean))
    } else {
      z <- get_PF_mean(filter)
    }
    colnames(z) <- c("X1", "X2", "V1", "V2")
    
    # Compute push and potential gradient
    push_mat           <- matrix(NA, nrow = n, ncol = 2)
    potential_grad_mat <- matrix(NA, nrow = n, ncol = 2)
    
    for (t in 1:n) {
      X <- z[t, c("X1", "X2")]
      push_mat[t, ]           <- compute_push(X, polygon, lambda)
      potential_grad_mat[t, ] <- mix_gaussian_grad_cpp(
        X, x_star, list(B = B, alpha = alpha_pot), exclude = integer(0)
      )
    }
    
    # Compute gradient of pseudo-likelihood
    z_df      <- as.data.frame(z)
    z_df$time <- data$time
    
    grad_mat <- llk_gradient(
      data = z_df, push_mat = push_mat,
      potential_grad_mat = potential_grad_mat,
      tau = theta[k, "tau"], nu = theta[k, "nu"],
      omega = theta[k, "omega"], scheme = scheme,
      verbose = FALSE
    )
    
    v_k <- colSums(grad_mat)
    names(v_k) <- param_names
    if (length(fix_idx) > 0) v_k[fix_idx] <- 0
    
    # --- Three-phase learning rate schedule (Baey et al. 2023, Section 3.4.1) ---
    if (k < K_preheat) {
      # Phase 1: exponential growth from gamma0 to 1
      gamma_k <- gamma0^(1 - k / K_preheat)
    } else if (!heating_finished) {
      # Phase 2: heating — hold gamma at 1
      gamma_k <- 1
    } else {
      # Phase 3: Robbins-Monro decreasing schedule
      gamma_k <- (k - K_end_heating)^(-alpha)
    }
    if (verbose) message("gamma_k: ", gamma_k, "\n")
    
    # Stochastic approximation update of Delta and Fisher matrix estimate
    Delta <- (1 - gamma_k) * Delta + gamma_k * grad_mat
    I_k   <- matrix(0, ncol = p, nrow = p)
    for (i in 1:nrow(Delta)) {
      I_k <- I_k + Delta[i, ] %*% t(Delta[i, ])
    }
    
    # Preconditioner: stabilised during pre-heating, raw FIM estimate after
    if (k < K_preheat) {
      r_k <- max(1, sum(diag(I_k)))
      P_k <- (1 - gamma_k) * r_k * I_p + gamma_k * I_k
    } else {
      P_k <- I_k
    }
    
    # Compute step direction — regularise if P_k is singular
    direction <- tryCatch(
      solve(P_k, v_k),
      error = function(e) {
        if (verbose) cat("solve failed; regularizing matrix\n")
        P_k_reg <- 0.5 * (P_k + t(P_k)) + 1e-4 * diag(p)
        return(solve(P_k_reg, v_k))
      }
    )
    direction <- as.numeric(direction)
    if (length(fix_idx) > 0) direction[fix_idx] <- 0
    theta[k + 1, ] <- theta[k, ] + gamma_k * direction
    
    # --- Heating phase termination via 3rd-order exponential mean filter ---
    if (k >= K_preheat && !heating_finished) {
      grad_norm_now     <- sqrt(sum(v_k^2))
      grad_norm_history <- c(grad_norm_history, grad_norm_now)
      
      if (length(grad_norm_history) >= 3) {
        if (is.na(grad_norm_filtered)) {
          # Initialise filter from mean of first 3 values
          grad_norm_filtered <- mean(tail(grad_norm_history, 3))
        } else {
          grad_norm_filtered_prev <- grad_norm_filtered
          grad_norm_filtered      <- (1 - C_heating) * grad_norm_filtered +
            C_heating * grad_norm_now
          
          # End heating when filtered norm is no longer decreasing
          if (!is.na(grad_norm_filtered_prev) &&
              grad_norm_filtered >= grad_norm_filtered_prev) {
            heating_finished <- TRUE
            K_end_heating    <- k
            if (verbose) message("Heating phase ended at iteration ", k, "\n")
          }
        }
      }
    }
    
    if (verbose) {
      message("v_k: ",     paste(round(v_k, 4), collapse = " "), "\n")
      message("I_k: ",     paste(round(I_k, 4), collapse = " "), "\n")
      message("eigenvalues I_k: ",
              paste(round(eigen(I_k)$values, 4), collapse = " "), "\n")
      message("direction: ",
              paste(round(direction, 4), collapse = " "), "\n")
      message("theta: ",
              paste(round(theta[k + 1, ], 4), collapse = " "), "\n")
    }
  }
  
  return(list(
    theta         = theta,
    Delta         = Delta,
    K_end_heating = K_end_heating
  ))
}