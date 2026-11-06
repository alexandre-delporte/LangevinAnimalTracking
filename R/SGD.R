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
#'   optimisation (their gradients and updates are zeroed out). Applies to both
#'   movement parameters (\code{tau}, \code{nu}, \code{omega}) and, when
#'   \code{estimate_potential_params = TRUE}, potential parameters
#'   (\code{alpha_1}, \code{logL11_1}, \code{L21_1}, \code{logL22_1}, ...) where
#'   \code{logL11_k = log(L_k[1,1])} and \code{logL22_k = log(L_k[2,2])} are the
#'   log-diagonal entries of the lower Cholesky factor of \eqn{B_k}.
#'   Default \code{NULL} means all active parameters are estimated.
#' @param estimate_potential_params Logical. If \code{TRUE}, the potential
#'   parameters (\eqn{\alpha_k} and the log-diagonal Cholesky parametrisation of
#'   \eqn{B_k}) are jointly estimated alongside the movement parameters.
#'   Diagonal entries are stored as \code{logL11_k = log(L_k[1,1])} and
#'   \code{logL22_k = log(L_k[2,2])}, keeping \eqn{B_k} positive definite.
#'   By default the centres \eqn{x_k^*} are held fixed; see
#'   \code{estimate_centers} to also estimate them. Individual parameters can
#'   be held fixed via \code{fixpar}. Default \code{FALSE}.
#' @param estimate_centers Logical. If \code{TRUE}, the attraction centres
#'   \eqn{x_k^*} of the potential are also estimated, alongside \eqn{\alpha_k}
#'   and \eqn{B_k} (columns \code{xstar_1_1}, \code{xstar_1_2}, ...,
#'   \code{xstar_J_1}, \code{xstar_J_2} in the output). Requires
#'   \code{estimate_potential_params = TRUE}. This gives three increasingly
#'   flexible estimation modes: movement parameters only (both flags
#'   \code{FALSE}); movement parameters plus potential shape with fixed
#'   centres (\code{estimate_potential_params = TRUE},
#'   \code{estimate_centers = FALSE}); and the full potential including
#'   centres (both \code{TRUE}). Default \code{FALSE}.
#' @param SGD_iter A positive integer giving the total number of SGD
#'   iterations to run (across all three phases). Default \code{2000}.
#' @param potential_params A named list of parameters for the mixture of
#'   Gaussian potential \eqn{H}, with elements:
#'   \describe{
#'     \item{\code{alpha}}{Numeric vector of mixture weights.}
#'     \item{\code{B}}{List of precision matrices for each Gaussian component.}
#'     \item{\code{x_star}}{Matrix of attraction centres (one row per
#'       component).}
#'   }
#'   Default \code{NULL} means no potential (\eqn{H \equiv 0}), encoded as a
#'   single component with zero weight, identity precision matrix and centre
#'   \eqn{(0, 0)}. Must be provided when
#'   \code{estimate_potential_params = TRUE}.
#' @param error_dist A character string specifying the observation error
#'   distribution: \code{"normal"}, \code{"scaled_t"} or \code{"argos"}.
#' @param error_params A named list of parameters for the observation error
#'   distribution (e.g. \code{list(sigma_obs = 0.2)} for \code{"normal"}).
#'   Can be \code{NULL} only when \code{error_params_per_obs} is provided.
#' @param scheme A character string specifying the numerical splitting scheme
#'   used to approximate the SDE transition density. Either
#'   \code{"Lie-Trotter"} (default) or \code{"Strang"}.
#' @param polygon An \code{sp::Polygon} object defining the spatial domain
#'   \eqn{\mathcal{D}} within which animal movement is constrained. Can be
#'   \code{NULL} (the default) only when \code{lambda = Inf}, i.e. when there
#'   is no domain penalisation.
#' @param U0 A numeric vector of length 4 giving the initial state
#'   \eqn{(x_1, x_2, v_1, v_2)} (position and velocity) at time 0.
#' @param lambda A positive scalar giving the domain penalisation parameter
#'   \eqn{\lambda} used in the reflected Langevin push. Default \code{Inf}
#'   (no penalisation).
#' @param num_particles A positive integer giving the number of particles
#'   used in the particle filter at the start of SGD (iteration 1). Used for
#'   both the forward particle filter/FFBS and the conditional particle
#'   filter (when \code{smoothing_method = "CPF"}).
#' @param num_particles_min A positive integer giving the floor of the
#'   particle count, reached at iteration \code{num_particles_decay_iter} and
#'   held constant afterwards. Number of particles decreases linearly from
#'   \code{num_particles} to \code{num_particles_min} across iterations. For
#'   CPF, this floor can be as low as 5-10 since the conditional particle
#'   filter with ancestor sampling mixes well with very few particles once
#'   warm-started. Default is \code{num_particles} (i.e. no decay).
#' @param num_particles_decay_iter A positive integer giving the iteration at
#'   which the linear decay of the particle count reaches
#'   \code{num_particles_min}. Default \code{SGD_iter} (decay over the whole
#'   run).
#' @param split_around_fixed_point Logical. If \code{TRUE}, the potential is
#'   split around its fixed point. Not implemented yet: the gradient of the
#'   log-likelihood would require differentiating the matrix exponential of
#'   the drift around the fixed point w.r.t. the movement and potential
#'   parameters, so this throws an error. Default \code{FALSE}.
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
#' @param alpha_decay A scalar in \eqn{(1/2, 1]} controlling the decay rate of the
#'   learning rate in the decreasing phase: \eqn{\gamma_k = (k -
#'   K_{\text{end-heating}})^{-\alpha_{\text{decay}}}}. Must satisfy \eqn{\alpha_{\text{decay}} > 1/2} to
#'   ensure \eqn{\sum \gamma_k^2 < \infty}. Default \code{2/3}.
#' @param C_heating A small positive scalar giving the smoothing constant of
#'   the 3rd-order exponential mean filter (three exponential means in
#'   cascade) applied to the gradient vector during the heating phase. The
#'   heating phase ends when the norm of the filtered gradient stops
#'   decreasing. The authors recommend \eqn{1/1000}. Default \code{1/100}.
#' @param n_smooth_samples A non-negative integer. If greater than 0, the
#'   smoothed trajectories from the last \code{n_smooth_samples} iterations
#'   are saved and returned. These can be used to reconstruct plausible latent
#'   tracks consistent with the observations. Default \code{0}.
#' @param error_params_per_obs Optional list of length \code{N} of per-observation
#'   error parameter lists, as produced by \code{\link{make_argos_obs_params}}.
#'   Passed through to \code{particle_filter2D_cpp} and
#'   \code{forward_filtering_backward_sampling} at every SGD iteration.
#'   Only used when \code{error_dist = "argos"}. Default \code{NULL}.
#' @param smoothing_method A character string, either \code{"FFBS"}
#'   or \code{"CPF"} (default). \code{"FFBS"} draws a fresh smoothed trajectory at every
#'   SGD iteration via forward-filtering backward-sampling (\code{\link{forward_filtering_backward_sampling}}).
#'   \code{"CPF"} instead uses the conditional particle filter with ancestor
#'   sampling (\code{\link{conditional_particle_filter_cpp}}, Svensson et al.
#'   2015) as a persistent MCMC chain over trajectories: the trajectory
#'   sampled at iteration \code{k} is used as the conditioning
#'   \code{reference_trajectory} at iteration \code{k+1} (Markovian stochastic
#'   approximation), which typically mixes well with far fewer particles than
#'   FFBS.
#' @param n_sweeps_cpf A positive integer giving the number of CPF-AS sweeps
#'   run per SGD iteration when \code{smoothing_method = "CPF"}. Default
#'   \code{1} (one sweep per iteration, relying on warm-starting across SGD
#'   iterations for mixing).
#' @param progress_every A non-negative integer. If greater than 0, a one-line
#'   progress message (iteration, phase of the step-size schedule, elapsed
#'   time and current parameters) is printed every \code{progress_every}
#'   iterations and at the last iteration, independently of \code{verbose}.
#'   Useful to follow long runs in a log file. Default \code{0} (no progress
#'   messages).
#' @param label An optional character string prepended to the progress
#'   messages, e.g. \code{"seed 42"}, to tell apart replicates run in
#'   parallel. Default \code{NULL}.
#' @param init_trajectory An optional \code{n x 4} matrix giving the initial
#'   conditioning trajectory for the CPF-AS chain (columns \code{X1, X2, V1,
#'   V2}). Only used when \code{smoothing_method = "CPF"}. If \code{NULL}
#'   (default), an initial trajectory is drawn via a single FFBS sample before
#'   the first SGD iteration.
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
#'   \item{\code{smooth_samples}}{A list of length \code{n_smooth_samples}
#'     containing the smoothed trajectories from the last iterations. Each
#'     element is an \code{n x 4} matrix with columns
#'     \code{c("X1", "X2", "V1", "V2")}. \code{NULL} if
#'     \code{n_smooth_samples = 0}.}
#' }
#'
#' @export
SGD_Fisher <- function(data, sde_params, fixpar = NULL, SGD_iter = 2000,
                       potential_params = NULL, estimate_potential_params = FALSE,
                       estimate_centers = FALSE,
                       error_dist, error_params = NULL,
                       scheme = "Lie-Trotter", polygon = NULL, U0, lambda = Inf,
                       num_particles, num_particles_min = num_particles,
                       num_particles_decay_iter = SGD_iter,
                       split_around_fixed_point = FALSE,
                       verbose = FALSE, gamma0 = 1e-4, K_preheat = 1000,
                       alpha_decay = 2/3, C_heating = 1/100, n_smooth_samples = 0,
                       error_params_per_obs = NULL, smoothing_method = "CPF",
                       n_sweeps_cpf = 1,
                       init_trajectory = NULL,
                       progress_every = 0, label = NULL) {

  if (split_around_fixed_point) {
    stop("Not implemented yet with splitting around fixed point")
  }

  if (is.null(error_params) && is.null(error_params_per_obs)) {
    stop("error_params must be provided when error_params_per_obs is NULL")
  }

  if (estimate_centers && !estimate_potential_params) {
    stop("estimate_centers = TRUE requires estimate_potential_params = TRUE")
  }

  if (is.null(potential_params)) {
    if (estimate_potential_params) {
      stop("potential_params must be provided when estimate_potential_params = TRUE")
    }
    # No potential: single component with zero weight, so H and its gradient
    # are identically zero
    potential_params <- list(x_star = matrix(0, 1, 2), B = list(diag(2)), alpha = 0)
  }

  smoothing_method <- match.arg(smoothing_method, c("FFBS", "CPF"))

  if (is.null(polygon) && is.finite(lambda)) {
    stop("polygon must be provided when lambda is finite")
  }
  # Empty coordinate matrix when there is no domain: the push is zero for
  # lambda = Inf and the coordinates are never read
  polygon_coords <- if (is.null(polygon)) matrix(numeric(0), 0, 2) else polygon@coords
  ref_traj_cpf <- init_trajectory

  x_star    <- potential_params$x_star
  param_names <- names(sde_params)
  if (is.null(param_names)) stop("sde_params must be a named list")

  if (estimate_potential_params) {
    J <- length(potential_params$alpha)
    xi_names <- c(paste0("alpha_", 1:J),
                  unlist(lapply(1:J, function(j) paste0(c("logL11","L21","logL22"), "_", j))))
    xi_init <- c(potential_params$alpha,
                 unlist(lapply(potential_params$B, function(Bk) {
                   Lk <- t(chol(Bk))
                   c(log(Lk[1, 1]), Lk[2, 1], log(Lk[2, 2]))
                 })))
    if (estimate_centers) {
      xstar_names <- unlist(lapply(1:J, function(j) paste0("xstar_", j, "_", 1:2)))
      xi_names <- c(xi_names, xstar_names)
      xi_init  <- c(xi_init, as.vector(t(x_star)))
    }
    param_names <- c(param_names, xi_names)
    theta_init <- c(unlist(sde_params), xi_init)
  } else {
    J <- 0
    theta_init <- unlist(sde_params)
  }

  p       <- length(param_names)
  fix_idx <- which(param_names %in% fixpar)

  theta <- matrix(0, nrow = SGD_iter + 1, ncol = p)
  colnames(theta) <- param_names
  theta[1, ] <- theta_init

  n   <- nrow(data)
  Delta <- matrix(0, ncol = p, nrow = n - 1)
  I_p   <- diag(1, p)
  
  # Heating phase tracking (3rd-order exponential mean filter on gradient vector)
  heating_finished        <- FALSE
  K_end_heating           <- NA
  grad_filter             <- NULL
  grad_norm_filtered      <- NA
  grad_norm_filtered_prev <- NA
  
  # Storage for the last n_smooth_samples trajectories
  smooth_samples <- list()
  
  t_start <- Sys.time()

  for (k in 1:SGD_iter) {
    if (verbose) message("SGD iter ", k, "\n")

    # Linearly decay the particle count from num_particles (k = 1) to
    # num_particles_min (k >= num_particles_decay_iter)
    frac_k <- min(1, (k - 1) / max(1, num_particles_decay_iter - 1))
    num_particles_k <- round(num_particles + frac_k * (num_particles_min - num_particles))
    if (verbose) message("num_particles_k: ", num_particles_k, "\n")

    sde_params_k <- as.list(theta[k, c("tau", "nu", "omega")])

    if (estimate_potential_params) {
      alpha_k <- theta[k, paste0("alpha_", 1:J)]
      B_k <- lapply(1:J, function(j) {
        l  <- theta[k, paste0(c("logL11", "L21", "logL22"), "_", j)]
        Lj <- matrix(c(exp(l[1]), l[2], 0, exp(l[3])), 2, 2)
        Lj %*% t(Lj)
      })
      if (estimate_centers) {
        x_star_k <- matrix(NA_real_, nrow = J, ncol = 2)
        for (j in 1:J) {
          x_star_k[j, ] <- theta[k, paste0("xstar_", j, "_", 1:2)]
        }
      } else {
        x_star_k <- x_star
      }
      potential_params_k <- list(alpha = alpha_k, B = B_k, x_star = x_star_k)
    } else {
      potential_params_k <- potential_params
      x_star_k <- x_star
    }

    if (smoothing_method == "FFBS") {

      # Run particle filter
      filter <- particle_filter2D_cpp(
        observations     = as.matrix(data[, c("time", "Y1", "Y2")]),
        sde_params       = sde_params_k,
        potential_params = potential_params_k,
        error_params     = error_params,
        error_dist       = error_dist,
        polygon_coords   = polygon_coords,
        U0 = U0, lambda = lambda,
        num_particles    = num_particles_k,
        split_around_fixed_point = split_around_fixed_point,
        scheme = scheme, ESS_threshold = 0.8,
        verbose = FALSE, print_timing = FALSE,
        error_params_per_obs = error_params_per_obs
      )

      # Obtain latent trajectory via forward-filtering backward-sampling (FFBS)
      backward_samples <- forward_filtering_backward_sampling(
        data, 1, forward_filter = filter,
        sde_params_k, potential_params_k,
        error_params, error_dist, polygon,
        U0, lambda, num_particles_k,
        scheme = scheme,
        split_around_fixed_point = FALSE,
        verbose = FALSE,
        error_params_per_obs = error_params_per_obs
      )
      z <- t(apply(backward_samples, c(2, 3), mean))
      colnames(z) <- c("X1", "X2", "V1", "V2")

    } else { # smoothing_method == "CPF"

      if (is.null(ref_traj_cpf)) {
        # Seed the CPF-AS chain with a single FFBS draw at the initial parameters
        filter0 <- particle_filter2D_cpp(
          observations     = as.matrix(data[, c("time", "Y1", "Y2")]),
          sde_params       = sde_params_k,
          potential_params = potential_params_k,
          error_params     = error_params,
          error_dist       = error_dist,
          polygon_coords   = polygon_coords,
          U0 = U0, lambda = lambda,
          num_particles    = num_particles_k,
          split_around_fixed_point = split_around_fixed_point,
          scheme = scheme, ESS_threshold = 0.8,
          verbose = FALSE, print_timing = FALSE,
          error_params_per_obs = error_params_per_obs
        )
        seed_sample <- forward_filtering_backward_sampling(
          data, 1, forward_filter = filter0,
          sde_params_k, potential_params_k,
          error_params, error_dist, polygon,
          U0, lambda, num_particles_k,
          scheme = scheme,
          split_around_fixed_point = FALSE,
          verbose = FALSE,
          error_params_per_obs = error_params_per_obs
        )
        ref_traj_cpf <- t(seed_sample[1, , ])
      }

      cpf_out <- conditional_particle_filter_cpp(
        observations     = as.matrix(data[, c("time", "Y1", "Y2")]),
        sde_params       = sde_params_k,
        potential_params = potential_params_k,
        error_params     = error_params,
        error_dist       = error_dist,
        polygon_coords   = polygon_coords,
        U0 = U0, lambda = lambda,
        num_particles    = num_particles_k,
        scheme = scheme,
        split_around_fixed_point = split_around_fixed_point,
        reference_trajectory = ref_traj_cpf,
        n_sweeps = n_sweeps_cpf,
        ESS_threshold = 0.8,
        error_params_per_obs = error_params_per_obs
      )

      z <- cpf_out$trajectory
      ref_traj_cpf <- z  # warm-start next SGD iteration
      colnames(z) <- c("X1", "X2", "V1", "V2")
    }

    # Compute push and potential gradient
    push_mat           <- matrix(NA, nrow = n, ncol = 2)
    potential_grad_mat <- matrix(NA, nrow = n, ncol = 2)

    for (t in 1:n) {
      X <- z[t, c("X1", "X2")]
      push_mat[t, ]           <- compute_push_cpp(as.numeric(X), polygon_coords, lambda)
      potential_grad_mat[t, ] <- mix_gaussian_grad_cpp(
        X, x_star_k, list(B = potential_params_k$B, alpha = potential_params_k$alpha),
        exclude = integer(0)
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
      potential_params = if (estimate_potential_params) potential_params_k else NULL,
      estimate_centers = estimate_centers,
      verbose = FALSE
    )

    # A fixed parameter is treated as known, not estimated: its per-step
    # score is zeroed here, at the source, so neither the gradient nor the
    # running Fisher information estimate (Delta / I_k below, both built
    # from grad_mat) carry any contribution from it — including its
    # covariance with the free parameters, which a naive v_k[fix_idx] <- 0
    # after the fact would not remove from I_k.
    if (length(fix_idx) > 0) grad_mat[, fix_idx] <- 0

    v_k <- colSums(grad_mat)
    names(v_k) <- param_names

    # Three-phase learning rate schedule (Baey et al. 2023, Section 3.4.1) 
    if (k < K_preheat) {
      # Phase 1: exponential growth from gamma0 to 1
      gamma_k <- gamma0^(1 - k / K_preheat)
    } else if (!heating_finished) {
      # Phase 2: heating — keep gamma at 1
      gamma_k <- 1
    } else {
      # Phase 3: Robbins-Monro decrease
      gamma_k <- (k - K_end_heating)^(-alpha_decay)
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
    
    # Compute step direction over the free parameters only. fix_idx rows/
    # columns of P_k are now exactly zero (grad_mat was zeroed there above),
    # not just small, so they are dropped from the solve entirely rather
    # than left in and patched up by regularisation — otherwise every
    # fixed parameter would make P_k structurally singular from here on,
    # and its (regularised, near-zero) row/column would still formally
    # enter the linear system used to compute the free parameters' step.
    free_idx <- if (length(fix_idx) > 0) setdiff(seq_len(p), fix_idx) else seq_len(p)
    P_k_free <- P_k[free_idx, free_idx, drop = FALSE]
    v_k_free <- v_k[free_idx]
    direction_free <- tryCatch(
      solve(P_k_free, v_k_free),
      error = function(e) {
        if (verbose) cat("solve failed; regularizing matrix\n")
        P_k_free_reg <- 0.5 * (P_k_free + t(P_k_free)) + 1e-4 * diag(length(free_idx))
        return(solve(P_k_free_reg, v_k_free))
      }
    )
    direction <- numeric(p)
    direction[free_idx] <- as.numeric(direction_free)
    theta[k + 1, ] <- theta[k, ] + gamma_k * direction
    
    # --- Heating phase termination via 3rd-order exponential mean filter ---
    # (Baey et al. 2023, Algorithm 3): three exponential means in cascade
    # applied to the gradient vector v_k, the heating ends when the norm of
    # the filtered gradient stops decreasing
    if (k >= K_preheat && !heating_finished) {
      if (is.null(grad_filter)) {
        # Initialise all three stages at the first gradient
        grad_filter <- list(v_k, v_k, v_k)
      } else {
        grad_filter[[1]] <- (1 - C_heating) * grad_filter[[1]] + C_heating * v_k
        grad_filter[[2]] <- (1 - C_heating) * grad_filter[[2]] + C_heating * grad_filter[[1]]
        grad_filter[[3]] <- (1 - C_heating) * grad_filter[[3]] + C_heating * grad_filter[[2]]
      }
      grad_norm_filtered_prev <- grad_norm_filtered
      grad_norm_filtered      <- sqrt(sum(grad_filter[[3]]^2))

      # End heating when filtered gradient norm is no longer decreasing
      if (!is.na(grad_norm_filtered_prev) &&
          grad_norm_filtered >= grad_norm_filtered_prev) {
        heating_finished <- TRUE
        K_end_heating    <- k
        if (verbose) message("Heating phase ended at iteration ", k, "\n")
      }
    }
    
    if (progress_every > 0 && (k %% progress_every == 0 || k == SGD_iter)) {
      phase <- if (k < K_preheat) "pre-heating" else if (!heating_finished) "heating" else "decreasing"
      sgd_progress_message(label, k, SGD_iter, phase, theta[k + 1, ], t_start)
    }

    # Save the trajectory from this iteration if we are in the last
    # n_smooth_samples iterations
    if (n_smooth_samples > 0 && k > SGD_iter - n_smooth_samples) {
      smooth_samples <- c(smooth_samples, list(z))
    }
    
    if (verbose) {
      
      ev <- eigen(I_k)
      
      message("v_k: ",     paste(round(v_k, 4), collapse = " "), "\n")
      message("I_k: ",     paste(round(I_k, 4), collapse = " "), "\n")
      message("eigenvalues I_k: ",
              paste(round(ev$values,4), collapse = " "), "\n")
      message("direction: ",
              paste(round(direction, 4), collapse = " "), "\n")
      message("theta: ",
              paste(round(theta[k + 1, ], 4), collapse = " "), "\n")
      
      null_dir<- ev$vectors[,which.min(ev$values)]
      names(null_dir) <- param_names
      message("null direction: ",paste(round(null_dir, 4), collapse = " "), "\n")
    }
  }
  
  return(list(
    theta          = theta,
    Delta          = Delta,
    K_end_heating  = K_end_heating,
    smooth_samples = if (n_smooth_samples > 0) smooth_samples else NULL
  ))
}

#' Stochastic Gradient Ascent with Adam for SDE parameter estimation
#'
#' Same stochastic approximation scheme as \code{\link{SGD_Fisher}} (a latent
#' trajectory is drawn at each iteration by FFBS or CPF-AS, and the gradient
#' of the complete-data log-likelihood is computed on it), but the parameters
#' are updated with the Adam optimiser (Kingma and Ba, 2015) instead of the
#' Fisher-preconditioned step. The Fisher information matrix is still estimated
#' across iterations, with the same per-step score stochastic approximation
#' \eqn{I = \Delta^\top \Delta} as in \code{SGD_Fisher}, but it never enters
#' the updates, so a poorly estimated (near singular) FIM cannot make the
#' algorithm diverge. It is returned for computing standard errors at the end.
#'
#' The Adam step size schedule has two phases: a constant phase at \code{lr},
#' followed by a Robbins-Monro decreasing phase
#' \eqn{lr_k = lr\,(k - K_{\text{decay}} + 1)^{-\alpha_{\text{decay}}}}. The
#' start of the decreasing phase \eqn{K_{\text{decay}}} is either given by
#' \code{K_decay}, or detected automatically, after at least
#' \code{K_constant_min} iterations, as the iteration at which the norm of the
#' 3rd-order exponentially filtered gradient stops decreasing (same criterion
#' as the end of the heating phase in \code{SGD_Fisher}).
#'
#' Since Adam normalises each coordinate of the gradient, the step taken on
#' each parameter is of order \code{lr} in the units of that parameter, whatever
#' the scale of the gradient. \code{lr} should therefore be chosen relative to
#' the scale of the parameters, and can be given per parameter.
#'
#' @inheritParams SGD_Fisher
#' @param lr A positive scalar, or a named numeric vector with one entry per
#'   estimated parameter (names as in the columns of the output \code{theta},
#'   e.g. \code{tau}, \code{nu}, \code{omega}, \code{alpha_1}, \code{logL11_1},
#'   ...), giving the Adam step size during the constant phase. Default
#'   \code{0.01}.
#' @param beta1,beta2 Exponential decay rates of the Adam first and second
#'   moment estimates. Defaults \code{0.9} and \code{0.999}.
#' @param eps Small constant added to the denominator of the Adam step.
#'   Default \code{1e-8}.
#' @param K_constant_min A positive integer giving the minimum number of
#'   iterations at constant step size before the automatic detection of the
#'   start of the decreasing phase is switched on. Ignored if \code{K_decay}
#'   is given. Default \code{200}.
#' @param K_decay Optional positive integer giving the iteration at which the
#'   decreasing phase starts. If \code{NULL} (default), it is detected
#'   automatically with the filtered gradient criterion.
#' @param alpha_decay A scalar in \eqn{(1/2, 1]} giving the decay rate of the
#'   step size in the decreasing phase. Default \code{2/3}.
#' @param C_heating Smoothing constant of the 3rd-order exponential mean
#'   filter applied to the gradient to detect the start of the decreasing
#'   phase. Default \code{1/100}.
#' @param fisher_rate A scalar in \eqn{(0, 1]} giving the weight of the
#'   stochastic approximation of the per-step scores \eqn{\Delta} during the
#'   constant phase (exponential moving average with memory of about
#'   \code{1/fisher_rate} iterations). During the decreasing phase the weight
#'   decreases as \eqn{1/(k - K_{\text{decay}} + 1/\text{fisher\_rate})}, so
#'   that \eqn{\Delta} becomes a running average of the scores over the
#'   decreasing phase. Default \code{0.05}.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{theta}}{A matrix of dimension \code{(SGD_iter + 1)}
#'     \eqn{\times} \code{p} of the parameter estimates at each iteration.}
#'   \item{\code{Delta}}{The \code{(n - 1)} \eqn{\times} \code{p} matrix of
#'     the stochastic approximation of the per-step scores at the final
#'     iteration.}
#'   \item{\code{fisher}}{The final estimate of the Fisher information
#'     matrix \eqn{\Delta^\top \Delta} (rows and columns of fixed parameters
#'     are zero).}
#'   \item{\code{fisher_min_eigen}}{A vector of length \code{SGD_iter} with
#'     the smallest eigenvalue of the FIM estimate (restricted to the free
#'     parameters) at each iteration, to monitor its stabilisation.}
#'   \item{\code{grad}}{A matrix of dimension \code{SGD_iter} \eqn{\times}
#'     \code{p} with the stochastic gradient at each iteration.}
#'   \item{\code{K_decay}}{The iteration at which the decreasing phase
#'     started, \code{NA} if it had not started when the algorithm
#'     terminated.}
#'   \item{\code{smooth_samples}}{As in \code{\link{SGD_Fisher}}.}
#' }
#'
#' @references Kingma, D. P. and Ba, J. (2015). Adam: A method for
#'   stochastic optimization. ICLR.
#'
#' @export
SGD_adam <- function(data, sde_params, fixpar = NULL, SGD_iter = 2000,
                     potential_params = NULL, estimate_potential_params = FALSE,
                     estimate_centers = FALSE,
                     error_dist, error_params = NULL,
                     scheme = "Lie-Trotter", polygon = NULL, U0, lambda = Inf,
                     num_particles, num_particles_min = num_particles,
                     num_particles_decay_iter = SGD_iter,
                     split_around_fixed_point = FALSE,
                     verbose = FALSE, lr = 0.01, beta1 = 0.9, beta2 = 0.999,
                     eps = 1e-8, K_constant_min = 200, K_decay = NULL,
                     alpha_decay = 2/3, C_heating = 1/100, fisher_rate = 0.05,
                     n_smooth_samples = 0,
                     error_params_per_obs = NULL, smoothing_method = "CPF",
                     n_sweeps_cpf = 1,
                     init_trajectory = NULL,
                     progress_every = 0, label = NULL) {

  if (split_around_fixed_point) {
    stop("Not implemented yet with splitting around fixed point")
  }

  if (is.null(error_params) && is.null(error_params_per_obs)) {
    stop("error_params must be provided when error_params_per_obs is NULL")
  }

  if (estimate_centers && !estimate_potential_params) {
    stop("estimate_centers = TRUE requires estimate_potential_params = TRUE")
  }

  if (is.null(potential_params)) {
    if (estimate_potential_params) {
      stop("potential_params must be provided when estimate_potential_params = TRUE")
    }
    potential_params <- list(x_star = matrix(0, 1, 2), B = list(diag(2)), alpha = 0)
  }

  smoothing_method <- match.arg(smoothing_method, c("FFBS", "CPF"))

  if (is.null(polygon) && is.finite(lambda)) {
    stop("polygon must be provided when lambda is finite")
  }
  polygon_coords <- if (is.null(polygon)) matrix(numeric(0), 0, 2) else polygon@coords
  ref_traj_cpf <- init_trajectory

  x_star    <- potential_params$x_star
  param_names <- names(sde_params)
  if (is.null(param_names)) stop("sde_params must be a named list")

  if (estimate_potential_params) {
    J <- length(potential_params$alpha)
    xi_names <- c(paste0("alpha_", 1:J),
                  unlist(lapply(1:J, function(j) paste0(c("logL11","L21","logL22"), "_", j))))
    xi_init <- c(potential_params$alpha,
                 unlist(lapply(potential_params$B, function(Bk) {
                   Lk <- t(chol(Bk))
                   c(log(Lk[1, 1]), Lk[2, 1], log(Lk[2, 2]))
                 })))
    if (estimate_centers) {
      xstar_names <- unlist(lapply(1:J, function(j) paste0("xstar_", j, "_", 1:2)))
      xi_names <- c(xi_names, xstar_names)
      xi_init  <- c(xi_init, as.vector(t(x_star)))
    }
    param_names <- c(param_names, xi_names)
    theta_init <- c(unlist(sde_params), xi_init)
  } else {
    J <- 0
    theta_init <- unlist(sde_params)
  }

  p        <- length(param_names)
  fix_idx  <- which(param_names %in% fixpar)
  free_idx <- setdiff(seq_len(p), fix_idx)

  # Per-parameter base step size
  if (length(lr) == 1 && is.null(names(lr))) {
    lr_vec <- rep(lr, p)
  } else {
    missing_lr <- setdiff(param_names[free_idx], names(lr))
    if (length(missing_lr) > 0) {
      stop("lr is missing entries for: ", paste(missing_lr, collapse = ", "))
    }
    lr_vec <- rep(0, p)
    lr_vec[free_idx] <- lr[param_names[free_idx]]
  }
  names(lr_vec) <- param_names

  theta <- matrix(0, nrow = SGD_iter + 1, ncol = p)
  colnames(theta) <- param_names
  theta[1, ] <- theta_init

  grad_hist <- matrix(NA_real_, nrow = SGD_iter, ncol = p)
  colnames(grad_hist) <- param_names
  fisher_min_eigen <- rep(NA_real_, SGD_iter)

  n     <- nrow(data)
  Delta <- matrix(0, ncol = p, nrow = n - 1)

  # Adam moments
  m_adam <- numeric(p)
  v_adam <- numeric(p)

  # Decreasing phase detection (3rd-order exponential mean filter on gradient)
  decay_started           <- FALSE
  K_decay_used            <- NA
  grad_filter             <- NULL
  grad_norm_filtered      <- NA
  grad_norm_filtered_prev <- NA

  smooth_samples <- list()

  t_start <- Sys.time()

  for (k in 1:SGD_iter) {
    if (verbose) message("SGD iter ", k, "\n")

    frac_k <- min(1, (k - 1) / max(1, num_particles_decay_iter - 1))
    num_particles_k <- round(num_particles + frac_k * (num_particles_min - num_particles))

    sde_params_k <- as.list(theta[k, c("tau", "nu", "omega")])

    if (estimate_potential_params) {
      alpha_k <- theta[k, paste0("alpha_", 1:J)]
      B_k <- lapply(1:J, function(j) {
        l  <- theta[k, paste0(c("logL11", "L21", "logL22"), "_", j)]
        Lj <- matrix(c(exp(l[1]), l[2], 0, exp(l[3])), 2, 2)
        Lj %*% t(Lj)
      })
      if (estimate_centers) {
        x_star_k <- matrix(NA_real_, nrow = J, ncol = 2)
        for (j in 1:J) {
          x_star_k[j, ] <- theta[k, paste0("xstar_", j, "_", 1:2)]
        }
      } else {
        x_star_k <- x_star
      }
      potential_params_k <- list(alpha = alpha_k, B = B_k, x_star = x_star_k)
    } else {
      potential_params_k <- potential_params
      x_star_k <- x_star
    }

    if (smoothing_method == "FFBS") {

      filter <- particle_filter2D_cpp(
        observations     = as.matrix(data[, c("time", "Y1", "Y2")]),
        sde_params       = sde_params_k,
        potential_params = potential_params_k,
        error_params     = error_params,
        error_dist       = error_dist,
        polygon_coords   = polygon_coords,
        U0 = U0, lambda = lambda,
        num_particles    = num_particles_k,
        split_around_fixed_point = split_around_fixed_point,
        scheme = scheme, ESS_threshold = 0.8,
        verbose = FALSE, print_timing = FALSE,
        error_params_per_obs = error_params_per_obs
      )

      backward_samples <- forward_filtering_backward_sampling(
        data, 1, forward_filter = filter,
        sde_params_k, potential_params_k,
        error_params, error_dist, polygon,
        U0, lambda, num_particles_k,
        scheme = scheme,
        split_around_fixed_point = FALSE,
        verbose = FALSE,
        error_params_per_obs = error_params_per_obs
      )
      z <- t(apply(backward_samples, c(2, 3), mean))
      colnames(z) <- c("X1", "X2", "V1", "V2")

    } else { # smoothing_method == "CPF"

      if (is.null(ref_traj_cpf)) {
        filter0 <- particle_filter2D_cpp(
          observations     = as.matrix(data[, c("time", "Y1", "Y2")]),
          sde_params       = sde_params_k,
          potential_params = potential_params_k,
          error_params     = error_params,
          error_dist       = error_dist,
          polygon_coords   = polygon_coords,
          U0 = U0, lambda = lambda,
          num_particles    = num_particles_k,
          split_around_fixed_point = split_around_fixed_point,
          scheme = scheme, ESS_threshold = 0.8,
          verbose = FALSE, print_timing = FALSE,
          error_params_per_obs = error_params_per_obs
        )
        seed_sample <- forward_filtering_backward_sampling(
          data, 1, forward_filter = filter0,
          sde_params_k, potential_params_k,
          error_params, error_dist, polygon,
          U0, lambda, num_particles_k,
          scheme = scheme,
          split_around_fixed_point = FALSE,
          verbose = FALSE,
          error_params_per_obs = error_params_per_obs
        )
        ref_traj_cpf <- t(seed_sample[1, , ])
      }

      cpf_out <- conditional_particle_filter_cpp(
        observations     = as.matrix(data[, c("time", "Y1", "Y2")]),
        sde_params       = sde_params_k,
        potential_params = potential_params_k,
        error_params     = error_params,
        error_dist       = error_dist,
        polygon_coords   = polygon_coords,
        U0 = U0, lambda = lambda,
        num_particles    = num_particles_k,
        scheme = scheme,
        split_around_fixed_point = split_around_fixed_point,
        reference_trajectory = ref_traj_cpf,
        n_sweeps = n_sweeps_cpf,
        ESS_threshold = 0.8,
        error_params_per_obs = error_params_per_obs
      )

      z <- cpf_out$trajectory
      ref_traj_cpf <- z
      colnames(z) <- c("X1", "X2", "V1", "V2")
    }

    push_mat           <- matrix(NA, nrow = n, ncol = 2)
    potential_grad_mat <- matrix(NA, nrow = n, ncol = 2)
    for (t in 1:n) {
      X <- z[t, c("X1", "X2")]
      push_mat[t, ]           <- compute_push_cpp(as.numeric(X), polygon_coords, lambda)
      potential_grad_mat[t, ] <- mix_gaussian_grad_cpp(
        X, x_star_k, list(B = potential_params_k$B, alpha = potential_params_k$alpha),
        exclude = integer(0)
      )
    }

    z_df      <- as.data.frame(z)
    z_df$time <- data$time

    grad_mat <- llk_gradient(
      data = z_df, push_mat = push_mat,
      potential_grad_mat = potential_grad_mat,
      tau = theta[k, "tau"], nu = theta[k, "nu"],
      omega = theta[k, "omega"], scheme = scheme,
      potential_params = if (estimate_potential_params) potential_params_k else NULL,
      estimate_centers = estimate_centers,
      verbose = FALSE
    )

    # Fixed parameters: zero their per-step scores at the source (see SGD_Fisher)
    if (length(fix_idx) > 0) grad_mat[, fix_idx] <- 0

    g_k <- colSums(grad_mat)
    names(g_k) <- param_names
    grad_hist[k, ] <- g_k

    # --- Start of the decreasing phase ---
    if (!decay_started) {
      if (!is.null(K_decay)) {
        if (k >= K_decay) {
          decay_started <- TRUE
          K_decay_used  <- k
        }
      } else if (k >= K_constant_min) {
        if (is.null(grad_filter)) {
          grad_filter <- list(g_k, g_k, g_k)
        } else {
          grad_filter[[1]] <- (1 - C_heating) * grad_filter[[1]] + C_heating * g_k
          grad_filter[[2]] <- (1 - C_heating) * grad_filter[[2]] + C_heating * grad_filter[[1]]
          grad_filter[[3]] <- (1 - C_heating) * grad_filter[[3]] + C_heating * grad_filter[[2]]
        }
        grad_norm_filtered_prev <- grad_norm_filtered
        grad_norm_filtered      <- sqrt(sum(grad_filter[[3]]^2))
        if (!is.na(grad_norm_filtered_prev) &&
            grad_norm_filtered >= grad_norm_filtered_prev) {
          decay_started <- TRUE
          K_decay_used  <- k
        }
      }
      if (decay_started && verbose) message("Decreasing phase started at iteration ", k, "\n")
    }

    # --- Step size schedule ---
    if (!decay_started) {
      lr_scale_k <- 1
      gamma_F_k  <- max(1 / k, fisher_rate)
    } else {
      lr_scale_k <- (k - K_decay_used + 1)^(-alpha_decay)
      gamma_F_k  <- 1 / (k - K_decay_used + 1 / fisher_rate)
    }

    # --- Fisher information estimate (not used in the update) ---
    Delta <- (1 - gamma_F_k) * Delta + gamma_F_k * grad_mat
    I_k   <- crossprod(Delta)
    fisher_min_eigen[k] <- min(eigen(I_k[free_idx, free_idx, drop = FALSE],
                                     symmetric = TRUE, only.values = TRUE)$values)

    # --- Adam update (gradient ascent on the log-likelihood) ---
    m_adam <- beta1 * m_adam + (1 - beta1) * g_k
    v_adam <- beta2 * v_adam + (1 - beta2) * g_k^2
    m_hat  <- m_adam / (1 - beta1^k)
    v_hat  <- v_adam / (1 - beta2^k)
    step   <- lr_scale_k * lr_vec * m_hat / (sqrt(v_hat) + eps)
    step[fix_idx] <- 0
    theta[k + 1, ] <- theta[k, ] + step

    if (progress_every > 0 && (k %% progress_every == 0 || k == SGD_iter)) {
      phase <- if (decay_started) "decreasing" else "constant"
      sgd_progress_message(label, k, SGD_iter, phase, theta[k + 1, ], t_start)
    }

    if (n_smooth_samples > 0 && k > SGD_iter - n_smooth_samples) {
      smooth_samples <- c(smooth_samples, list(z))
    }

    if (verbose) {
      message("grad: ",  paste(round(g_k, 4), collapse = " "), "\n")
      message("step: ",  paste(signif(step, 4), collapse = " "), "\n")
      message("min eigenvalue I_k: ", signif(fisher_min_eigen[k], 4), "\n")
      message("theta: ", paste(round(theta[k + 1, ], 4), collapse = " "), "\n")
    }
  }

  fisher <- crossprod(Delta)
  dimnames(fisher) <- list(param_names, param_names)

  return(list(
    theta            = theta,
    Delta            = Delta,
    fisher           = fisher,
    fisher_min_eigen = fisher_min_eigen,
    grad             = grad_hist,
    K_decay          = K_decay_used,
    smooth_samples   = if (n_smooth_samples > 0) smooth_samples else NULL
  ))
}


#' One-line progress message for the SGD algorithms
#'
#' Prints the label, iteration, phase of the step-size schedule, elapsed time
#' and current parameters (4 significant digits). Used by
#' \code{\link{SGD_Fisher}} and \code{\link{SGD_adam}} when
#' \code{progress_every > 0}.
#'
#' @noRd
sgd_progress_message <- function(label, k, SGD_iter, phase, theta_k, t_start) {
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  message(sprintf("%siter %d/%d (%s, %.1f min) theta: %s",
                  if (is.null(label)) "" else paste0("[", label, "] "),
                  k, SGD_iter, phase, elapsed,
                  paste(names(theta_k), signif(theta_k, 4), sep = "=", collapse = " ")))
}
