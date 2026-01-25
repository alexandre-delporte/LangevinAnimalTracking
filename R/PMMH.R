
#' Particle Marginal Metropolis-Hastings for penalised Langevin SDE . See Andrieu (2010)
#'  for the general algorithm
#' @param data A data frame containing observations with columns: time, Y1, Y2
#' @param theta_init named list with initial values for parameters tau, nu, omega
#' @param num_iterations number of PMMH iterations
#' @param num_particles number of particles to use in the particle filter
#' @param potential_params parameters for the potential gradient function
#' @param error_params parameters for the measurement error distribution
#' @param error_dist type of measurement error : "normal", "scaled_t" or "argos"
#' @param polygon Polygon defining the ecological domain of interest
#' @param U0 initial state (vector of length 4)
#' @param lambda penalty parameter for the inward push
#' @param proposal_cov proposal covariance matrix for the random walk on log-scale
#' @param prior_log_density function to compute log-density of the prior distribution
#' @param adapt_proposal logical. If TRUE, adapt the proposal covariance during burn-in
#' @param burnin number of burn-in iterations
#' @param fixpar character vector with names of parameters to keep fixed during the inference
#' @param verbose logical. If TRUE, print information during the computation
#' # @return A list with the following elements:
#' #' - theta: matrix of parameter values at each iteration
#' #' - log_likelihood: vector of log-likelihood values at each iteration
#' #' - acceptance_rate: overall acceptance rate
#' #' - proposal_cov: final proposal covariance matrix
#' #' - posterior_samples: matrix of posterior samples after burn-in
#'  - posterior_mean: posterior mean estimates for each parameter
#' - posterior_sd: posterior standard deviation estimates for each parameter
#' - posterior_quantiles: posterior quantiles (2.5%, 25%, 50%, 75%, 97.5%) for each parameter
#' - effective_sample_size: effective sample size for each parameter
#'
#' @importFrom MASS mvrnorm
#' @export

PMMH_langevin <- function(data, theta_init, num_iterations, num_particles,
                          potential_params = NULL,
                          error_params, error_dist,
                          polygon, U0, lambda,
                          prior_tau, prior_nu, prior_omega,
                          fixpar = NULL,
                          # AM-specific parameters
                          initial_cov = NULL,
                          t0 = 100,  # Initial period length
                          epsilon = 1e-5,  # Small constant for regularization
                          sd_scaling = NULL,  # Will use (2.4)^2/d by default
                          verbose = TRUE) {
  
  p <- length(theta_init)
  param_names <- names(theta_init)
  
  freepar <- setdiff(param_names, fixpar)
  p_free <- length(freepar)
  
  # Set default scaling parameter (Gelman et al. 1996)
  if (is.null(sd_scaling)) {
    sd_scaling <- (2.4)^2 / p_free
  }
  
  # Initialize proposal covariance
  if (is.null(initial_cov)) {
    # Default: diagonal matrix with reasonable initial variances
    initial_cov <- diag(c(tau = 0.1, nu = 0.1, omega = 0.1)[freepar])
  }
  
  # Storage
  theta_init<-unlist(theta_init)
  theta_chain <- matrix(0, nrow = num_iterations, ncol = p)
  colnames(theta_chain) <- param_names
  theta_chain[1, ] <- theta_init
  
  log_lik_chain <- numeric(num_iterations)
  accept_count <- 0
  
  # Transformed parameters
  theta_transformed_chain <- matrix(0, nrow = num_iterations, ncol = p_free)
  colnames(theta_transformed_chain) <- freepar
  
  # Transform initial parameters
  theta_trans <- numeric(p_free)
  names(theta_trans) <- freepar
  if ("tau" %in% freepar) theta_trans["tau"] <- log(theta_init["tau"])
  if ("nu" %in% freepar) theta_trans["nu"] <- log(theta_init["nu"])
  if ("omega" %in% freepar) theta_trans["omega"] <- theta_init["omega"]
  theta_transformed_chain[1, ] <- theta_trans
  
  # Current proposal covariance
  C_t <- initial_cov
  
  # Running sum for efficient covariance computation
  sum_theta <- theta_trans
  sum_theta_sq <- outer(theta_trans, theta_trans)
  
  # Log normal prior
  prior_log_ratio <- function(theta, theta_new) {
    if (theta["tau"] <= 0 || theta["nu"] <= 0) return(-Inf)
    
    log_tau_ratio <- -((log(theta_new["tau"]) - log(prior_tau$mean))^2 - 
                         (log(theta["tau"]) - log(prior_tau$mean))^2) / (2 * prior_tau$sd^2)
    log_nu_ratio <- -((log(theta_new["nu"]) - log(prior_nu$mean))^2 - 
                        (log(theta["nu"]) - log(prior_nu$mean))^2) / (2 * prior_nu$sd^2)
    log_omega_ratio <- -(theta_new["omega"]^2 - theta["omega"]^2) / (2 * prior_omega$sd^2)
    
    return(log_tau_ratio + log_nu_ratio + log_omega_ratio)
  }
  
  # Run initial particle filter
  if (verbose) cat("Running initial particle filter...\n")
  pf_init <- particle_filter2D(
    data = data,
    sde_params = as.list(theta_init),
    potential_params = potential_params,
    error_params = error_params,
    error_dist = error_dist,
    polygon = polygon,
    U0 = U0,
    lambda = lambda,
    num_particles = num_particles,
    scheme = "Lie-Trotter",
    split_around_fixed_point = FALSE,
    verbose = FALSE
  )
  
  log_lik_current <- pf_init$loglik
  log_lik_chain[1] <- log_lik_current
  
  if (verbose) {
    cat("Initial log-likelihood:", log_lik_current, "\n")
    cat("Using Adaptive Metropolis with sd_scaling =", sd_scaling, "\n")
    cat("Initial period: iterations 1 to", t0, "\n")
  }
  
  # MCMC loop
  for (iter in 2:num_iterations) {
    if (verbose && iter %% 100 == 0) {
      cat("\n========================================\n")
      cat("PMMH iteration", iter, "/", num_iterations, "\n")
      cat("  Current theta:", round(theta_chain[iter-1, ], 4), "\n")
      cat("  Current log-lik:", round(log_lik_current, 2), "\n")
      cat("  Acceptance rate:", round(accept_count / (iter-1), 3), "\n")
    }
    
    theta_current <- theta_chain[iter - 1, ]
    theta_trans_current <- theta_transformed_chain[iter - 1, ]
    
    # Update proposal covariance using Adaptive Metropolis
    if (iter > t0) {
      # Compute empirical covariance from all previous samples
      mean_theta <- sum_theta / (iter - 1)
      
      # Empirical covariance: E[XX^T] - E[X]E[X]^T
      emp_cov <- (sum_theta_sq / (iter - 1)) - outer(mean_theta, mean_theta)
      
      # Apply AM formula: C_t = s_d * Cov + s_d * epsilon * I
      C_t <- sd_scaling * emp_cov + sd_scaling * epsilon * diag(p_free)
      
      if (verbose && iter %% 500 == 0) {
        cat("  Updated proposal covariance (diag):", round(diag(C_t), 6), "\n")
      }
    }
    
    # Propose new parameters from multivariate normal
    # Use Cholesky decomposition for sampling
    L <- tryCatch(chol(C_t), error = function(e) NULL)
    
    if (is.null(L)) {
      # Fallback if covariance is not positive definite
      if (verbose) cat("  Warning: Proposal covariance not positive definite, using diagonal\n")
      L <- diag(sqrt(pmax(diag(C_t), epsilon)))
    }
    
    # Sample from N(theta_current, C_t) in transformed space
    z <- rnorm(p_free)
    theta_trans_prop <- theta_trans_current + as.vector(t(L) %*% z)
    names(theta_trans_prop) <- freepar
    
    # Transform back to original space
    theta_prop <- theta_current
    if ("tau" %in% freepar) theta_prop["tau"] <- exp(theta_trans_prop["tau"])
    if ("nu" %in% freepar) theta_prop["nu"] <- exp(theta_trans_prop["nu"])
    if ("omega" %in% freepar) theta_prop["omega"] <- theta_trans_prop["omega"]
    
    if (verbose && iter %% 100 == 0) {
      cat("  Proposed theta:", round(theta_prop, 4), "\n")
    }
    
    # Check for invalid proposals
    if (theta_prop["tau"] <= 0 || theta_prop["nu"] <= 0 || 
        any(is.na(theta_prop)) || any(is.infinite(theta_prop))) {
      if (verbose && iter %% 100 == 0) {
        cat("  REJECTED (invalid proposal)\n")
      }
      theta_chain[iter, ] <- theta_current
      theta_transformed_chain[iter, ] <- theta_trans_current
      log_lik_chain[iter] <- log_lik_current
      next
    }
    
    # Run particle filter with proposed parameters
    pf_prop <- tryCatch({
      particle_filter2D(
        data = data,
        sde_params = as.list(theta_prop),
        potential_params = potential_params,
        error_params = error_params,
        error_dist = error_dist,
        polygon = polygon,
        U0 = U0,
        lambda = lambda,
        num_particles = num_particles,
        scheme = "Lie-Trotter",
        split_around_fixed_point = FALSE,
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) cat("  Particle filter error:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(pf_prop)) {
      theta_chain[iter, ] <- theta_current
      theta_transformed_chain[iter, ] <- theta_trans_current
      log_lik_chain[iter] <- log_lik_current
      next
    }
    
    log_lik_prop <- pf_prop$loglik
    
    if (verbose && iter %% 100 == 0) {
      cat("  Proposed log-lik:", round(log_lik_prop, 2), "\n")
    }
    
    # Acceptance probability (no Jacobian needed - proposal is symmetric in transformed space)
    log_alpha <- (log_lik_prop - log_lik_current) + 
      prior_log_ratio(theta_current, theta_prop)
    
    if (verbose && iter %% 100 == 0) {
      cat("  Log acceptance prob:", round(log_alpha, 4), "\n")
    }
    
    if (!is.finite(log_alpha)) {
      accept <- FALSE
      if (verbose && iter %% 100 == 0) {
        cat("  REJECTED (non-finite log-alpha)\n")
      }
    } else {
      u <- runif(1)
      accept <- (log(u) < log_alpha)
      if (verbose && iter %% 100 == 0) {
        if (accept) {
          cat("  *** ACCEPTED ***\n")
        } else {
          cat("  REJECTED\n")
        }
      }
    }
    
    if (accept) {
      theta_chain[iter, ] <- theta_prop
      theta_transformed_chain[iter, ] <- theta_trans_prop
      log_lik_current <- log_lik_prop
      accept_count <- accept_count + 1
    } else {
      theta_chain[iter, ] <- theta_current
      theta_transformed_chain[iter, ] <- theta_trans_current
    }
    
    log_lik_chain[iter] <- log_lik_current
    
    # Update running sums for covariance computation
    theta_trans_accepted <- theta_transformed_chain[iter, ]
    sum_theta <- sum_theta + theta_trans_accepted
    sum_theta_sq <- sum_theta_sq + outer(theta_trans_accepted, theta_trans_accepted)
  }
  
  # Post-processing
  burnin <- max(t0, round(num_iterations * 0.2))  # Use at least t0 as burnin
  post_burnin <- theta_chain[(burnin+1):num_iterations, ]
  
  if (verbose) {
    cat("\n========================================\n")
    cat("MCMC complete!\n")
    cat("Final acceptance rate:", round(accept_count / num_iterations, 3), "\n")
    cat("Burnin period:", burnin, "iterations\n")
    cat("Posterior samples:", nrow(post_burnin), "\n")
  }
  
  return(list(
    theta = theta_chain,
    theta_transformed = theta_transformed_chain,
    log_likelihood = log_lik_chain,
    acceptance_rate = accept_count / num_iterations,
    proposal_cov=C_t,
    posterior_samples = post_burnin,
    posterior_mean = colMeans(post_burnin),
    posterior_sd = apply(post_burnin, 2, sd),
    posterior_quantiles = apply(post_burnin, 2, quantile, 
                                probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
    burnin = burnin,
    final_proposal_cov = C_t
  ))
}
