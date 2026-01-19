
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
#'  @export
#'
#' @importFrom MASS mvrnorm
#'
PMMH_langevin <- function(data, theta_init, num_iterations, num_particles,
                          potential_params = NULL,
                          error_params, error_dist,
                          polygon, U0, lambda,
                          proposal_cov = NULL,
                          prior_log_density = NULL,
                          adapt_proposal = TRUE,
                          burnin = 2000,fixpar=NULL,
                          verbose = TRUE) {
  
  p <- length(theta_init)
  param_names <- names(theta_init)
  
  freepar <- setdiff(param_names, fixpar)
  p_free <- length(freepar)
  
  # Store markov chain
  theta_chain <- matrix(0, nrow = num_iterations, ncol = p)
  colnames(theta_chain) <- param_names
  theta_chain[1, ] <- unlist(theta_init)
  
  log_lik_chain <- numeric(num_iterations)
  accept_count <- 0
  
  # Default proposal covariance (random walk on log scale)
  if (is.null(proposal_cov)) {
    proposal_sd <- c(tau = 0.2, nu = 0.05, omega = 0.1)
    proposal_sd <- proposal_sd[freepar]
    proposal_cov <- diag(proposal_sd^2)
  }
  
  # Prior uniform on log scale
  if (is.null(prior_log_density)) {
    prior_log_density <- function(theta) {
      # Uniform prior for scale parameters: p(theta) propto 1/theta
      if (any(theta <= 0)) return(-Inf)
      return(-sum(log(theta)))
    }
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
  log_prior_current <- prior_log_density(theta_init[freepar])
  
  
  if (verbose) {
    cat("Initial log-likelihood:", log_lik_current, "\n")
    cat("Initial log-prior:", log_prior_current, "\n")
  }
  
  # MCMC loop
  for (iter in 2:num_iterations) {
    if (verbose && iter %% 10 == 0) {
      cat("\n========================================\n")
      cat("PMMH iteration", iter, "/", num_iterations, "\n")
      cat("  Current theta:", round(theta_chain[iter-1, ], 4), "\n")
      cat("  Current log-lik:", round(log_lik_current, 2), "\n")
      cat("  Acceptance rate:", round(accept_count / (iter-1), 3), "\n")
    }
    
    theta_current <- theta_chain[iter - 1, ]
    
    theta_free_current <- theta_current[freepar]
    theta_log_free_current <- log(theta_free_current)
    
    
    # Propose new parameters 
    
    theta_log_free_prop <- theta_log_free_current +
      mvrnorm(1, mu = rep(0, p_free), Sigma = proposal_cov)
    
    theta_free_prop <- exp(theta_log_free_prop)
    
    theta_prop <- theta_current
    theta_prop[freepar] <- theta_free_prop
    names(theta_prop) <- param_names
    
    
    if (verbose && iter %% 10 == 0) {
      cat("  Proposed theta:", round(theta_prop, 4), "\n")
    }
    
    # Check for invalid proposals
    if (any(theta_prop <= 0) || any(is.na(theta_prop)) || any(is.infinite(theta_prop))) {
      if (verbose && iter %% 10 == 0) {
        cat("  REJECTED (invalid proposal)\n")
      }
      theta_chain[iter, ] <- theta_current
      log_lik_chain[iter] <- log_lik_current
      next
    }
    
    # run particle filter with proposed parameters
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
      # reject if error in particlefilter 
      theta_chain[iter, ] <- theta_current
      log_lik_chain[iter] <- log_lik_current
      next
    }
    
    log_lik_prop <- pf_prop$loglik
    
    if (verbose && iter %% 10 == 0) {
      cat("  Proposed log-lik:", round(log_lik_prop, 2), "\n")
    }
    
    
    log_prior_prop <- prior_log_density(theta_prop[freepar])
    
    log_jacobian <- sum(theta_log_free_prop - theta_log_free_current)
    
    
    # acceptance probability
    log_alpha <- (log_lik_prop + log_prior_prop) - 
      (log_lik_current + log_prior_current) +
      log_jacobian
    
    if (verbose && iter %% 10 == 0) {
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
      if (verbose && iter %% 10 == 0) {
        cat("  u =", round(u, 4), ", log(u) =", round(log(u), 4), "\n")
        if (accept) {
          cat("  *** ACCEPTED ***\n")
        } else {
          cat("  REJECTED\n")
        }
      }
    }
    
    if (accept) {
      theta_chain[iter, ] <- theta_prop
      log_lik_current <- log_lik_prop
      log_prior_current <- log_prior_prop
      accept_count <- accept_count + 1
    } else {
      theta_chain[iter, ] <- theta_current
    }
    
    log_lik_chain[iter] <- log_lik_current
    
    # Adapt proposal covariance during burn-in
    if (adapt_proposal && iter <= burnin && iter %% 10 == 0) {
      accept_rate <- accept_count / iter
      
      # Target acceptance rate between 0.15 and 0.4 
      if (accept_rate < 0.15) {
        proposal_cov <- proposal_cov * 0.8
        if (verbose) cat("\n  >>> Decreasing proposal variance (accept rate =", 
                         round(accept_rate, 3), ")\n")
      } else if (accept_rate > 0.4) {
        proposal_cov <- proposal_cov * 1.2
        if (verbose) cat("\n  >>> Increasing proposal variance (accept rate =", 
                         round(accept_rate, 3), ")\n")
      }
      
      
    }
  }
  
  # Post-processing
  post_burnin <- theta_chain[(burnin+1):num_iterations, ]

  
  return(list(
    theta = theta_chain,
    log_likelihood = log_lik_chain,
    acceptance_rate = accept_count / num_iterations,
    proposal_cov = proposal_cov,
    posterior_samples = post_burnin,
    posterior_mean = colMeans(post_burnin),
    posterior_sd = apply(post_burnin, 2, sd),
    posterior_quantiles = apply(post_burnin, 2, quantile, 
                                probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
    burnin = burnin
  ))
}