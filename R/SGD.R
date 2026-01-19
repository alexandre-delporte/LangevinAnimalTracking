

#' Stochastic Gradient Descent with Fisher Information Preconditioning for SDE parameter estimation 
#' (Baey et .al 2023). Only works with Naive (split around 0) Lie-Trotter scheme for now.
#' @param data Data frame with observed data
#' @param sde_params Named list of initial SDE parameters
#' @param fixpar Character vector of parameter names to fix during optimization
#' @param SGD_iter Number of SGD iterations
#' @param potential_params Parameters ofr the mixture of gaussian potential. Only
#' needed when split_around_fixed_point=TRUE
#' @param error_dist Error distribution for the observation model
#' @param error_params Parameters for this error distribution
#' @param scheme Splitting scheme ("Lie-Trotter" or "Strang")
#' @param polygon Polygon defining the domain of interest
#' @param U0 Initial position and velocity vector
#' @param lambda Penlization parameter
#' @param num_particles Number of particles for the particle filter
#' @param split_around_fixed_point Logical indicating whether to split the potential around its fixed point
#' @param verbose Logical indicating whether to print progress messages
#' @param gamma0 Initial step size parameter
#' @param K_preheat Number of preheating iterations
#' @param alpha Step size decay parameter
#' @param smoothing Logical indicating whether to use smoothing when inputting a latent trajectory
#' @return A list containing:
#' #' - theta: Matrix of parameter estimates at each iteration
#' #' - Delta: Matrix of running gradient estimates
#' @export

SGD_Fisher <- function(data, sde_params, fixpar = NULL, SGD_iter,
                       potential_params = NULL, error_dist = NULL, error_params = NULL,
                       scheme = "Lie-Trotter", polygon, U0, lambda,
                       num_particles, split_around_fixed_point = FALSE,
                       verbose = FALSE, gamma0 = 1e-4, K_preheat = 1000,
                       alpha = 2/3,smoothing=FALSE) {
  
  if (scheme %in% c( "Strang")) {
    stop("'Strang' implemented yet")
  }
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
 
  Delta <- matrix(0, ,ncol=p,nrow=n-1)
  I_p   <- diag(1, p)
  
  alpha<-potential_params$alpha;B<-potential_params$B
  x_star<-potential_params$x_star
  
  #pre compute true potential and true inward push to later compare
  # true and estimate gradients
  push_mat_true <- matrix(NA, nrow = n, ncol = 2)
  potential_grad_mat_true <- matrix(NA, nrow = n, ncol = 2) 
  for (t in 1:n) { 
    X <- as.numeric(data[t, c("X1", "X2")])
    push_mat_true[t, ] <- compute_push(X, polygon, lambda)
    potential_grad_mat_true[t, ] <- mix_gaussian_grad_cpp(X,x_star,list(B=B,alpha=alpha),
                                                          exclude=integer(0))
  } 
  
  
  for (k in 1:SGD_iter) {
    if (verbose) message("SGD iter", k, "\n")

    # Run particle filter
    filter <- particle_filter2D(data,sde_params = as.list(theta[k, ]),
                                potential_params = potential_params,
                                error_params = error_params, error_dist = error_dist,
                                polygon = polygon, U0 = U0, lambda = lambda,
                                num_particles = num_particles, scheme = scheme,
                                split_around_fixed_point = split_around_fixed_point,
                                verbose = FALSE)
    
    if (smoothing) {
      
      backward_samples <- forward_filtering_backward_sampling(data,1,forward_filter=filter,
                                                as.list(theta[k, ]),potential_params=NULL,
                                               error_params,error_dist,polygon,
                                               U0,lambda,num_particles,
                                               scheme="Lie-Trotter",
                                               split_around_fixed_point=FALSE,
                                               verbose=FALSE)
      z<-t(apply(backward_samples, c(2, 3), mean))
    }
    else {
      # Get mean trajectory from particle filter
      z<-get_PF_mean(filter)
    }
 
    colnames(z) <- c("X1", "X2", "V1", "V2") 
      
    compare_trajectories_in_polygon(data, as.data.frame(z), polygon)
    
    n_obs = nrow(data)
    times = data$time[2:n_obs]
    plot(data$time,z[,"V1"],type='l',main=paste("Smoothed VS True V1 at iter",k))
    lines(data$time,data[,"V1"],col="blue")
    lines(times,diff(data$Y1)/diff(data$time),type='l',col="red",
          lwd=0.1)
    plot(data$time,z[,"V2"],type='l',main=paste("Smoothed VS True V2 at iter",k))
    lines(data$time,data[,"V2"],col="blue")
    lines(times,diff(data$Y2)/diff(data$time),type='l',col="red",
          lwd=0.1)
    # Compute push and potential gradient for this trajectory 
    push_mat <- matrix(NA, nrow = n, ncol = 2)
    potential_grad_mat <- matrix(NA, nrow = n, ncol = 2) 
    
    for (t in 1:n) { 
      X <- z[t, c("X1", "X2")]
      push_mat[t, ] <- compute_push(X, polygon, lambda)
      potential_grad_mat[t, ] <- potential_grad(X)
      } 
    
    # Compute gradient for this trajectory 
    z_df <- as.data.frame(z) 
    z_df$time <- data$time 
    
    grad_mat <- llk_gradient(data = z_df, push_mat = push_mat, 
                        potential_grad_mat = potential_grad_mat,
                        tau = theta[k, "tau"], nu = theta[k, "nu"],
                        omega = theta[k, "omega"], error_dist = NULL,
                        error_params = NULL, scheme = scheme, 
                        ind_fixed_point = NULL, verbose = FALSE )
    
    v_k = colMeans(grad_mat)
    names(v_k) <- param_names
    if (length(fix_idx) > 0) v_k[fix_idx] <- 0
    
    #compute true gradient at current parameters for comparison
    true_grad_mat = llk_gradient(data, push_mat = push_mat_true, 
                                            potential_grad_mat = potential_grad_mat_true,
                                            tau = theta[k, "tau"], nu = theta[k, "nu"],
                                            omega = theta[k, "omega"], error_dist = NULL,
                                            error_params = NULL, scheme = scheme, 
                                            ind_fixed_point = NULL, verbose = FALSE )

    # compute gamma_k
    if (k <= K_preheat) {
      gamma_k <- gamma0 + (1 - gamma0) * (k / K_preheat)
    } else {
      gamma_k <- k^(-alpha)
    }
    if (verbose) message("gamma_k:", gamma_k, "\n")
    
    # stochastic approximation updates
    Delta <- (1 - gamma_k) * Delta + gamma_k * grad_mat
    I_k=matrix(0,ncol=p,nrow=p)
    for (i in 1:nrow(Delta)) {
      I_k= I_k+Delta[i,]%*% t(Delta[i,])
    }
    
    I_k <- I_k/n
    
    # Preconditioner with pre-heating phase
    if (k < K_preheat) {
      r_k <- max(1, sum(diag(I_k)))
      P_k <- (1 - gamma_k) * r_k * I_p + gamma_k * I_k
    } else {
      P_k <- I_k
    }
    
    # Compute step direction and update parameters
    direction <- tryCatch(
      solve(P_k, v_k),
      error = function(e) {
        if (verbose) cat("solve failed; regularizing matrix \n")
        I_k <- 0.5 * (P_k + t(P_k)) + 1e-4 * diag(p)
        return(solve(P_k, v_k))
      }
    )
    direction <- as.numeric(direction)
    if (length(fix_idx) > 0) direction[fix_idx] <- 0
    theta[k + 1, ] <- theta[k, ] + gamma_k * direction
    
    if (verbose) {
      message(cat("v_k:", round(v_k, 4), "\n"))
      message(cat("v_k_true:", round(colMeans(true_grad_mat), 4),"\n"))
      message("I_k:", round(I_k, 4), "\n")
      message("eigenvalues I_k:", round(eigen(I_k)$values, 6), "\n")
      message("direction:", round(direction, 6), "\n")
      message("theta:", round(theta[k + 1, ], 6), "\n")
      
    }
  } 
  
  return(list(
    theta = theta,
    Delta = Delta))
}
