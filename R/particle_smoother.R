
#' Forward filtering backward Sampling for penalised Langevin SDE 
#' see Doucet et .al (2011) for the general algorithm.
#' @param data A data.frame containing observations with columns: time, Y1, Y2
#' @param n_samples number of smoothed trajectories to sample
#' @param forward_filter output of the particle_filter2D_cpp function. If NULL, the function
#'  will run the particle filter first.
#' @param sde_params list with elements tau, nu, omega
#' @param potential_params parameters for the potential gradient function.
#' NULL (default) means no potential, encoded as a single component with zero
#' weight, identity precision matrix and centre (0, 0)
#' @param error_params parameters for the measurement error distribution
#' @param error_dist type of measurement error : "normal", "scaled_t" or "argos"
#' @param polygon \code{sp::Polygon} object defining the polygonal domain.
#' Can be NULL when lambda = Inf (no domain penalisation)
#' @param U0 initial state (vector of length 4)
#' @param lambda penalty parameter for the inward push
#' @param num_particles number of particles to use in the particle filter
#' @param scheme type of splitting scheme to apply "Lie-Trotter" or "Strang"
#' @param split_around_fixed_point logical. If TRUE, split the SDE around the fixed point.
#' @param ESS_threshold threshold for effective sample size to trigger
#' resampling in the particle filter
#' @param verbose logical. If TRUE, print information during the computation
#' @param print_timing logical. If TRUE, print timing information for each step
#' @param error_params_per_obs Optional list of length N of per-observation error
#' parameter lists (e.g. from \code{make_argos_obs_params}), passed to
#' \code{particle_filter2D_cpp}. Default NULL
#' @return array of dimension (n_samples, num_states, N)
#' with the smoothed trajectories
#' @importFrom mvtnorm dmvnorm
#' @export
forward_filtering_backward_sampling <-
  function(data,n_samples,forward_filter=NULL,sde_params,
           potential_params=NULL,error_params,error_dist,polygon,U0,lambda,
           num_particles,scheme="Lie-Trotter",split_around_fixed_point=FALSE,
           ESS_threshold=1,
           verbose = FALSE,print_timing = FALSE,
           error_params_per_obs = NULL) {
  if (is.null(polygon) && is.finite(lambda)) {
    stop("polygon must be provided when lambda is finite")
  }

  # No potential: single component with zero weight, so H and its gradient
  # are identically zero
  if (is.null(potential_params)) {
    potential_params <- list(x_star = matrix(0, 1, 2), B = list(diag(2)), alpha = 0)
  }

  # forward filtering

  if (is.null(forward_filter)) {

    polygon_coords <- if (is.null(polygon)) matrix(numeric(0), 0, 2) else polygon@coords

    forward_filter = particle_filter2D_cpp(observations=as.matrix(data[,c("time","Y1","Y2")]),
                                           sde_params = sde_params, potential_params = potential_params,
                                           error_params = error_params, error_dist = error_dist,
                                           polygon_coords = polygon_coords, U0 = U0, lambda = lambda,
                                           num_particles = num_particles,
                                           split_around_fixed_point = split_around_fixed_point,
                                           scheme = scheme,ESS_threshold=ESS_threshold,
                                           verbose = verbose,print_timing = print_timing,
                                           error_params_per_obs = error_params_per_obs)
  }
    
  particles =  forward_filter$particles
  weights = forward_filter$weights
  push_array = forward_filter$push # inward push at each time step
  ind_fixed_point = forward_filter$ind_fixed_point
  N <- dim(particles)[3] #number of observations
  num_states <- dim(particles)[2] #number of states
  
  # movement parameters
  tau<-sde_params$tau
  nu <- sde_params$nu;omega<-sde_params$omega
  
  smoothed_trajectories <- array(NA, dim = c(n_samples, num_states, N))
  
  backward_weights <- matrix(0, nrow = num_particles, ncol = N)
  
  for (m in 1:n_samples) {
    
    # Sample final state from the particle set at time N
    smoothed_trajectories[m, , N] <- 
      particles[sample(1:num_particles, 1, prob = weights[, N]), , N]
    
    for (j in (N-1):1) {
      
      for (k in 1:num_particles) {
        
        U_prev <- particles[k, , j]
        delta <- data$time[j + 1] - data$time[j]
        U_next<-smoothed_trajectories[m, , j + 1]
        
        if (scheme=="Lie-Trotter") {
          
          if (split_around_fixed_point) {
              ind_fixed_point_current <- ind_fixed_point[k,j]
            } else {
              ind_fixed_point_current <- NULL
            }
          #ODE step 
          U_hat<-solve_ODE_cpp(U_prev,delta,push_array[k,,j],potential_params,
                           nu,ind_fixed_point_current)
          
          #SDE mean and covariance
          OU_solution<-solve_SDE_cpp(U_hat,delta,tau,nu,omega,potential_params,
                                 ind_fixed_point_current,NULL,NULL)
          Q<-OU_solution$Q
          mean<-OU_solution$mean
          
          backward_weights[k,j] <- weights[k, j]*exp(log_dmvnorm_chol_cpp(U_next,
                                                         mean = mean,
                                                         cholSigma = chol_cpp(Q)))
        }
        
        if (scheme=="Strang") {
          
          if (split_around_fixed_point) {
            ind_fixed_point_current <- ind_fixed_point[k,j]
          } else {
            ind_fixed_point_current <- NULL
          }
          
          alpha <- potential_params$alpha; B <- potential_params$B
          x_star <- potential_params$x_star
          
          
          #ODE solution
          U_hat<-solve_ODE_cpp(U_prev,delta/2,push_array[k,,j],potential_params,
                               nu,ind_fixed_point_current)
          #SDE mean and covariance
          OU_solution<-solve_SDE_cpp(U_hat,delta,tau,nu,omega,potential_params,
                                 ind_fixed_point_current,NULL,NULL)
          Q<-OU_solution$Q
          mean<-OU_solution$mean
          
          X_next <- U_next[1:2]
          V_next <- U_next[3:4]
          
          push_next<-compute_push(X_next,polygon,lambda)
          # Push and gradient scaled by 2*nu^2/pi — see llk_gradient.R docs.
          nu_scale <- 2*nu^2/pi
          if (split_around_fixed_point) {
            l <- ind_fixed_point_current
            B_l <- B[[l]]; alpha_l <- alpha[l]; x_star_l <- x_star[l,]
            e_l <- as.numeric(exp(-t(X_next-x_star_l)%*%B_l%*%(X_next-x_star_l)))
            grad_next <- mix_gaussian_grad_cpp(X_next, x_star,
                                               list(B=B,alpha=alpha),
                                               exclude=l)+
              2*alpha_l*(e_l-1)*B_l%*%(X_next-x_star_l)
          } else {
            grad_next <- mix_gaussian_grad_cpp(X_next, x_star,
                                               list(B=B,alpha=alpha),
                                               exclude=integer(0))
          }

          V_tilde <- V_next + (delta/2) * nu_scale * (push_next + grad_next)
          U_tilde_next <- c(X_next, V_tilde)
          
          
          backward_weights[k,j] <- weights[k, j]*exp(log_dmvnorm_chol_cpp(U_tilde_next,
                                                         mean = mean,
                                                         cholSigma = chol_cpp(Q)))
        }
      }
      
      bw <- backward_weights[,j]
      smoothed_trajectories[m, , j] <- 
        particles[sample(1:num_particles, 1, prob = bw/sum(bw)), , j]
    }
  }
  return (smoothed_trajectories)
}
    