
#' Forward filtering backward Sampling for penalised Langevin SDE 
#' see Doucet et .al (2011) for the general algorithm.
#' Currently, works for Strang but only with naive splitting
#' @param data A data.frame containing observations with columns: time, Y1, Y2
#' @param n_samples number of smoothed trajectories to sample
#' @param forward_filter output of the particle_filter2D function. If NULL, the function
#'  will run the particle filter first.
#' @param sde_params list with elements tau, nu, omega
#' @param potential_params parameters for the potential gradient function
#' @param error_params parameters for the measurement error distribution
#' @param error_dist type of measurement error : "normal", "scaled_t" or "argos"
#' @param polygon matrix with two columns defining the vertices of the polygonal domain
#' @param U0 initial state (vector of length 4)
#' @param lambda penalty parameter for the inward push
#' @param num_particles number of particles to use in the particle filter
#' @param scheme type of splitting scheme to apply "Lie-Trotter" or "Strang"
#' @param split_around_fixed_point logical. If TRUE, split the SDE around the fixed point
#' @param verbose logical. If TRUE, print information during the computation
#' @return array of dimension (n_samples, num_states, N) with the smoothed trajectories7
#' @importFrom mvtnorm dmvnorm
#' @export
forward_filtering_backward_sampling <-
  function(data,n_samples,forward_filter=NULL,sde_params,
           potential_params=NULL,error_params,error_dist,polygon,U0,lambda,
           num_particles,scheme="Lie-Trotter",split_around_fixed_point=FALSE,
           verbose=FALSE) {
  
  # forward filtering
  if (is.null(forward_filter)) {
    forward_filter = particle_filter2D(data,sde_params,potential_params,
                                     error_params,error_dist,
                                     polygon,U0,lambda,num_particles,scheme,
                                     split_around_fixed_point,verbose)
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
          U_hat<-solve_ODE(U_prev,delta,push_array[k,,j],potential_params,
                           ind_fixed_point_current)
          
          #SDE mean and covariance
          OU_solution<-solve_SDE(U_hat,delta,tau,nu,omega,potential_params,
                                 ind_fixed_point_current)
          Q<-OU_solution$Q
          mean<-OU_solution$mean
          
          backward_weights[k,j] <- weights[k, j]*dmvnorm(U_next,
                                                         mean = mean,
                                                         sigma = Q)
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
                               ind_fixed_point_current)
          #SDE mean and covariance
          OU_solution<-solve_SDE(U_hat,delta,tau,nu,omega,potential_params,
                                 ind_fixed_point_current)
          Q<-OU_solution$Q
          mean<-OU_solution$mean
          
          X_next <- U_next[1:2]
          V_next <- U_next[3:4]
          
          push_next<-compute_push(X_next,polygon,lambda)
          potential_grad_next<-mix_gaussian_grad_cpp(X_next, x_star,
                                                     list(B=B,alpha=alpha), 
                                                     exclude=integer(0))
          
          V_tilde <- V_next + (delta/2) * (push_next + potential_grad_next)
          U_tilde_next <- c(X_next, V_tilde)
          
          
          backward_weights[k,j] <- weights[k, j]*dmvnorm(U_tilde_next,
                                                         mean = mean,
                                                         sigma = Q)
          
          
        }
      }
      
      bw <- backward_weights[,j]
      smoothed_trajectories[m, , j] <- 
        particles[sample(1:num_particles, 1, prob = bw/sum(bw)), , j]
    }
  }
  return (smoothed_trajectories)
}
    