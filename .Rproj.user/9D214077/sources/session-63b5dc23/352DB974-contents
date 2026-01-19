



#' Joint log-likelihood for one step of the SDE
#' @param Y numeric vector of observation (length 2)
#' @param U_next numeric vector c(X1,X2,V1,V2) of next hidden state
#' @param U_prev numeric vector c(X1,X2,V1,V2) of previous hidden state
#' @param delta time step between U_prev and U_next
#' @param push_prev vector of length 2 with the value of the inward push at time 
#' of U_prev
#' @param push_next vector of length 2 with the value of the inward push at time
#' of U_next. Only necessary for Strang scheme
#' @param potential_params Parameters for a mixture of Gaussian potential.
#' @param tau persistence parameter tau of the SDE
#' @param nu speed parameter nu of the SDE
#' @param omega rotation parameter omega of the SDE
#' @param error_dist Type of measurement error : "normal", "scaled_t" or "argos"
#' @param error_params List of parameters for the measurement error distribution
#' @param scheme Type of splitting scheme to apply "Lie-Trotter" or "Strang"
#' @param ind_fixed_point Default NULL. 
#' If not NULL, index of the fixed point around which to split
#' @param verbose logical. If TRUE, print information during the computation

#' @return log-likelihood value (numeric) 
#' @importFrom mvtnorm dmvnorm
#'  
llk_one_step<-function(Y,U_next,U_prev,delta,push_prev,push_next=NULL,
                       potential_params=NULL,tau,nu,omega,
                       error_dist=NULL,error_params=NULL,scheme="Lie-Trotter",
                       ind_fixed_point=NULL,
                       verbose=FALSE) {

  
  if (scheme=="Lie-Trotter") {
    
      U_hat<-solve_ODE(U_prev,delta,push_prev,
                       potential_params,
                       ind_fixed_point)
      
      #SDE mean and covariance
      OU_solution<-solve_SDE(U_hat,delta,tau,nu,omega,potential_params,
                             ind_fixed_point)
      Q<-OU_solution$Q
      mean<-OU_solution$mean
      
      latent_llk = dmvnorm(U_next,mean=mean,sigma=Q,log=TRUE)
      obs_llk = switch(error_dist,
                       normal   = dmvnorm(Y,mean=U_next, sigma=error_params$sigma_obs^2*diag(2),
                                          log=TRUE),
                       scaled_t = dscaledt(Y[1],mean=U_next[1], error_params$scale, error_params$df,log=TRUE)
                       +dscaledt(Y[2], mean=U_next[2],error_params$scale, error_params$df,log=TRUE),
                       argos  = dmvt_mixture(Y,U_next,error_params,log=TRUE),
                       stop("Unknown error_dist"))
  }
  if (scheme=="Strang") {
    
    U_hat<-solve_ODE(U_prev,delta/2,push_prev,
                     potential_params,
                     ind_fixed_point)
    
    #SDE mean and covariance
    OU_solution<-solve_SDE(U_hat,delta,tau,nu,omega,potential_params,
                           ind_fixed_point)
    Q<-OU_solution$Q
    mean<-OU_solution$mean
    
    latent_llk = dmvnorm(U_next+c(0,0,delta/2*(push_next+potential_grad(U_next[1:2]))),
                         mean=mean,sigma=Q,log=TRUE)
    obs_llk = switch(error_dist,
                     normal = dmvnorm(Y,mean=U_next, sigma=error_params$sigma_obs^2*diag(2),
                                        log=TRUE),
                     scaled_t = dscaledt(Y[1],mean=U_next[1], error_params$scale, error_params$df,log=TRUE)
                     +dscaledt(Y[2], mean=U_next[2],error_params$scale, error_params$df,log=TRUE),
                     argos    = dmvt_mixture(Y,U_next,error_params,log=TRUE),
                     stop("Unknown error_dist"))
  }
  
  return (latent_llk+obs_llk)
}


#' Joint log-likelihood for the full data set
#' @param data data frame with columns time, Y1, Y2, X1, X2, V1, V2
#' @param push_mat matrix with two columns with the inward push values at each time step
#' @param potential_params Parameters for a mixture of Gaussian potential. Only necessary
#' if split_around_fixed_point is TRUE.
#' @param tau persistence parameter tau of the SDE
#' @param nu speed parameter nu of the SDE
#' @param omega rotation parameter omega of the SDE
#' @param error_dist Type of measurement error : "normal", "scaled_t" or "argos"
#' @param error_params List of parameters for the measurement error distribution
#' @param scheme Type of splitting scheme to apply "Lie-Trotter" or "Strang"
#' @param split_around_fixed_point logical. If TRUE, split the SDE around the fixed point
#' @param verbose logical. If TRUE, print information during the computation
#' 
#' @return vector of log-likelihood values (numeric) for each time step
#' @export
#' @importFrom mvtnorm dmvnorm
llk<-function(data,push_mat,
             potential_params=NULL,tau,nu,omega,
             error_dist=NULL,error_params=NULL,scheme="Lie-Trotter",
             split_around_fixed_point=TRUE,
             verbose=FALSE) {
  
  n_steps <- nrow(data)
  total_llk <- numeric(n_steps-1)
  
  for (j in 1:(n_steps-1)) {
    
    Y = U_prev <- as.numeric(data[j,c("Y1","Y2")])
    U_prev <- as.numeric(data[j,c("X1","X2","V1","V2")])
    U_next <- as.numeric(data[j+1,c("X1","X2","V1","V2")])
    delta <- as.numeric(data[j+1,"time"] - data[j,"time"])
    push_prev <- push_mat[j,]
    push_next <- push_mat[j+1,]
    if (verbose) {
      cat("---- Step", j, "----\n")
      cat("delta:", delta, "\n")
      cat("U_prev:", U_prev, "\n")
      cat("U_next:", U_next, "\n")
      cat("push:", push, "\n")
    }
    
    ind_fixed_point<-if (split_around_fixed_point)
      choose_center(U_prev[1:2], potential_params$x_star, 
                    list(B=potential_params$B,
                         alpha=potential_params$alpha)) 
    else NULL
    step_llk <- llk_one_step(Y,U_next,U_prev,delta,push_prev,push_next,
                                       potential_params,
                                       tau,nu,omega,
                                       error_dist,error_params,scheme,
                                       ind_fixed_point,
                                       verbose)
    total_llk[j] <- step_llk
  }
  
  
  return(total_llk)
}