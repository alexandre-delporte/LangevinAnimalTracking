

#' Kalman filter algorithm for penalised Langevin SDE based on 
#' Lie-Trotter approximation
#' @param data data frame with columns Y1,Y2 for observations, time for time
#' @param sde_params List of parameters tau,nu,omega
#' @param potential_params Parameters for a mixture of Gaussian potential.
#' @param sigma_obs Standard deviation of the (isotropic) measurement error
#' @param polygon polygon defining the ecological domain of interest
#' @param lambda Penalisation parameter. Default : Inf (no penalisation)
#' @param U0 Initial guess for the hidden state. A vector of length 4 with position and velocity for RACVM, 
#' and a vector of length for OU
#' @param R0 Initial guess for the covariance matrix of the hidden state. Dimensions 
#' 4x4 for RACVM and 2x2 for OU
#' @param split_around_fixed_point Boolean. If TRUE, split around fixed point of 
#' the potential. Only works for mixture of Gaussian potential
#' 
#' @return matrix with filtered states (position and velocity). 
#' Dimension (\code{nrow(data)},4) 
#' @export
#' @importFrom expm expm
#'

kalman_filter <- function(data,sde_params,potential_params=NULL,
                          sigma_obs,polygon,lambda=Inf,U0, R0,split_around_fixed_point=FALSE) {
  
  # Check data
  required_cols <- c("Y1", "Y2", "time")
  if (!all(required_cols %in% colnames(data))) {
    stop("`data` must contain columns: Y1, Y2, and time.")
  }
  
   # Check sde parameters
  required_params <- c("tau", "nu", "omega")
  missing_params <- setdiff(required_params, names(sde_params))
  if (length(missing_params) > 0) {
    stop(paste("Missing parameters in `sde_params`: ", paste(missing_params, collapse = ", ")))
  }
  
  # Check initial state
  if (length(U0) != 4) {
    stop("`U0` must be a numeric vector of length 4 for RACVM.")
  }
  
  # Check initial covariance matrix
  expected_dim <- length(U0)
  if (!all(dim(R0) == c(expected_dim, expected_dim))) {
    stop(paste("`R0` must be a", expected_dim, "x", expected_dim, "matrix."))
  }
  
  # Check potential parameters if splitting around fixed point
  if (split_around_fixed_point && is.null(potential_params)) {
    stop("split_around_fixed_point = TRUE requires non-null potential_params")
  }
  
  
  tau<-sde_params$tau;nu <- sde_params$nu;omega<-sde_params$omega
  
  n_steps <- nrow(data)
  U_hat <- matrix(0, nrow = n_steps, ncol = length(U0)) 
  U_hat[1, ] <- U0
  R <- R0 
  Y<-as.matrix(data[,c("Y1","Y2")])
  H <- cbind(diag(2), matrix(0, 2, 2)) # Observation matrix
 
  alpha<-potential_params$alpha;B<-potential_params$B
  x_star<-potential_params$x_star
  Gamma <- rbind(
    cbind(matrix(0,2,2), matrix(0,2,2)),
    cbind(matrix(0,2,2), 4*nu^2/pi/tau*diag(2)))
    
  
  for (k in 2:n_steps) {
    
    dt<-data$time[k]-data$time[k-1]
    
    X_prev <- U_hat[k - 1, 1:2]
    
    push<-compute_push(X_prev,polygon,lambda)
    
    
    if (split_around_fixed_point) {
      
      #compute mahalanobis distances
      quad <- sapply(seq_along(alpha), function(idx) {
        t(X_prev - x_star[idx, ]) %*% B[[idx]] %*% (X_prev - x_star[idx, ])
      })
      
      #choose center
      l<-choose_center(X_prev, x_star, list(B=B,alpha=alpha))
      
      #extract parameters for the chosen center
      B_l<-B[[l]];alpha_l<-alpha[l];x_star_l<-x_star[l,]
      
      #compute non linear term
      e_l<-as.numeric(exp(-quad[l]))
      gv<-push+mix_gaussian_grad_cpp(X_prev,x_star,
                                 list(B=B,alpha=alpha),
                                 exclude=l)+
        2*alpha_l*(e_l-1)*B_l%*%(X_prev-x_star_l)
      
      
      Atilde<-rbind(cbind(matrix(0,ncol=2,nrow=2),diag(2)),
                    cbind(-2*alpha_l*B_l,-matrix(c(1/tau,-omega,omega,1/tau),
                                                 nrow=2,byrow=TRUE)))
      
      expAh<-expm(Atilde*dt)
      Qtilde<-OU_cov_exact(Atilde, Gamma,dt,expAh)
      
      U_pred <- expAh %*%(U_hat[k - 1, ]-c(x_star_l,0,0) -dt*c(0,0,gv))+c(x_star_l,0,0)
      R_pred <- expAh %*% R %*% t(expAh) + Qtilde
      
    }
    
    else {
      
      L <- RACVM_link(tau,omega,dt)
      Q<- RACVM_cov(tau,nu,omega, dt)
      
      # Prediction step
      fv <- push+mix_gaussian_grad_cpp(X_prev,x_star,list(B=B,alpha=alpha),
                                       exclude=integer(0))
      
      
      U_pred <- L %*%U_hat[k - 1, ] -dt*L%*% c(0,0,fv)
      R_pred <- L %*% R %*% t(L) + Q
    }
      
    # Correction step
    K <- R_pred %*% t(H) %*% solve(H %*% R_pred %*% t(H) + sigma_obs^2 * diag(2))
    U_hat[k, ] <- U_pred + K %*% (Y[k, ] - H %*% U_pred)
    R <- (diag(length(U0)) - K %*% H) %*% R_pred
  }
  
  return(U_hat)
}


#' Extended Kalman filter algorithm for penalised Langevin SDE 
#' based on Lie-Trotter approximation
#' @param data data frame with columns Y1,Y2 for observations, time for time
#' @param sde_params List of parameters tau,nu,omega 
#' @param potential_hessian Function to compute Hessian matrix of the potential surface
#' @param potential_params Parameters for a mixture of gaussian potential. Only necessary
#' if split_around_fixed_point is TRUE.
#' @param sigma_obs Standard deviation of the (isotropic) measurement error
#' @param polygon polygon defining the ecological domain of interest
#' @param lambda Penalisation parameter. Default : Inf (no penalisation)
#' @param U0 Initial guess for the hidden state. A vector of length 4 with position and velocity for RACVM, 
#' and a vector of length for OU
#' @param R0 Initial guess for the covariance matrix of the hidden state. Dimensions 
#' 4x4 for RACVM and 2x2 for OU
#' @param split_around_fixed_point Boolean. If TRUE, split around fixed point of
#' the potential. Only works for mixture of mixture of gaussian potential
#' 
#' @return matrix with filtered states. Dimension (\code{nrow(data)},4)
#' @export
#' @importFrom expm expm
extended_kalman_filter <- function(data,sde_params,potential_hessian,
                                   potential_params=NULL,sigma_obs,
                                   polygon, lambda=Inf, U0, R0,split_around_fixed_point=FALSE) {
  
  # Check data
  required_cols <- c("Y1", "Y2", "time")
  if (!all(required_cols %in% colnames(data))) {
    stop("`data` must contain columns: Y1, Y2, and time.")
  }
  
  # Check sde_params
  required_params <- c("tau", "nu", "omega")
  
  missing_params <- setdiff(required_params, names(sde_params))
  if (length(missing_params) > 0) {
    stop(paste("Missing parameters in `sde_params`: ", paste(missing_params, collapse = ", ")))
  }
  
  # Check initial state
  if (length(U0)!=4) {
    stop("`U0` must be a numeric vector of length 4 for RACVM.")
    }
  
  
  # Check initial covariance matrix
  expected_dim <- length(U0)
  if (!all(dim(R0) == c(expected_dim, expected_dim))) {
    stop(paste("`R0` must be a", expected_dim, "x", expected_dim, "matrix."))
  }
  
  n_steps <- nrow(data)
  U_hat <- matrix(0, nrow = n_steps, ncol = length(U0)) 
  U_hat[1, ] <- U0
  R <- R0 
  Y<-as.matrix(data[,c("Y1","Y2")])
  H <- cbind(diag(2), matrix(0, 2, 2)) # Observation matrix
  
  tau<-sde_params$tau;nu <- sde_params$nu;omega<-sde_params$omega
  
  if (split_around_fixed_point) {
    
    alpha<-potential_params$alpha;B<-potential_params$B
    x_star<-potential_params$x_star
    Gamma <- rbind(
      cbind(matrix(0,2,2), matrix(0,2,2)),
      cbind(matrix(0,2,2), 4*nu^2/pi/tau*diag(2)))
  }
  
  for (k in 2:n_steps) {
    
    dt<-data$time[k]-data$time[k-1]
    X_prev <- U_hat[k - 1, 1:2]
    
    if (is_point_inside_polygon(X_prev, polygon)) {
      projection<-list(point=X_prev,gradient=diag(2))
      push<-0
    }
    else {
      projection <-closest_point_on_boundary(X_prev, polygon)
      push <- (X_prev-projection$point)/lambda
    }
    
    if (is.null(projection$gradient)) {
      projection_grad <- diag(2)
    } else {
      projection_grad <- projection$gradient
    }
    
    if (split_around_fixed_point) {
      
      #compute mahalanobis distances
      quad<-sapply(1:length(alpha),function(k) {
        t(X_prev-x_star[k,])%*%B[[k]]%*%(X_prev-x_star[k,])})
      
      #choose center
      l<-choose_center(X_prev, x_star, list(B=B,alpha=alpha))
      
      #extract parameters for the chosen center
      B_l<-B[[l]];alpha_l<-alpha[l];x_star_l<-x_star[l,]
      
      #compute non linear term
      e_l<-as.numeric(exp(-quad[l]))
      gv<-push+mix_gaussian_grad(X_prev,x_star,
                                 list(B=B,alpha=alpha),
                                 exclude=l)+
        2*alpha_l*(e_l-1)*B_l%*%(X_prev-x_star_l)
      
      
      Atilde<-rbind(cbind(matrix(0,ncol=2,nrow=2),diag(2)),
                    cbind(-2*alpha_l*B_l,-matrix(c(1/tau,-omega,omega,1/tau),
                                                 nrow=2,byrow=TRUE)))
      
      expAh<-expm(Atilde*dt)
      Qtilde<-OU_cov_exact(Atilde, Gamma,dt,expAh)
      
      Dxgv<-(diag(2)-projection_grad)/lambda+ 
        mix_gaussian_hessian(X_prev, x_star,list(B=B,alpha=alpha),l)+
        2*alpha_l*((e_l-1)*B_l-2*e_l*B_l%*%(X_prev-x_star_l)%*%t(X_prev-x_star_l)%*%B_l)
      
      G_prev <- rbind(
        cbind(matrix(0, 2, 2), matrix(0, 2, 2)),
        cbind(Dxgv,matrix(0, 2, 2)))
      
      #prediction
      U_pred <- expAh %*%(U_hat[k - 1, ]-c(x_star_l,0,0) -dt*c(0,0,gv))+c(x_star_l,0,0)
      R_pred <- expAh %*% R %*% t(expAh)+dt^2* expAh%*%G_prev%*%R%*%t(G_prev)%*%t(expAh)+
        Qtilde
    }
    
    else {
      
      L<-RACVM_link(tau,omega,dt)
      Q<-RACVM_cov(tau,nu,omega,dt)
      
      f_prev <- c(rep(0, 2),push+mix_gaussian_grad_cpp(X,x_star,list(B=B,alpha=alpha),
                                                       exclude=integer(0)))
      
      F_prev <- rbind(
        cbind(matrix(0, 2, 2), matrix(0, 2, 2)),
        cbind((diag(2)-projection_grad)/lambda+potential_hessian(X_prev),matrix(0, 2, 2)))
    
      #prediction
      U_pred <- L %*% U_hat[k - 1, ] -dt*L%*% f_prev
      R_pred <- L %*% R %*% t(L) +dt^2*L%*% F_prev %*% R %*% t(F_prev) %*% t(L) + Q
    }
    
      # correction
    K <- R_pred %*% t(H) %*% solve(H %*% R_pred %*% t(H) + sigma_obs^2 * diag(2))
    U_hat[k, ] <- U_pred + K %*% (Y[k, ] - H %*% U_pred)
    R <- (diag(4) - K %*% H) %*% R_pred
  }
  
  return(U_hat)
}

