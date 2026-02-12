
#' Function to compute partial derivatives of RACVM covariance matrix Q (Delporte et .al 2025)
#' @param tau RACVM parameter tau
#' @param nu RACVM parameter nu
#' @param omega RACVM parameter omega
#' @param h time step
#' @return List with three elements: derivatives of Q with respect to
#' #' tau, nu, and omega : three 4x4 matrices: dQ_dtau, dQ_dnu, dQ_domega
#' @export
dRACVM_cov <- function(tau, nu, omega,h) {
  

  C <- 1/tau^2 + omega^2
  sigma <- 2*nu / sqrt(pi * tau)
  
  exp_h_tau <- exp(-h/tau)
  exp_2h_tau <- exp(-2*h/tau)
  
  cos_omegah <- cos(omega * h)
  sin_omegah <- sin(omega * h)
  
  # Compute q1, q2, gamma1, gamma2
  phi <- omega * sin_omegah - (1/tau) * cos_omegah
  
  f1 <- h - 2 * (phi / C) * exp_h_tau + 
    (tau/2) * ((omega^2 - 3/tau^2) / C - exp_2h_tau)
  
  q1 <- (sigma^2 / C) * f1
  
  q2 <- (2*nu^2/pi) * (1 - exp_2h_tau)
  
  psi1 <- 1 + exp_2h_tau - 2*exp_h_tau*cos_omegah
  gamma1 <- (sigma^2 / (2*C)) * psi1
  
  psi2 <- exp_h_tau * sin_omegah - (omega*tau/2) * (1 - exp_2h_tau)
  gamma2 <- (sigma^2 / C) * psi2
  
  #  Derivative with respect to tau 
  
  # dC/dtau
  dC_dtau <- -2/tau^3
  
  # dsigma/dtau
  dsigma_dtau <- -sigma / (2*tau)
  
  # d(sigma^2/C)/dtau
  d_sigma2_C_dtau <- (sigma^2/C) * (-1/tau + 2/(tau^3 * C))
  
  # Derivative of phi
  dphi_dtau <- cos_omegah / tau^2
  
  # Derivative of g1 = 2 * (phi/C) * exp(-h/tau)
  dg1_dtau <- 2 * exp_h_tau * (
    (C * cos_omegah/tau^2 + phi * 2/tau^3) / C^2 +
      (phi * h/tau^2) / C)
  
  # Derivative of k = (omega^2 - 3/tau^2)/C - exp(-2h/tau)
  dk_dtau <- (6*C + 2*omega^2 - 6/tau^2) / (tau^3 * C^2) - 
    (2*h/tau^2) * exp_2h_tau
  
  # Derivative of g2 = (tau/2) * k
  dg2_dtau <- (1/2) * ((omega^2 - 3/tau^2) / C - exp_2h_tau) + 
    (tau/2) * dk_dtau
  
  # Derivative of f1
  df1_dtau <- -dg1_dtau + dg2_dtau
  
  # dq1/dtau
  dq1_dtau <- d_sigma2_C_dtau * f1 + (sigma^2/C) * df1_dtau
  
  # dq2/dtau
  dq2_dtau <- -(4*nu^2*h) / (pi * tau^2) * exp_2h_tau
  
  # dpsi1/dtau
  dpsi1_dtau <- (2*h/tau^2) * (exp_2h_tau - cos_omegah * exp_h_tau)
  
  # dgamma1/dtau
  dgamma1_dtau <- (sigma^2/(2*C)) * (
    (-1/tau + 2/(tau^3 * C)) * psi1 + dpsi1_dtau
  )
  
  # dpsi2/dtau
  dpsi2_dtau <- (h/tau^2) * sin_omegah * exp_h_tau - 
    (omega/2) * (1 - exp_2h_tau) + 
    (omega*h/tau) * exp_2h_tau
  
  # dgamma2/dtau
  dgamma2_dtau <- (sigma^2/C) * (
    (-1/tau + 2/(tau^3 * C)) * psi2 + dpsi2_dtau
  )
  
  dQ_dtau <- matrix(c(
    dq1_dtau, 0, dgamma1_dtau, dgamma2_dtau,
    0, dq1_dtau, -dgamma2_dtau, dgamma1_dtau,
    dgamma1_dtau, -dgamma2_dtau, dq2_dtau, 0,
    dgamma2_dtau, dgamma1_dtau, 0, dq2_dtau
  ), nrow=4, ncol=4, byrow=TRUE)
  
  
  # Derivative with respect to nu 
  
  #  dQ/dnu = (2/nu) * Q
  dQ_dnu <- (2/nu) * matrix(c(
    q1, 0, gamma1, gamma2,
    0, q1, -gamma2, gamma1,
    gamma1, -gamma2, q2, 0,
    gamma2, gamma1, 0, q2
  ), nrow=4, ncol=4, byrow=TRUE)
  
  # Derivative with respect to omega 
  
  # dC/domega
  dC_domega <- 2*omega
  
  # d(sigma^2/C)/domega
  d_sigma2_C_domega <- -2*omega*sigma^2 / C^2
  
  # Derivative of phi
  dphi_domega <- sin_omegah * (1 + h/tau) + omega*h*cos_omegah
  
  # Derivative of g1
  dg1_domega <- (2*exp_h_tau/C^2) * (C*dphi_domega - 2*omega*phi)
  
  # Derivative of g2
  dg2_domega <- 4*omega / (tau * C^2)
  
  # Derivative of f1
  df1_domega <- -dg1_domega + dg2_domega
  
  # dq1/domega
  dq1_domega <- d_sigma2_C_domega * f1 + (sigma^2/C) * df1_domega
  
  # dq2/domega
  dq2_domega <- 0
  
  # dpsi1/domega
  dpsi1_domega <- 2*h*exp_h_tau * sin_omegah
  
  # dgamma1/domega
  dgamma1_domega <- (sigma^2/C^2) * (
    h*C*exp_h_tau*sin_omegah - omega*psi1
  )
  
  # dpsi2/domega
  dpsi2_domega <- h*exp_h_tau*cos_omegah - (tau/2)*(1 - exp_2h_tau)
  
  # dgamma2/domega
  dgamma2_domega <- (sigma^2/C^2) * (
    h*C*exp_h_tau*cos_omegah - (C*tau/2)*(1 - exp_2h_tau) - 2*omega*psi2
  )
  
  dQ_domega <- matrix(c(
    dq1_domega, 0, dgamma1_domega, dgamma2_domega,
    0, dq1_domega, -dgamma2_domega, dgamma1_domega,
    dgamma1_domega, -dgamma2_domega, dq2_domega, 0,
    dgamma2_domega, dgamma1_domega, 0, dq2_domega
  ), nrow=4, ncol=4, byrow=TRUE)
  
  
  return(list(
    tau = dQ_dtau,
    nu = dQ_dnu,
    omega = dQ_domega
  ))
}



#' Function to compute partial derivatives of RACVM link matrix T (Delporte et .al 2025)
#' @param tau RACVM parameter tau
#' @param omega RACVM parameter omega
#' @param h time step
#' @return List with three elements: derivatives of Q with respect to
#' #' tau, nu, and omega : three 4x4 matrices: dT_dtau, dT_dnu = 0, dQ_domega
#' @export
dRACVM_link=function(tau,omega,h) {
  
  C <- 1/tau^2 + omega^2
  
  R <- matrix(c(cos(omega*h),sin(omega*h),-sin(omega*h),cos(omega*h)),
           byrow=TRUE,nrow=2)
  
  M <- matrix(c(1/tau,omega,-omega,1/tau),nrow=2,byrow=TRUE)
  
  # dT/dtau
  
 dT_dtau <- matrix(0,ncol=4,nrow=4)
 
 dT_dtau[1:2,3:4] <- 2/tau^3/C^2*M%*%(diag(2)-exp(-h/tau)*R)-
   1/C/tau^2*(diag(2)-exp(-h/tau)*R)-h/tau^2*exp(-h/tau)/C*M%*%R
 
 dT_dtau[3:4,3:4] <- h/tau^2*exp(-h/tau)*R
  
 # dT/dnu
 dT_dnu <- matrix(0,ncol=4,nrow=4)
 
 # dT_domega 
 dT_domega <- matrix(0,ncol=4,nrow=4)
 dT_domega[3:4,3:4] <- h* exp(-h/tau)*matrix(c(-sin(omega*h),cos(omega*h),
                                              -cos(omega*h),-sin(omega*h)),byrow=TRUE,nrow=2)
 
 dT_domega[1:2,3:4] <- -1/C^2*2*omega*M%*%(diag(2)-exp(-h/tau)*R)+
   +1/C* matrix(c(0,1,-1,0),ncol=2,byrow=TRUE)%*%(diag(2)-exp(-h/tau)*R)-
   h/C*exp(-h/tau)*M%*%matrix(c(-sin(omega*h),cos(omega*h),
                                -cos(omega*h),-sin(omega*h)),byrow=TRUE,nrow=2)
 
 return(list(
   tau = dT_dtau,
   nu = matrix(0,ncol=4,nrow=4),
   omega = dT_domega
 ))
 
}


#' Compute gradient of log-likelihood for one time step
#' Supports both Lie-Trotter and Strang schemes
#' 
#' @param U_next State vector at next time step (X1, X2, V1, V2)
#' @param U_prev State vector at previous time step (X1, X2, V1, V2)
#' @param delta Time step size
#' @param push Push vector at previous time step (length 2)
#' @param potential_grad Gradient of potential at previous time step (length 2)
#' @param tau RACVM parameter tau
#' @param nu RACVM parameter nu
#' @param omega RACVM parameter omega
#' @param scheme Integration scheme: "Lie-Trotter" or "Strang"
#' @param push_next Push vector at next time step (length 2) - REQUIRED for Strang scheme
#' @param potential_grad_next Gradient of potential at next time step (length 2) 
#' - REQUIRED for Strang scheme
#' @param verbose If TRUE, print detailed computation steps
#' @return Numeric vector of length 3 with gradients w.r.t. tau, nu, omega
#'    
llk_gradient_one_step <- function(U_next, U_prev, delta, push, potential_grad,
                                  tau, nu, omega, scheme = "Lie-Trotter", 
                                  push_next = NULL, potential_grad_next = NULL,
                                  verbose = FALSE) {
  
  # Get covariance and link matrices and their derivatives
  Q_mat <- RACVM_cov(tau, nu, omega, delta)
  T_mat <- RACVM_link(tau, omega, delta)
  Q_derivs <- dRACVM_cov(tau, nu, omega, delta)
  T_derivs <- dRACVM_link(tau, omega, delta)
  
  # Compute inverse of Q
  Q_inv <- solve(Q_mat)
  
  # invQ derivatives
  Q_inv_derivs <- lapply(Q_derivs, function(dQ) {
    Q_inv %*% dQ %*% Q_inv 
  })
  
  #determinant derivative
  log_det_Q_derivs <- sapply(Q_derivs, function(dQ) {
    sum(diag(Q_inv %*% dQ))
  })
  
  if (scheme == "Strang") {
    
    # Check that required arguments are provided
    if (is.null(push_next) || is.null(potential_grad_next)) {
      stop("For Strang scheme, push_next and potential_grad_next must be provided")
    }
    
    X_next <- U_next[1:2]
    V_next <- U_next[3:4]
    
    V_tilde <- V_next + (delta/2) * (push_next + potential_grad_next)
    U_tilde_next <- c(X_next, V_tilde)
    
    g_prev <- c(0, 0, push + potential_grad)
    
    U_hat_half <- U_prev - (delta/2) * g_prev
    
    mu <- as.vector(T_mat %*% U_hat_half)
    
    mu_derivs <- lapply(T_derivs, function(dT) {
      as.vector(dT %*% U_hat_half)
    })
    
    r <- U_tilde_next - mu
    
  } else if (scheme == "Lie-Trotter") {
    
    g <- c(0, 0, push + potential_grad)
    
    mu <- as.vector(T_mat %*% (U_prev - delta * g))
    
    mu_derivs <- lapply(T_derivs, function(dT) {
      as.vector(dT %*% (U_prev - delta * g))
    })
    
    r <- U_next - mu
    
  } else {
    stop("scheme must be either 'Lie-Trotter' or 'Strang'")
  }
  
  grad <- numeric(3)
  names(grad) <- c("tau", "nu", "omega")
  
  for (i in 1:3) {
    term1 <- -0.5 * log_det_Q_derivs[i]
    
    term2 <- 0.5 * t(r) %*% Q_inv_derivs[[i]] %*% r
    
    term3 <- t(r) %*% Q_inv %*% mu_derivs[[i]]
    
    grad[i] <- term1 + term2 + term3
    
    if (verbose) {
      cat(sprintf("Parameter %s: term1=%.6f, term2=%.6f, term3=%.6f, total=%.6f\n",
                  names(grad)[i], term1, term2, term3, grad[i]))
    }
  }
  
  return(grad)
}

#' Compute gradient of log-likelihood over entire trajectory
#' 
#' @param data Data frame with columns X1, X2, V1, V2, time
#' containing "true" positions and velocities at each time step
#' @param push_mat Matrix of shape (n_steps-1, 2) containing the push
#' at each time step
#' @param potential_grad_mat Matrix of shape (n_steps-1, 2) containing the
#' gradient of the potential at each time step
#' @param tau RACVM parameter tau
#' @param nu RACVM parameter nu
#' @param omega RACVM parameter omega
#' @param scheme Integration scheme
#' @param verbose If TRUE, print detailed computation steps
#' @return Numeric vector of length 3 with gradients w.r.t. tau, nu, omega
#' 
llk_gradient = function(data,push_mat,potential_grad_mat,tau,nu,omega,
                        scheme="Lie-Trotter",verbose=FALSE) {
  
  n_steps <- nrow(data)
  total_grad <- matrix(0,ncol=3,nrow=n_steps-1)
  
  for (j in 1:(n_steps-1)) {
    
    U_prev <- as.numeric(data[j,c("X1","X2","V1","V2")])
    U_next <- as.numeric(data[j+1,c("X1","X2","V1","V2")])
    delta <- as.numeric(data[j+1,"time"] - data[j,"time"])
    push <- push_mat[j,]
    push_next <- push_mat[j,]
    potential_grad <- potential_grad_mat[j,]
    potential_grad_next <- potential_grad_mat[j,]
    if (verbose) {
      cat("---- Step", j, "----\n")
      cat("delta:", delta, "\n")
      cat("U_prev:", U_prev, "\n")
      cat("U_next:", U_next, "\n")
      cat("push:", push, "\n")
      cat("potential_grad:", potential_grad, "\n")
    }
    step_grad <- llk_gradient_one_step(U_next,U_prev,delta,push,
                                       potential_grad,
                                       tau,nu,omega,scheme,
                                       push_next,potential_grad_next,
                                       verbose)
    total_grad[j,] <- step_grad
  }

  
  return(total_grad)
}
