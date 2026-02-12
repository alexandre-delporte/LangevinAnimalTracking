


test_that("Gradient of one-step log-likelihood matches finite differences", {
  
  library(sp)
  source("set_up.R")
  tau <- sde_params$tau
  nu <- sde_params$nu
  omega <- sde_params$omega
  
  dt <- 1/60  # Time step
  Tmax <- 2 # Total time
  N <- Tmax / dt  # Number of time steps
  q<-0.8
  lambda=dt^q

  data <- simulate_2D_trajectory(1,sde_params,potential_params,
                                 error_params=list(sigma_obs=1),error_dist="normal",
                                 polygon,lambda,U0,N, dt,scheme="Lie-Trotter",
                                 split_around_fixed_point=FALSE,
                                 seed=2025)
  
  for (i in 2:nrow(data)) {
    U_prev <- as.numeric(data[i-1,c("X1","X2","V1","V2")])
    U_next <- as.numeric(data[i,c("X1","X2","V1","V2")])
    X_prev<-U_prev[1:2]
    X_next<-U_next[1:2]
    V_next<-U_next[3:4]
    
    push <- compute_push(X_prev,polygon,lambda)
    push_next<-compute_push(X_next,polygon,lambda)
    potential_grad <- mix_gaussian_grad_cpp(X_prev,x_star,params,exclude=integer(0))
    potential_grad_next<-mix_gaussian_grad_cpp(X_next,x_star,params,exclude=integer(0))
  
    lie_trotter_grad <- llk_gradient_one_step(
      U_next, U_prev, delta, push, potential_grad,
      tau = tau, nu = nu, omega = omega)
  
    strang_grad<-llk_gradient_one_step(U_next, U_prev, delta, push, potential_grad,
                                    tau, nu, omega, 
                                    scheme = "Strang",
                                    push_next = push_next, potential_grad_next = 
                                    potential_grad_next) 
  
    # Numerical log-likelihood
    lie_trotter_llk <- function(tau, nu, omega) {
    
      Q <- RACVM_cov(tau, nu, omega, delta)
      T_mat <- RACVM_link(tau, omega, delta)
    
      g_U <- c(0, 0, push + potential_grad)
      mu <- as.vector(T_mat %*% (U_prev - delta * g_U))
      r <- U_next - mu
    
      log_det_Q <- determinant(Q, logarithm = TRUE)$modulus[1]
    
      -0.5 * log_det_Q - 0.5 * drop(t(r) %*% solve(Q, r))
    }
  
    strang_llk<-function(tau,nu,omega) {
    
      Q <- RACVM_cov(tau, nu, omega, delta)
      T_mat <- RACVM_link(tau, omega, delta)
    
      g_U <- c(0, 0, push + potential_grad)
      mu <- as.vector(T_mat %*% (U_prev - delta/2 * g_U))
      V_tilde <- V_next + (delta/2) * (push_next + potential_grad_next)
      U_tilde_next <- c(X_next, V_tilde)
      r <- U_tilde_next - mu
      
      log_det_Q <- determinant(Q, logarithm = TRUE)$modulus[1]
      
      -0.5 * log_det_Q - 0.5 * drop(t(r) %*% solve(Q, r))
    }
  
    eps <- 1e-4
  
    lie_trotter_grad_num <- c(
      tau=(lie_trotter_llk(tau + eps, nu, omega) -
                 lie_trotter_llk(tau - eps, nu, omega)) / (2 * eps),
      
      nu=(lie_trotter_llk(tau, nu + eps, omega) -
                 lie_trotter_llk(tau, nu - eps, omega)) / (2 * eps),
      
      omega=(lie_trotter_llk(tau, nu, omega + eps) -
                 lie_trotter_llk(tau, nu, omega - eps)) / (2 * eps)
    )
  
    strang_grad_num <- c(tau=
      (strang_llk(tau + eps, nu, omega) -
         strang_llk(tau - eps, nu, omega)) / (2 * eps),
      
      nu=(strang_llk(tau, nu + eps, omega) -
         strang_llk(tau, nu - eps, omega)) / (2 * eps),
      
      omega=(strang_llk(tau, nu, omega + eps) -
         strang_llk(tau, nu, omega - eps)) / (2 * eps)
    )
  
    expect_equal(
      lie_trotter_grad,
      lie_trotter_grad_num,
      tolerance = 1e-5,
      scale = 1
    )
  
    expect_equal(
      strang_grad,
      strang_grad_num,
      tolerance = 1e-5,
      scale = 1
    )
  }
  
})


