test_that("Gradient of one-step log-likelihood matches finite differences", {
  
  tau <- 1.0
  nu <- 4
  omega <- 1
  delta <- 1/60
  
  set.seed(123)
  U_prev <- rnorm(4)
  U_next <- rnorm(4)
  push <- c(0.1, 0.1)
  potential_grad <- c(0.2, 0.2)
  
  grad_analytical <- llk_gradient_one_step(
    U_next, U_prev, delta, push, potential_grad,
    tau = tau, nu = nu, omega = omega
  )
  
  # Numerical log-likelihood
  single_llk <- function(tau, nu, omega) {
    
    Q <- RACVM_cov(tau, nu, omega, delta)
    T_mat <- RACVM_link(tau, omega, delta)
    
    g_U <- c(0, 0, push + potential_grad)
    mu <- as.vector(T_mat %*% (U_prev - delta * g_U))
    r <- U_next - mu
    
    log_det_Q <- determinant(Q, logarithm = TRUE)$modulus[1]
    
    -0.5 * log_det_Q - 0.5 * drop(t(r) %*% solve(Q, r))
  }
  
  eps <- 1e-4
  
  grad_num <- c(
    tau   = (single_llk(tau + eps, nu, omega) -
               single_llk(tau - eps, nu, omega)) / (2 * eps),
    
    nu    = (single_llk(tau, nu + eps, omega) -
               single_llk(tau, nu - eps, omega)) / (2 * eps),
    
    omega = (single_llk(tau, nu, omega + eps) -
               single_llk(tau, nu, omega - eps)) / (2 * eps)
  )
  
  expect_equal(
    grad_analytical,
    grad_num,
    tolerance = 1e-5,
    scale = 1
  )
})
