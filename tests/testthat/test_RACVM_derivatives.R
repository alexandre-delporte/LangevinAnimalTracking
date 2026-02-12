

test_that("RACVM analytical derivatives match finite differences", {
  
  h <- 0.02
  tau <- 1.0
  nu <- 5
  omega <- 0.1
  epsilon <- 1e-8
  
  Q_analytical <- dRACVM_cov(tau, nu, omega, h)
  T_analytical <- dRACVM_link(tau, omega, h)
  
  # dQ/dtau
  dQ_dtau_num <- (RACVM_cov(tau + epsilon, nu, omega, h) -
                    RACVM_cov(tau - epsilon, nu, omega, h)) / (2 * epsilon)
  
  expect_lt(max(abs(Q_analytical$tau - dQ_dtau_num)), 1e-5)
  
  # dQ/dnu
  dQ_dnu_num <- (RACVM_cov(tau, nu + epsilon, omega, h) -
                   RACVM_cov(tau, nu - epsilon, omega, h)) / (2 * epsilon)
  
  expect_lt(max(abs(Q_analytical$nu - dQ_dnu_num)), 1e-5)
  
  # dQ/domega
  dQ_domega_num <- (RACVM_cov(tau, nu, omega + epsilon, h) -
                      RACVM_cov(tau, nu, omega - epsilon, h)) / (2 * epsilon)
  
  expect_lt(max(abs(Q_analytical$omega - dQ_domega_num)), 1e-5)
  
  # dT/dtau
  dT_dtau_num <- (RACVM_link(tau + epsilon, omega, h)-
               RACVM_link(tau - epsilon, omega, h)) / (2 * epsilon)
  
  expect_lt(max(abs(T_analytical$tau - dT_dtau_num)), 1e-5)
  
  # dT/domega
  dT_domega_num <- (RACVM_link(tau, omega + epsilon, h) -
                            RACVM_link(tau, omega - epsilon, h)) / (2 * epsilon)
  expect_lt(max(abs(T_analytical$omega - dT_domega_num)), 1e-5)
  
})

