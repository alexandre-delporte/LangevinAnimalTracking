test_that("log multivariate normal density via Cholesky matches dmvnorm", {
  
  skip_if_not_installed("mvtnorm")
  
    

  set.seed(123)
  
  for (d in c(1, 2, 5, 10)) {
    
    A <- matrix(rnorm(d * d), d, d)
    Sigma <- crossprod(A) + diag(d)
    cholSigma <- chol(Sigma)
    
    x <- rnorm(d)
    mu <- rnorm(d)
    
    ref <- mvtnorm::dmvnorm(
      x,
      mean = mu,
      sigma = Sigma,
      log = TRUE
    )
    
    r_val <- as.numeric(log_dmvnorm_chol(x, mu, cholSigma))
    cpp_val <- log_dmvnorm_chol_cpp(x, mu, cholSigma)
    
    expect_equal(
      r_val,
      ref,
      tolerance = 1e-10,
      scale = 1,
      info = paste("R implementation failed for d =", d)
    )
    
    cat("Dimension:", d, "\n")
    cat("Log density (R):", r_val, "\n")
    cat("Log density (C++):", cpp_val, "\n")
    cat("Log density (Reference):", ref, "\n")
    
    expect_equal(
      cpp_val,
      ref,
      tolerance = 1e-10,
      scale = 1,
      info = paste("C++ implementation failed for d =", d)
    )
  }
})


test_that("predictive and proposal densities via Cholesky matches dmvnorm 
          for Lie-Trotter scheme", {
  

library(sp)
source("set_up.R")
tau <- sde_params$tau
nu <- sde_params$nu
omega <- sde_params$omega

delta <- 1/60  # Time step
Tmax <- 2 # Total time
N <- Tmax / delta  # Number of time steps
q<-0.8
lambda=delta^q
error_params=list(scale=0.2,df=3)

data <- simulate_2D_trajectory(1,sde_params,potential_params,
                               error_params,error_dist="scaled_t",
                               polygon,lambda,U0,N, delta,scheme="Lie-Trotter",
                               split_around_fixed_point=FALSE,
                               seed=2025)

for (i in 2:nrow(data)) {
  
  U <- as.numeric(data[i-1,c("X1","X2","V1","V2")])
  U_pred <- as.numeric(data[i,c("X1","X2","V1","V2")])
  y <- as.numeric(data[i,c("Y1","Y2")])
  
  L<-cbind(diag(2), matrix(0,2,2))
  X<-U[1:2]
  
  push <- compute_push(X,polygon,lambda)
  U_hat<-solve_ODE(U,delta,push,potential_params,
                   ind_fixed_point=NULL)
  
  OU_solution<-solve_SDE(U_hat,delta,tau,nu,omega,potential_params,
                         ind_fixed_point=NULL)
  Q<-OU_solution$Q
  mean<-OU_solution$mean
  cholQ <- chol(Q)
  invQ  <- chol2inv(cholQ)
  
  scale<-error_params$scale;df<-error_params$df
  
  gaussian_proposal<-product_gaussian(invQ,(df-2)/df/scale^2*diag(2),
                                           mean,y,L)
  
  log_pred_density_chol <- 
    as.numeric(log_dmvnorm_chol(U_pred, mean, cholQ))
  
  log_prop_density_chol <-as.numeric(log_dmvnorm_chol(
    U_pred,gaussian_proposal$mean,
    gaussian_proposal$chol))
  
  log_pred_density<-dmvnorm(U_pred, mean = mean, sigma = Q,
                            log=TRUE)
  log_prop_density<-dmvnorm(U_pred,mean=gaussian_proposal$mean,
                            sigma=gaussian_proposal$cov,
                            log=TRUE)
  cat("Time step:", i, "\n")
  cat("Log predictive density (Cholesky):", log_pred_density_chol, "\n")
  cat("Log predictive density (dmvnorm):", log_pred_density, "\n")
  cat("Log proposal density (Cholesky):", log_prop_density_chol, "\n")
  cat("Log proposal density (dmvnorm):", log_prop_density, "\n")
  
  expect_equal(
    log_pred_density_chol,
    log_pred_density,
    tolerance = 1e-1,
    scale = 1,
    info = paste("Prediction density mismatch at time step", i)
  )
  
  expect_equal(
    log_prop_density_chol,
    log_prop_density,
    tolerance = 1e-1,
    scale = 1,
    info = paste("Proposal density mismatch at time step", i)
  )
}
})
  
  
  
  
