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
    
    r_val <- log_dmvnorm_chol(x, mu, cholSigma)
    cpp_val <- log_dmvnorm_chol_cpp(x, mu, cholSigma)
    
    expect_equal(
      r_val,
      ref,
      tolerance = 1e-10,
      scale = 1,
      info = paste("R implementation failed for d =", d)
    )
    
    expect_equal(
      cpp_val,
      ref,
      tolerance = 1e-10,
      scale = 1,
      info = paste("C++ implementation failed for d =", d)
    )
  }
})
