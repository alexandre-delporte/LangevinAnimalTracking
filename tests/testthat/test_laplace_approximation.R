#' Test Suite for Laplace Approximation Student-t Function
#' 
#' This script provides comprehensive tests to verify that the Laplace
#' approximation is working correctly.

library(mvtnorm)
library(testthat)
devtools::load_all("../../LangevinAnimalTracking")

# Source your function
# source("laplace_approximation_student.R")

#' Helper function: Log-posterior density
#' @param U State vector (length 4)
#' @param y Observation (length 2)
#' @param mean Prior mean (length 4)
#' @param Q_inv Inverse prior covariance
#' @param M Observation matrix (2x4)
#' @param scale Student-t scale parameter
#' @param df Student-t degrees of freedom
log_posterior <- function(U, y, mean, Q_inv, M, scale, df) {
  r <- y - M %*% U
  r_squared <- sum(r^2)
  
  # Student-t log-likelihood (sum of two independent univariate t)
  log_lik <- -(df + 1) / 2 * sum(log(1 + r^2 / (df * scale^2)))
  
  # Gaussian prior
  log_prior <- -0.5 * t(U - mean) %*% Q_inv %*% (U - mean)
  
  return(as.numeric(log_lik + log_prior))
}

#' Helper function: Numerical gradient
#' @param U State vector
#' @param ... Additional arguments passed to log_posterior
numerical_gradient <- function(U, ..., eps = 1e-8) {
  grad <- numeric(length(U))
  for (i in seq_along(U)) {
    U_plus <- U_minus <- U
    U_plus[i] <- U[i] + eps
    U_minus[i] <- U[i] - eps
    grad[i] <- (log_posterior(U_plus, ...) - log_posterior(U_minus, ...)) / (2 * eps)
  }
  return(grad)
}

#' Helper function: Numerical Hessian
#' @param U State vector
#' @param ... Additional arguments passed to log_posterior
numerical_hessian <- function(U, ..., eps = 1e-5) {
  n <- length(U)
  H <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      U_pp <- U_pm <- U_mp <- U_mm <- U
      
      U_pp[i] <- U_pp[i] + eps
      U_pp[j] <- U_pp[j] + eps
      
      U_pm[i] <- U_pm[i] + eps
      U_pm[j] <- U_pm[j] - eps
      
      U_mp[i] <- U_mp[i] - eps
      U_mp[j] <- U_mp[j] + eps
      
      U_mm[i] <- U_mm[i] - eps
      U_mm[j] <- U_mm[j] - eps
      
      H[i,j] <- (log_posterior(U_pp, ...) - log_posterior(U_pm, ...) - 
                   log_posterior(U_mp, ...) + log_posterior(U_mm, ...)) / (4 * eps^2)
    }
  }
  
  return((H + t(H)) / 2)  # Force symmetry
}

# ==============================================================================
# TEST 1: Convergence Test
# ==============================================================================

test_that("Laplace approximation converges", {
  set.seed(123)
  
  # Setup
  mean <- c(25, 5, 0.5, 0.3)
  Q <- diag(c(0.1, 0.1, 0.5, 0.5))
  M <- cbind(diag(2), matrix(0, 2, 2))
  y <- c(25.3, 5.2)
  scale <- 0.2
  df <- 3
  
  # Run Laplace approximation
  result <- laplace_approximation_student(y, mean, Q, M, scale, df, return_cov = TRUE)
  
  # Check convergence
  expect_true(result$converged, info = "Algorithm should converge")
  expect_true(result$n_iter < 20, info = "Should converge in < 10 iterations")
  
  cat("\n✓ Test 1 PASSED: Convergence test\n")
  cat("  Converged:", result$converged, "\n")
  cat("  Iterations:", result$n_iter, "\n")
})

# ==============================================================================
# TEST 2: Mode is Actually a Maximum
# ==============================================================================

test_that("Found mode is a local maximum", {
  set.seed(123)
  
  # Setup
  mean <- c(25, 5, 0.5, 0.3)
  Q <- diag(c(0.1, 0.1, 0.5, 0.5))
  Q_inv <- solve(Q)
  M <- cbind(diag(2), matrix(0, 2, 2))
  y <- c(25.3, 5.2)
  scale <- 0.2
  df <- 3
  
  # Run Laplace approximation
  result <- laplace_approximation_student(y, mean, Q, M, scale, df)
  
  # Evaluate log-posterior at mode
  log_post_mode <- log_posterior(result$mean, y, mean, Q_inv, M, scale, df)
  
  # Evaluate at random perturbations
  n_tests <- 200
  all_lower <- TRUE
  
  for (i in 1:n_tests) {
    perturbation <- rnorm(4, 0, 0.001)
    U_perturb <- result$mean + perturbation
    log_post_perturb <- log_posterior(U_perturb, y, mean, Q_inv, M, scale, df)
    print(log_post_perturb)
    if (log_post_perturb > log_post_mode + 1e-6) {
      all_lower <- FALSE
      cat("  Found higher point! Difference:", log_post_perturb - log_post_mode, "\n")
    }
  }
  
  expect_true(all_lower, info = "All nearby points should have lower log-posterior")
  
  cat("\n✓ Test 2 PASSED: Mode verification\n")
  cat("  Log-posterior at mode:", log_post_mode, "\n")
  cat("  All", n_tests, "random perturbations had lower log-posterior\n")
})

# ==============================================================================
# TEST 3: Gradient is Zero at Mode
# ==============================================================================

test_that("Gradient is zero at the mode", {
  set.seed(123)
  
  # Setup
  mean <- c(25, 5, 0.5, 0.3)
  Q <- diag(c(0.1, 0.1, 0.5, 0.5))
  Q_inv <- solve(Q)
  M <- cbind(diag(2), matrix(0, 2, 2))
  y <- c(25.3, 5.2)
  scale <- 0.2
  df <- 3
  
  # Run Laplace approximation
  result <- laplace_approximation_student(y, mean, Q, M, scale, df)
  
  # Compute numerical gradient at mode
  grad_numerical <- numerical_gradient(result$mean, y, mean, Q_inv, M, scale, df)
  
  # Gradient should be close to zero
  grad_norm <- sqrt(sum(grad_numerical^2))
  
  expect_true(grad_norm < 1e-5, 
              info = paste("Gradient norm should be < 1e-5, got", grad_norm))
  
  cat("\n✓ Test 3 PASSED: Gradient at mode is zero\n")
  cat("  Gradient norm:", sprintf("%.2e", grad_norm), "\n")
  cat("  Gradient:", sprintf("%.2e", grad_numerical), "\n")
})

# ==============================================================================
# TEST 4: Hessian Matches Numerical Hessian
# ==============================================================================

test_that("Analytical Hessian matches numerical Hessian", {
  set.seed(123)
  
  # Setup
  mean <- c(25, 5, 0.5, 0.3)
  Q <- diag(c(0.1, 0.1, 0.5, 0.5))
  Q_inv <- solve(Q)
  M <- cbind(diag(2), matrix(0, 2, 2))
  y <- c(25.3, 5.2)
  scale <- 0.2
  df <- 3
  
  # Run Laplace approximation
  result <- laplace_approximation_student(y, mean, Q, M, scale, df, return_cov = TRUE)
  
  # Analytical Hessian (negative)
  H_analytical <- -solve(result$cov)
  
  # Numerical Hessian
  H_numerical <- numerical_hessian(result$mean, y, mean, Q_inv, M, scale, df)
  
  # Compare
  max_diff <- max(abs(H_analytical - H_numerical))
  rel_error <- max_diff / max(abs(H_numerical))
  
  expect_true(max_diff < 1e-3, 
              info = paste("Max Hessian difference should be < 1e-3, got", max_diff))
  
  cat("\n✓ Test 4 PASSED: Hessian verification\n")
  cat("  Max absolute difference:", sprintf("%.2e", max_diff), "\n")
  cat("  Relative error:", sprintf("%.2e", rel_error), "\n")
})

# ==============================================================================
# TEST 5: Covariance is Positive Definite
# ==============================================================================

test_that("Posterior covariance is positive definite", {
  set.seed(123)
  
  # Setup
  mean <- c(25, 5, 0.5, 0.3)
  Q <- diag(c(0.1, 0.1, 0.5, 0.5))
  M <- cbind(diag(2), matrix(0, 2, 2))
  y <- c(25.3, 5.2)
  scale <- 0.2
  df <- 3
  
  # Run Laplace approximation
  result <- laplace_approximation_student(y, mean, Q, M, scale, df, return_cov = TRUE)
  
  # Check eigenvalues
  eigenvalues <- eigen(result$cov, only.values = TRUE)$values
  
  expect_true(all(eigenvalues > 0), 
              info = "All eigenvalues should be positive")
  
  # Check Cholesky exists
  expect_true(!is.null(result$chol), 
              info = "Cholesky decomposition should exist")
  
  # Verify Cholesky
  cov_reconstructed <- t(result$chol) %*% result$chol
  max_diff <- max(abs(cov_reconstructed - result$cov))
  
  expect_true(max_diff < 1e-10,
              info = "Cholesky should reconstruct covariance")
  
  cat("\n✓ Test 5 PASSED: Positive definiteness\n")
  cat("  Min eigenvalue:", min(eigenvalues), "\n")
  cat("  Max eigenvalue:", max(eigenvalues), "\n")
  cat("  Cholesky reconstruction error:", sprintf("%.2e", max_diff), "\n")
})

# ==============================================================================
# TEST 6: Consistency Across Different df Values
# ==============================================================================

test_that("Works for different degrees of freedom", {
  set.seed(123)
  
  # Setup
  mean <- c(25, 5, 0.5, 0.3)
  Q <- diag(c(0.1, 0.1, 0.5, 0.5))
  M <- cbind(diag(2), matrix(0, 2, 2))
  y <- c(25.3, 5.2)
  scale <- 0.2
  
  df_values <- c(3, 5, 10, 20, 100)
  
  cat("\n✓ Test 6: Testing different df values\n")
  
  for (df in df_values) {
    result <- laplace_approximation_student(y, mean, Q, M, scale, df)
    
    expect_true(result$converged, 
                info = paste("Should converge for df =", df))
    
    cat(sprintf("  df = %3d: converged in %d iterations\n", df, result$n_iter))
  }
  
  cat("  All df values converged successfully\n")
})

# ==============================================================================
# TEST 7: Sampling and Density Evaluation
# ==============================================================================

test_that("Can sample from proposal and evaluate density", {
  set.seed(123)
  
  # Setup
  mean <- c(24, 4.9, 0.5, 0.3)
  Q <- diag(c(0.1, 0.1, 0.5, 0.5))
  M <- cbind(diag(2), matrix(0, 2, 2))
  y <- c(25.3, 5.2)
  scale <- 0.2
  df <- 3
  
  # Run Laplace approximation
  result <- laplace_approximation_student(y, mean, Q, M, scale, df, return_cov = TRUE)
  
  # Sample from proposal
  n_samples <- 10000
  samples <- matrix(0, n_samples, 4)
  
  for (i in 1:n_samples) {
    Z <- rnorm(4)
    samples[i, ] <- result$mean + result$chol %*% Z
  }
  
  # Check sample mean is close to mode
  sample_mean <- colMeans(samples)
  mean_diff <- sqrt(sum((sample_mean - result$mean)^2))
  
  expect_true(mean_diff < 0.1,
              info = "Sample mean should be close to mode")
  
  # Check sample covariance is close to analytical covariance
  sample_cov <- cov(samples)
  cov_diff <- max(abs(sample_cov - result$cov))
  
  expect_true(cov_diff < 0.05,
              info = "Sample covariance should be close to analytical")
  
  cat("\n✓ Test 7 PASSED: Sampling test\n")
  cat("  Sample mean error:", sprintf("%.3f", mean_diff), "\n")
  cat("  Sample cov error:", sprintf("%.3f", cov_diff), "\n")
})

# ==============================================================================
# TEST 8: Comparison with Naive Variance-Matching
# ==============================================================================

test_that("Laplace approximation improves over naive approach", {
  set.seed(123)
  
  # Setup
  mean <- c(25, 5, 0.5, 0.3)
  Q <- diag(c(0.1, 0.1, 0.5, 0.5))
  Q_inv <- solve(Q)
  M <- cbind(diag(2), matrix(0, 2, 2))
  y <- c(25.3, 5.2)
  scale <- 0.2
  df <- 3
  
  # Laplace approximation
  result_laplace <- laplace_approximation_student(y, mean, Q, M, scale, df, 
                                                  return_cov = TRUE)
  
  # Naive variance-matching approach
  P2 <- (df - 2) / (df * scale^2) * diag(2)
  Sigma_inv_naive <- Q_inv + t(M) %*% P2 %*% M
  Sigma_naive <- solve(Sigma_inv_naive)
  m_naive <- Sigma_naive %*% (Q_inv %*% mean + t(M) %*% P2 %*% y)
  
  # Compare modes
  mode_diff <- sqrt(sum((result_laplace$mean - m_naive)^2))
  
  cat("\n✓ Test 8: Comparison with naive approach\n")
  cat("  Laplace mode:", sprintf("%.3f", result_laplace$mean[1:2]), "\n")
  cat("  Naive mode:  ", sprintf("%.3f", m_naive[1:2]), "\n")
  cat("  Difference:  ", sprintf("%.3f", mode_diff), "\n")
  
  # The modes should be different (showing improvement)
  expect_true(mode_diff > 1e-6,
              info = "Laplace and naive should give different results")
  
  # Evaluate log-posterior at both modes
  log_post_laplace <- log_posterior(result_laplace$mean, y, mean, Q_inv, M, scale, df)
  log_post_naive <- log_posterior(m_naive, y, mean, Q_inv, M, scale, df)
  
  cat("  Log-posterior at Laplace mode:", sprintf("%.3f", log_post_laplace), "\n")
  cat("  Log-posterior at naive mode:  ", sprintf("%.3f", log_post_naive), "\n")
  
  # Laplace should find higher log-posterior (better approximation)
  expect_true(log_post_laplace >= log_post_naive - 1e-6,
              info = "Laplace mode should have at least as high log-posterior")
})

# ==============================================================================
# TEST 9: Edge Case - Observation Equals Prior Mean
# ==============================================================================

test_that("Handles edge case when observation equals prior mean", {
  set.seed(123)
  
  # Setup where observation = prior mean
  mean <- c(25, 5, 0.5, 0.3)
  Q <- diag(c(0.1, 0.1, 0.5, 0.5))
  M <- cbind(diag(2), matrix(0, 2, 2))
  y <- mean[1:2]  # Observation equals prior mean for position
  scale <- 0.2
  df <- 3
  
  # Run Laplace approximation
  result <- laplace_approximation_student(y, mean, Q, M, scale, df)
  
  expect_true(result$converged,
              info = "Should converge even when y = prior mean")
  
  # Mode should be close to prior mean
  diff <- sqrt(sum((result$mean - mean)^2))
  
  expect_true(diff < 0.1,
              info = "Mode should be close to prior mean when observation matches")
  
  cat("\n✓ Test 9 PASSED: Edge case - observation = prior mean\n")
  cat("  Distance from prior mean:", sprintf("%.3f", diff), "\n")
})

# ==============================================================================
# RUN ALL TESTS
# ==============================================================================

run_all_tests <- function() {
  cat("\n")
  cat("================================================================================\n")
  cat("  LAPLACE APPROXIMATION TEST SUITE\n")
  cat("================================================================================\n")
  
  test_results <- test_dir(".", reporter = "summary")
  
  cat("\n")
  cat("================================================================================\n")
  cat("  TEST SUMMARY\n")
  cat("================================================================================\n")
  print(test_results)
  
  return(test_results)
}

# To run all tests:
# run_all_tests()