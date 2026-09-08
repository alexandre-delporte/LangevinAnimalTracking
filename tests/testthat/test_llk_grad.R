


test_that("Gradient of one-step log-likelihood w.r.t. movement parameters
          matches finite differences", {
  
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
      U_next, U_prev, dt, push, potential_grad,
      tau = tau, nu = nu, omega = omega)
  
    strang_grad<-llk_gradient_one_step(U_next, U_prev, dt, push, potential_grad,
                                    tau, nu, omega, 
                                    scheme = "Strang",
                                    push_next = push_next, potential_grad_next = 
                                    potential_grad_next) 
  
    # Numerical log-likelihood
    lie_trotter_llk <- function(tau, nu, omega) {
    
      Q <- RACVM_cov(tau, nu, omega, dt)
      T_mat <- RACVM_link(tau, omega, dt)
    
      g_U <- c(0, 0, push + potential_grad)
      mu <- as.vector(T_mat %*% (U_prev - dt * g_U))
      r <- U_next - mu
    
      log_det_Q <- determinant(Q, logarithm = TRUE)$modulus[1]
    
      -0.5 * log_det_Q - 0.5 * drop(t(r) %*% solve(Q, r))
    }
  
    strang_llk<-function(tau,nu,omega) {
    
      Q <- RACVM_cov(tau, nu, omega, dt)
      T_mat <- RACVM_link(tau, omega, dt)
    
      g_U <- c(0, 0, push + potential_grad)
      mu <- as.vector(T_mat %*% (U_prev - dt/2 * g_U))
      V_tilde <- V_next + (dt/2) * (push_next + potential_grad_next)
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


test_that("Gradient of one-step log-likelihood w.r.t. potential parameters 
          matches finite differences", {

  tau   <- sde_params$tau
  nu    <- sde_params$nu
  omega <- sde_params$omega

  dt     <- 1/60
  Tmax   <- 2
  N      <- Tmax / dt
  q      <- 0.8
  lambda <- dt^q

  data <- simulate_2D_trajectory(1, sde_params, potential_params,
                                 error_params = list(sigma_obs = 1), error_dist = "normal",
                                 polygon, lambda, U0, N, dt, scheme = "Lie-Trotter",
                                 split_around_fixed_point = FALSE, seed = 2025)

  eps <- 1e-6
  J   <- length(potential_params$alpha)

  for (i in 2:nrow(data)) {
    U_prev <- as.numeric(data[i-1, c("X1","X2","V1","V2")])
    U_next <- as.numeric(data[i,   c("X1","X2","V1","V2")])
    X_prev <- U_prev[1:2]
    X_next <- U_next[1:2]
    V_next <- U_next[3:4]

    push        <- compute_push(X_prev, polygon, lambda)
    push_next   <- compute_push(X_next, polygon, lambda)
    grad_H_prev <- mix_gaussian_grad_cpp(X_prev, x_star, params, exclude = integer(0))
    grad_H_next <- mix_gaussian_grad_cpp(X_next, x_star, params, exclude = integer(0))

    lt_grad <- llk_gradient_one_step(
      U_next, U_prev, dt, push, grad_H_prev,
      tau, nu, omega, scheme = "Lie-Trotter",
      potential_params = potential_params, x_star = x_star,
      estimate_centers = TRUE
    )
    strang_grad <- llk_gradient_one_step(
      U_next, U_prev, dt, push, grad_H_prev,
      tau, nu, omega, scheme = "Strang",
      push_next = push_next, potential_grad_next = grad_H_next,
      potential_params = potential_params, x_star = x_star,
      estimate_centers = TRUE
    )

    lt_llk_xi <- function(pp) {
      Q     <- RACVM_cov(tau, nu, omega, dt)
      T_mat <- RACVM_link(tau, omega, dt)
      gH    <- mix_gaussian_grad_cpp(X_prev, x_star, list(alpha = pp$alpha, B = pp$B), exclude = integer(0))
      mu    <- as.vector(T_mat %*% (U_prev - dt * c(0, 0, push + gH)))
      r     <- U_next - mu
      -0.5 * determinant(Q, logarithm = TRUE)$modulus[1] - 0.5 * drop(t(r) %*% solve(Q, r))
    }

    strang_llk_xi <- function(pp) {
      Q      <- RACVM_cov(tau, nu, omega, dt)
      T_mat  <- RACVM_link(tau, omega, dt)
      gH_j   <- mix_gaussian_grad_cpp(X_prev, x_star, list(alpha = pp$alpha, B = pp$B), exclude = integer(0))
      gH_jp1 <- mix_gaussian_grad_cpp(X_next, x_star, list(alpha = pp$alpha, B = pp$B), exclude = integer(0))
      mu     <- as.vector(T_mat %*% (U_prev - (dt/2) * c(0, 0, push + gH_j)))
      r      <- c(X_next, V_next + (dt/2) * (push_next + gH_jp1)) - mu
      -0.5 * determinant(Q, logarithm = TRUE)$modulus[1] - 0.5 * drop(t(r) %*% solve(Q, r))
    }

    for (k in 1:J) {
      pp_fwd <- potential_params; pp_fwd$alpha[k] <- potential_params$alpha[k] + eps
      pp_bwd <- potential_params; pp_bwd$alpha[k] <- potential_params$alpha[k] - eps

      expect_equal(as.numeric(lt_grad[paste0("alpha_", k)]),
                   (lt_llk_xi(pp_fwd) - lt_llk_xi(pp_bwd)) / (2 * eps),
                   tolerance = 1e-5, scale = 1)
      expect_equal(as.numeric(strang_grad[paste0("alpha_", k)]),
                   (strang_llk_xi(pp_fwd) - strang_llk_xi(pp_bwd)) / (2 * eps),
                   tolerance = 1e-5, scale = 1)
    }

    # b=1 -> logL11=(1,1), b=2 -> L21=(2,1), b=3 -> logL22=(2,2)
    # Diagonal entries (b=1,3) are on the log scale: perturb log(L_ii) by eps.
    L_names  <- c("logL11", "L21", "logL22")
    is_log   <- c(TRUE, FALSE, TRUE)
    row_idx  <- c(1, 2, 2)
    col_idx  <- c(1, 1, 2)

    for (k in 1:J) {
      Lk <- t(chol(potential_params$B[[k]]))

      for (b in 1:3) {
        ri <- row_idx[b]; ci <- col_idx[b]
        Lk_fwd <- Lk; Lk_bwd <- Lk
        if (is_log[b]) {
          Lk_fwd[ri, ci] <- exp(log(Lk[ri, ci]) + eps)
          Lk_bwd[ri, ci] <- exp(log(Lk[ri, ci]) - eps)
        } else {
          Lk_fwd[ri, ci] <- Lk[ri, ci] + eps
          Lk_bwd[ri, ci] <- Lk[ri, ci] - eps
        }

        pp_fwd <- potential_params; pp_fwd$B[[k]] <- Lk_fwd %*% t(Lk_fwd)
        pp_bwd <- potential_params; pp_bwd$B[[k]] <- Lk_bwd %*% t(Lk_bwd)

        param_name <- paste0(L_names[b], "_", k)
        expect_equal(as.numeric(lt_grad[param_name]),
                     (lt_llk_xi(pp_fwd) - lt_llk_xi(pp_bwd)) / (2 * eps),
                     tolerance = 1e-5, scale = 1)
        expect_equal(as.numeric(strang_grad[param_name]),
                     (strang_llk_xi(pp_fwd) - strang_llk_xi(pp_bwd)) / (2 * eps),
                     tolerance = 1e-5, scale = 1)
      }
    }

    # Gradient w.r.t. attraction centres x_star
    lt_llk_xstar <- function(xs) {
      Q     <- RACVM_cov(tau, nu, omega, dt)
      T_mat <- RACVM_link(tau, omega, dt)
      gH    <- mix_gaussian_grad_cpp(X_prev, xs, potential_params, exclude = integer(0))
      mu    <- as.vector(T_mat %*% (U_prev - dt * c(0, 0, push + gH)))
      r     <- U_next - mu
      -0.5 * determinant(Q, logarithm = TRUE)$modulus[1] - 0.5 * drop(t(r) %*% solve(Q, r))
    }

    strang_llk_xstar <- function(xs) {
      Q      <- RACVM_cov(tau, nu, omega, dt)
      T_mat  <- RACVM_link(tau, omega, dt)
      gH_j   <- mix_gaussian_grad_cpp(X_prev, xs, potential_params, exclude = integer(0))
      gH_jp1 <- mix_gaussian_grad_cpp(X_next, xs, potential_params, exclude = integer(0))
      mu     <- as.vector(T_mat %*% (U_prev - (dt/2) * c(0, 0, push + gH_j)))
      r      <- c(X_next, V_next + (dt/2) * (push_next + gH_jp1)) - mu
      -0.5 * determinant(Q, logarithm = TRUE)$modulus[1] - 0.5 * drop(t(r) %*% solve(Q, r))
    }

    for (k in 1:J) {
      for (l in 1:2) {
        xs_fwd <- x_star; xs_fwd[k, l] <- xs_fwd[k, l] + eps
        xs_bwd <- x_star; xs_bwd[k, l] <- xs_bwd[k, l] - eps

        param_name <- paste0("xstar_", k, "_", l)
        expect_equal(as.numeric(lt_grad[param_name]),
                     (lt_llk_xstar(xs_fwd) - lt_llk_xstar(xs_bwd)) / (2 * eps),
                     tolerance = 1e-5, scale = 1)
        expect_equal(as.numeric(strang_grad[param_name]),
                     (strang_llk_xstar(xs_fwd) - strang_llk_xstar(xs_bwd)) / (2 * eps),
                     tolerance = 1e-5, scale = 1)
      }
    }
  }
})
