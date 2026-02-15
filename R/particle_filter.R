

#' Combine two Gaussian distributions into a single Gaussian
#'
#' Computes the mean and covariance of the product of two Gaussian densities:  
#' \deqn{N(mean1, P1^{-1}) \times N(L mean2, P2^{-1}) \propto N(m, \Gamma).}
#'
#' @param P1 Precision matrix (inverse covariance) of dimension \eqn{n \times n}.
#' @param P2 Precision matrix (inverse covariance) of dimension \eqn{q \times q}.
#' @param mean1 Numeric vector of length \eqn{n}. Mean of the first Gaussian.
#' @param mean2 Numeric vector of length \eqn{q}. Mean of the second Gaussian.
#' @param M Link matrix of dimension \eqn{q \times n} that maps the state space.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{mean}}{Posterior mean vector.}
#'   \item{\code{cov}}{Posterior covariance matrix.}
#' }
product_gaussian<-function(P1,P2,mean1,mean2,M) { 
  Sigma_inv <- P1 + t(M) %*% P2 %*% M
  Sigma_inv <- (Sigma_inv + t(Sigma_inv)) / 2  # force symmetry
  
  b <- P1 %*% mean1 + t(M) %*% P2 %*% mean2
  R <- chol(Sigma_inv)
  z <- forwardsolve(t(R), b)
  m <- backsolve(R, z)
  
  chol_Sigma <- forwardsolve(t(R), diag(nrow(R)))
  
  list(
    mean = m,
    chol = chol_Sigma
  )
}

#' Laplace approximation for Student-t observation likelihood
#'
#' Finds the Gaussian approximation 
#' @param y Numeric vector of length 2. Observation.
#' @param mean Numeric vector of length 4
#' @param Q Matrix 4x4. Covariance
#' @param M Matrix 2x4. Observation matrix (typically [I_2, 0_2]).
#' @param scale Numeric. Scale parameter of Student-t.
#' @param df Numeric. Degrees of freedom of Student-t.
#' @param max_iter Integer. Maximum Newton iterations (default: 10).
#' @param tol Numeric. Convergence tolerance (default: 1e-6).
#' @param return_cov Logical. Whether to include 
#' the posterior covariance matrix in the output (default: FALSE).
#'
#' @return List with:
#'   - mean: Posterior mode (length 4)
#'   - chol: Cholesky factor of posterior covariance (for sampling)
#'   - cov: Posterior covariance matrix (optional, for diagnostics)
#'   - converged: Logical. Did Newton-Raphson converge?
#'   - n_iter: Number of iterations used.
#'
laplace_approximation_student <- function(y, mean, Q, M, scale, df,
                                          max_iter = 20, tol = 1e-12,
                                          return_cov = FALSE) {
  #browser()
  # Initialize at prior mean
  U <- as.vector(mean)
  Q_inv <- solve(Q)
  
  converged <- FALSE
  
  for (iter in 1:max_iter) {
    
    r <- as.vector(y - M %*% U)
    
    w <- (df + 1) / (df * scale^2 + r^2)        
    W <- diag(w, nrow = 2)
    
    # ---- Gradient: Student-t likelihood ----
    grad_student <- t(M) %*% (W %*% r)
    a <- (df + 1) * (df * scale^2 - r^2) / (df * scale^2 + r^2)^2
    A <- diag(a, nrow = 2)
    D2_student <- -t(M) %*% A %*% M
    
  
    grad_prior <- -Q_inv %*% (U - mean)
    D2_prior <- -Q_inv
    

    grad_total <- grad_student + grad_prior
    D2_total <- D2_student + D2_prior
    D2_total <- (D2_total + t(D2_total)) / 2  # numerical symmetry
    
 
    H_neg <- -D2_total
    
    delta <- tryCatch({
      R <- chol(H_neg)
      z <- forwardsolve(t(R), grad_total)
      backsolve(R, z)
    }, error = function(e) {
      solve(H_neg, grad_total)
    })
    
    U_new <- U + as.vector(delta)
    
    if (sqrt(sum(delta^2)) < tol) {
      converged <- TRUE
      U <- U_new
      break
    }
    
    U <- U_new
  }
  
  # ---- Final Hessian at mode ----
  r <- as.vector(y - M %*% U)
  w <- (df + 1) / (df * scale^2 + r^2)
  a <- (df + 1) * (df * scale^2 - r^2) / (df * scale^2 + r^2)^2
  
  A <- diag(a, nrow = 2)
  D2_student <- -t(M) %*% A %*% M
  D2_total <- D2_student - Q_inv
  D2_total <- (D2_total + t(D2_total)) / 2
  
  # Posterior covariance
  Sigma_post <- solve(-D2_total)
  Sigma_post <- (Sigma_post + t(Sigma_post)) / 2
  
  chol_Sigma <- tryCatch({
    chol(Sigma_post)
  }, error = function(e) {
    Sigma_post <- Sigma_post + 1e-6 * diag(nrow(Sigma_post))
    chol(Sigma_post)
  })
  
  result <- list(
    mean = as.vector(U),
    chol = t(chol_Sigma),
    converged = converged,
    n_iter = iter
  )
  
  if (return_cov) {
    result$cov <- Sigma_post
  }
  
  return(result)
}


#' Plot comparison between true target and Laplace approximation
#' 
#' @param y Numeric vector of length 2. Observation.
#' @param x_prev_mean Numeric vector of length 2. Prior mean for the state.
#' @param Qxx Matrix 2x2. Prior covariance for the state.
#' @param scale Numeric. Scale parameter of Student-t.
#' @param df Numeric. Degrees of freedom of Student-t.
#' @param laplace_mean Numeric vector of length 2. Laplace mode.
#' @param laplace_chol Matrix 2x2. Cholesky factor of Laplace covariance.
#' @param grid_radius Numeric. Radius of the grid around Laplace mode.
#' @param n_grid Integer. Number of grid points per dimension.
#' @return None. Produces a 2-panel plot.
#' 
plot_student_laplace_comparison <- function(
    y,
    x_prev_mean,
    Qxx,
    scale,
    df,
    laplace_mean,
    laplace_chol,
    grid_radius = 0.8,
    n_grid = 120) {
  
  #browser()
  stopifnot(length(y) == 2,
            length(x_prev_mean) == 2,
            length(laplace_mean) == 2)
  
  Sigma_laplace <- laplace_chol %*% t(laplace_chol)
  
  # Grid centered at Laplace mode
  x1 <- seq(laplace_mean[1] - grid_radius,
            laplace_mean[1] + grid_radius,
            length.out = n_grid)
  
  x2 <- seq(laplace_mean[2] - grid_radius,
            laplace_mean[2] + grid_radius,
            length.out = n_grid)
  
  grid <- expand.grid(x1 = x1, x2 = x2)
  
  log_prior <- mvtnorm::dmvnorm(
    as.matrix(grid),
    mean = x_prev_mean,
    sigma = Qxx,
    log = TRUE
  )
  
  log_like <-
    dscaledt(y[1], mean = grid$x1, scale = scale, df = df, log = TRUE) +
    dscaledt(y[2], mean = grid$x2, scale = scale, df = df, log = TRUE)
  
  log_target <- log_prior + log_like
  z_true <- matrix(exp(log_target), n_grid, n_grid)
  z_true <- z_true / sum(z_true)
  
  ## ---- Laplace Gaussian ----
  z_laplace <- mvtnorm::dmvnorm(
    as.matrix(grid),
    mean = laplace_mean,
    sigma = Sigma_laplace
  )
  z_laplace <- matrix(z_laplace, n_grid, n_grid)
  z_laplace <- z_laplace / sum(z_laplace)
  
  ## ---- Shared color scale ----
  zmax <- max(z_true, z_laplace)
  print(zmax)
  
  ## ---- Plot ----
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 5))
  on.exit(par(op))
  
  image(
    x1, x2, z_true,
    col = hcl.colors(100, "viridis"),
    zlim = c(0, zmax),
    main = "True target",
    xlab = expression(X[1]),
    ylab = expression(X[2])
  )
  points(x_prev_mean[1], x_prev_mean[2], col = "red", pch = 19)
  points(y[1], y[2], col = "orange", pch = 19)
  points(laplace_mean[1], laplace_mean[2], col = "white", pch = 4, lwd = 2)
  
  legend(
    "topright",
    legend = c(
      expression(" Next state mean"),
      expression("Observation"),
      "Laplace mode"),
    col = c("red", "orange", "white"),
    pch = c(19, 19, 4),
    pt.lwd = c(1, 1, 2),
    bty = "n",
    cex = 0.5
  )
  
  image(
    x1, x2, z_laplace,
    col = hcl.colors(100, "viridis"),
    zlim = c(0, zmax),
    main = "Laplace Gaussian approximation",
    xlab = expression(X[1]),
    ylab = expression(X[2])
  )
  points(x_prev_mean[1], x_prev_mean[2], col = "red", pch = 19)
  points(y[1], y[2], col = "orange", pch = 19)
  points(laplace_mean[1], laplace_mean[2], col = "white", pch = 4, lwd = 2)
  
}



#' Helper function to propagate particles for Langevin SDE
#' 
#' @param U Numeric vector of length 4. Current state vector, typically 
#'  \eqn{(x, y, v_x, v_y)} representing position and velocity.
#' @param y Numeric vector of length 2. Observed position at the next time step.
#' @param M Numeric matrix. Observation link matrix.
#' @param delta Numeric scalar. Time step size.
#' @param push Numeric vector. External push .
#' @param potential_params Parameters for mixture of gaussian potential. List
#' with elements x_star, B and alpha
#' @param tau Numeric scalar. Correlation timescale parameter.
#' @param nu Numeric scalar. Velocity  parameter.
#' @param omega Numeric scalar. Rotation parameter.
#' @param lambda penalization parameter
#' @param error_dist Character string. Distribution of observation error. 
#'   One of `"scaled_t"`, `"normal"`, or `"argos"`.
#' @param error_params List of error distribution parameters:
#'   \describe{
#'     \item{`scale`, `df`}{(for `"scaled_t"`) scale and degrees of freedom.}
#'     \item{`sigma_obs`}{(for `"normal"` and `"argos"`) observation error scale.}
#'     \item{`df`, `rho`, `a`, `p`}{(for `"argos"`) degrees of freedom, 
#'       correlation, anisotropy parameter, and mixture probability.}
#'   }
#' @param scheme Character string. Numerical scheme for propagation, 
#'   either `"Lie-Trotter"` (default) or `"Strang"`.
#' @param polygon Polygon defining the domain of interest
#' @param ind_fixed_point index of fixed point of the mixture of 
#' Gaussian potentials for splitting scheme. Default NULL
#' @param L Precomputed observation link matrix for solving the SDE. Default NULL
#' @param Q Precomputed process covariance matrix for solving the SDE. Default NULL
#' @param proposal_weight weight in the gaussian proposal mean and covariance
#' @param verbose Logical, whether to print progress messages (default: FALSE)
#' 
#'@return Numeric vector of length 4. The propagated state \eqn{U_{t+\Delta}}


propagate_langevin_particle<-function(U,y,M,delta,push,
                                      potential_params,tau,nu,omega,
                                      lambda,error_dist,error_params,
                                   scheme="Lie-Trotter",polygon,
                                   ind_fixed_point=NULL,L=NULL,Q=NULL,
                                   proposal_weight=0.5,verbose=FALSE) {
  
  
  
  X<-U[1:2]
  I2 <- diag(2)
  
  if (scheme=="Lie-Trotter")  {
    
    #ODE solution
    U_hat<-solve_ODE_cpp(U,delta,push,potential_params,
                     ind_fixed_point)
    
    #SDE mean and covariance
    OU_solution<-solve_SDE(U_hat,delta,tau,nu,omega,potential_params,
              ind_fixed_point,L,Q)
    Q<-OU_solution$Q
    mean<-OU_solution$mean
    cholQ <- chol_cpp(Q)
    invQ  <- chol2inv_cpp(cholQ)
  
    
    
    if (error_dist=="scaled_t") { 

      
      scale<-error_params$scale;df<-error_params$df
      gaussian_proposal<-product_gaussian_cpp(invQ,(df-2)/df/scale^2*I2,
                                          mean,y,M,proposal_weight)
      
    }
    
    else if (error_dist=="normal") {
      
      sigma_obs<-error_params$sigma_obs
      gaussian_proposal<-product_gaussian_cpp(invQ,1/sigma_obs^2*I2,
                                          mean,y,M,proposal_weight)
      
    } 
    else if (error_dist=="argos") {
      
      df<-error_params$df;sigma_obs<-error_params$sigma_obs
      rho<-error_params$rho;a<-error_params$a;p<-error_params$p
      
      Sigma1 <- sigma_obs^2*matrix(c(1,rho*sqrt(a),rho*sqrt(a),1),nrow=2)
      Sigma2 <- sigma_obs^2*matrix(c(1,-rho*sqrt(a),-rho*sqrt(a),1),nrow=2)
      
      S1<-df/(df-2)*Sigma1;S2<-df/(df-2)*Sigma2
      invS1<-solve(S1);invS2<-solve(S2)
    
      u<-runif(1)
      invS <- if (u < p) invS1 else invS2
      
      gaussian_proposal<-product_gaussian_cpp(invQ,invS,mean,y,M,proposal_weight)
  
    }
    
    if (verbose) {
      cat("proposal mean :",gaussian_proposal$mean,"\n")
    }
    
    # Sample next state 
    U_next <- matrix(gaussian_proposal$mean + gaussian_proposal$chol %*% rnorm(length(mean)),
                               ncol = 1)
  }
  else if (scheme=="Strang")  {
    
    alpha <- potential_params$alpha; B <- potential_params$B
    x_star <- potential_params$x_star
    
    
    #ODE solution
    U_hat<-solve_ODE_cpp(U,delta/2,push,potential_params,
                     ind_fixed_point)
    #SDE mean and covariance
    OU_solution<-solve_SDE(U_hat,delta,tau,nu,omega,potential_params,
                           ind_fixed_point,L,Q)
    Q<-OU_solution$Q
    mean<-OU_solution$mean
    
    Qxx<-Q[1:2,1:2]
    cholQxx <- chol_cpp(Qxx)
    invQxx  <- chol2inv_cpp(cholQxx)
    
    # Precompute conditional covariance for velocity
    Q_v_cond_x <- Q[3:4,3:4] - Q[3:4,1:2] %*% invQxx %*% Q[1:2,3:4]
    cholQ_v_cond_x <- chol_cpp(Q_v_cond_x)
    
    # Propagate position only
    if (error_dist == "scaled_t") {
      
      scale <- error_params$scale; df <- error_params$df
      gaussian_proposal<-product_gaussian_cpp(invQxx,(df-2)/df/scale^2*I2,
                                          mean[1:2],y,I2,proposal_weight)
      
      
    } else if (error_dist == "normal") {
      sigma_obs <- error_params$sigma_obs
      gaussian_proposal<-product_gaussian_cpp(invQxx,1/sigma_obs^2*I2,
                                          mean[1:2],y,I2,proposal_weight)
      X_next <- gaussian_proposal$mean + gaussian_proposal$chol %*% rnorm(2)
      
    } else if (error_dist == "argos") {
      df<-error_params$df;sigma_obs<-error_params$sigma_obs
      rho<-error_params$rho;a<-error_params$a;p<-error_params$p
      
      Sigma1 <- sigma_obs^2*matrix(c(1,rho*sqrt(a),rho*sqrt(a),1),nrow=2)
      Sigma2 <- sigma_obs^2*matrix(c(1,-rho*sqrt(a),-rho*sqrt(a),1),nrow=2)
      
      S1<-df/(df-2)*Sigma1;S2<-df/(df-2)*Sigma2
      invS1<-solve(S1);invS2<-solve(S2)
      
      u <- runif(1)
      invS <- if (u < p) invS1 else invS2
      gaussian_proposal<-product_gaussian_cpp(invQxx,invS,
                                          mean[1:2],y,I2,proposal_weight)
    }
    
    
    X_next <- gaussian_proposal$mean + gaussian_proposal$chol %*% rnorm(2)
    
    push_next<-compute_push(X_next,polygon,lambda)
    
    # Recompute non linear term
    if (!(is.null(ind_fixed_point))) {
    
      quad <- sapply(1:length(alpha), function(k) {
        t(X_next - x_star[k,]) %*% B[[k]] %*% (X_next - x_star[k,])
      })
      l <- ind_fixed_point
      
      B_l <- B[[l]]; alpha_l <- alpha[l]; x_star_l <- x_star[l,]
      e_l_next <- as.numeric(exp(-t(X_next-x_star_l) %*% B_l %*% (X_next-x_star_l)))
      grad_term <- push_next + mix_gaussian_grad_cpp(X_next, x_star,
                                                 list(B=B,alpha=alpha), exclude=l) +
        2*alpha_l*(e_l_next-1)*B_l %*% (X_next-x_star_l)
    } else {
      grad_term <- push_next + mix_gaussian_grad_cpp(X_next, x_star,
                                                     list(B=B,alpha=alpha), 
                                                     exclude=integer(0))
    }
    
    # Propagate velocity
    m_v_cond_x <- mean[3:4] + Q[3:4,1:2] %*% invQxx %*% (X_next - mean[1:2]) -
      delta/2 * grad_term
    V_next <- m_v_cond_x+cholQ_v_cond_x%*%rnorm(2)
    
    U_next <- matrix(c(X_next, V_next), ncol=1)
  }
  
  return(U_next)
}


#' Helper function to compute particle weights
#' 
#' @param U_pred Numeric vector of length 4. Predicted state vector
#'  \eqn{(x, y, v_x, v_y)} representing position and velocity.
#' @param U_prev Numeric vector of length 4. Previous state vector
#' @param y Numeric vector of length 2. Observed position at the next time step.
#' @param M Numeric matrix. Observation link matrix.
#' @param delta Numeric scalar. Time step size.
#' @param push Numeric vector. External push (control or forcing term).
#' @param potential_params Parameters for mixture of gaussian potential. List
#' with elements x_star, B and alpha
#' @param tau Numeric scalar. Correlation timescale parameter.
#' @param nu Numeric scalar. Velocity  parameter.
#' @param omega Numeric scalar. Rotation parameter.
#' @param error_dist Character string. Distribution of observation error. 
#'   One of `"scaled_t"`, `"normal"`, or `"argos"`.
#' @param error_params List of error distribution parameters:
#'   \describe{
#'     \item{`scale`, `df`}{(for `"scaled_t"`) scale and degrees of freedom.}
#'     \item{`sigma_obs`}{(for `"normal"` and `"argos"`) observation error scale.}
#'     \item{`df`, `rho`, `a`, `p`}{(for `"argos"`) degrees of freedom, 
#'       correlation, anisotropy parameter, and mixture probability.}
#'   }
#' @param scheme Character string. Numerical scheme for propagation, 
#'   either `"Lie-Trotter"` (default) or `"Strang"`.
#' @param ind_fixed_point index of fixed point of the mixture of
#' Gaussian potentials for splitting scheme. Default NULL
#' @param verbose Logical, whether to print progress messages (default: FALSE)
#' @param L Precomputed observation link matrix for solving the SDE. Default NULL
#' @param Q Precomputed process covariance matrix for solving the SDE. Default NULL
#' @param proposal_weight weight in the gaussian proposal mean and covariance
#' 
#'@return Numeric vector of length 4. The propagated state \eqn{U_{t+\Delta}}

compute_langevin_weight <- function(U_pred, U_prev, y,M, delta, push,
                                    potential_params = NULL, tau, nu, omega,
                                    error_dist, error_params,
                                    scheme = "Lie-Trotter",
                                    ind_fixed_point = NULL,L=NULL,Q=NULL,
                                    proposal_weight=0.5,
                                    verbose = FALSE) {
  
  I2 <- diag(2)
  if (scheme == "Lie-Trotter") {
    
    ## ODE solution
    U_hat <- solve_ODE_cpp(U_prev, delta, push,
                       potential_params, ind_fixed_point)
    
    ## SDE mean and covariance
    OU_solution <- solve_SDE(U_hat, delta, tau, nu, omega,
                             potential_params, ind_fixed_point,L,Q)
    Q    <- OU_solution$Q
    mean <- OU_solution$mean
    
    cholQ <- chol_cpp(Q)
    invQ  <- chol2inv_cpp(cholQ)
    
    
    if (error_dist == "scaled_t") {
      
      scale <- error_params$scale
      df    <- error_params$df
      
       gaussian_proposal <- product_gaussian_cpp(
         invQ,
         (df - 2) / df / scale^2 * I2,
        mean, y, M,proposal_weight)
      
      
      local_llk <-
        dscaledt(y[1], mean = U_pred[1], scale = scale, df = df, log = TRUE) +
        dscaledt(y[2], mean = U_pred[2], scale = scale, df = df, log = TRUE)
      
      log_pred_density <- 
        log_dmvnorm_chol_cpp(U_pred, mean, cholQ)
      
       log_prop_density <-log_dmvnorm_chol_cpp(
         U_pred,gaussian_proposal$mean,
         gaussian_proposal$chol)
      
      log_weight <- local_llk + log_pred_density - log_prop_density
    }
    
    else if (error_dist == "normal") {
      
      sigma_obs <- error_params$sigma_obs
      cholSigma_obs <- sigma_obs * I2
      
      gaussian_proposal <- product_gaussian_cpp(
        invQ, 1 / sigma_obs^2 * I2,
        mean, y, M,proposal_weight)
      
      local_llk <-
        log_dmvnorm_chol_cpp(y, U_pred[1:2], cholSigma_obs)
      
      log_pred_density <-
        log_dmvnorm_chol_cpp(U_pred, mean, cholQ)
      
      log_prop_density <-
        log_dmvnorm_chol_cpp(
          U_pred,
          gaussian_proposal$mean,
          gaussian_proposal$chol
        )
      
      log_weight <- local_llk + log_pred_density - log_prop_density
    }
    
    else if (error_dist == "argos") {
      
      df        <- error_params$df
      sigma_obs <- error_params$sigma_obs
      rho       <- error_params$rho
      a         <- error_params$a
      p         <- error_params$p
      
      Sigma1 <- sigma_obs^2 *
        matrix(c(1,  rho * sqrt(a),
                 rho * sqrt(a), 1), 2)
      Sigma2 <- sigma_obs^2 *
        matrix(c(1, -rho * sqrt(a),
                 -rho * sqrt(a), 1), 2)
      
      S1 <- df / (df - 2) * Sigma1
      S2 <- df / (df - 2) * Sigma2
      
      cholS1 <- chol_cpp(S1)
      cholS2 <- chol_cpp(S2)
      
      gaussian_proposal1 <- product_gaussian_cpp(invQ, chol2inv_cpp(cholS1),
                                             mean, y, M,proposal_weight)
      gaussian_proposal2 <- product_gaussian_cpp(invQ, chol2inv_cpp(cholS2),
                                             mean, y, M,proposal_weight)
      
      local_llk <- log(dmvt_mixture(y, U_pred[1:2], error_params))
      
      log_pred_density <-
        log_dmvnorm_chol_cpp(U_pred, mean, cholQ)
      
      l1 <- log(p) +
        log_dmvnorm_chol_cpp(U_pred,
                         gaussian_proposal1$mean,
                         gaussian_proposal1$chol)
      l2 <- log(1 - p) +
        log_dmvnorm_chol_cpp(U_pred,
                         gaussian_proposal2$mean,
                         gaussian_proposal2$chol)
      
      m <- max(l1, l2)
      log_prop_density <- m + log(exp(l1 - m) + exp(l2 - m))
      
      log_weight <- local_llk + log_pred_density - log_prop_density
    }
  }
  
  else if (scheme == "Strang") {
    
    ## ODE half-step
    U_hat <- solve_ODE_cpp(U_prev, delta / 2, push,
                       potential_params, ind_fixed_point)
    
    ## SDE
    OU_solution <- solve_SDE(U_hat, delta, tau, nu, omega,
                             potential_params, ind_fixed_point,L,Q)
    
    Q    <- OU_solution$Q
    mean <- OU_solution$mean
    
    Qxx      <- Q[1:2, 1:2]
    cholQxx  <- chol_cpp(Qxx)
    invQxx   <- chol2inv_cpp(cholQxx)
    
    if (error_dist == "normal") {
      
      I2 <- diag(2) 
      sigma_obs <- error_params$sigma_obs
      cholSigma_obs <- sigma_obs * I2
      
      gaussian_proposal <- product_gaussian_cpp(
         invQxx, 1 / sigma_obs^2 * I2,
         mean[1:2], y, I2,proposal_weight)
      
      
      local_llk <-
        log_dmvnorm_chol_cpp(y, U_pred[1:2], cholSigma_obs)
      
      log_pred_density <-
        log_dmvnorm_chol_cpp(U_pred[1:2], mean[1:2], cholQxx)
      
       log_prop_density <-
         log_dmvnorm_chol_cpp(
           U_pred[1:2],
           gaussian_proposal$mean,
           gaussian_proposal$chol)
      
      
      log_weight <- local_llk + log_pred_density - log_prop_density
    }
    
    else if (error_dist == "scaled_t") {
      
      scale <- error_params$scale
      df    <- error_params$df
      
       gaussian_proposal <- product_gaussian_cpp(
         invQxx,
         (df - 2) / df / scale^2 *I2,
         mean[1:2], y, I2,proposal_weight)
      
      local_llk <-
        dscaledt(y[1], mean = U_pred[1], scale = scale, df = df, log = TRUE) +
        dscaledt(y[2], mean = U_pred[2], scale = scale, df = df, log = TRUE)
      
      log_pred_density <-
        log_dmvnorm_chol_cpp(U_pred[1:2], mean[1:2], cholQxx)
      
      log_prop_density <-
         log_dmvnorm_chol_cpp(
           U_pred[1:2],
          gaussian_proposal$mean,
           gaussian_proposal$chol)
      
      
      log_weight <- local_llk + log_pred_density - log_prop_density
    }
    
    else if (error_dist == "argos") {
      
      df        <- error_params$df
      sigma_obs <- error_params$sigma_obs
      rho       <- error_params$rho
      a         <- error_params$a
      p         <- error_params$p
      
      Sigma1 <- sigma_obs^2 *
        matrix(c(1,  rho * sqrt(a),
                 rho * sqrt(a), 1), 2)
      Sigma2 <- sigma_obs^2 *
        matrix(c(1, -rho * sqrt(a),
                 -rho * sqrt(a), 1), 2)
      
      S1 <- df / (df - 2) * Sigma1
      S2 <- df / (df - 2) * Sigma2
      
      cholS1 <- chol_cpp(S1)
      cholS2 <- chol_cpp(S2)
      
      gaussian_proposal1 <- product_gaussian_cpp(invQxx, chol2inv_cpp(cholS1),
                                             mean[1:2], y, I2,proposal_weight)
      gaussian_proposal2 <- product_gaussian_cpp(invQxx, chol2inv_cpp(cholS2),
                                             mean[1:2], y, I2,proposal_weight)
      
      local_llk <- log(dmvt_mixture(y, U_pred[1:2], error_params))
      
      log_pred_density <-
        log_dmvnorm_chol_cpp(U_pred[1:2], mean[1:2], cholQxx)
      
      l1 <- log(p) +
        log_dmvnorm_chol_cpp(U_pred[1:2],
                         gaussian_proposal1$mean,
                         gaussian_proposal1$chol)
      l2 <- log(1 - p) +
        log_dmvnorm_chol_cpp(U_pred[1:2],
                         gaussian_proposal2$mean,
                         gaussian_proposal2$chol)
      
      m <- max(l1, l2)
      log_prop_density <- m + log(exp(l1 - m) + exp(l2 - m))
      
      log_weight <- local_llk + log_pred_density - log_prop_density
    }
  }
  
  if (verbose) {
    cat("local lk =", exp(local_llk),
        "pred_density =", exp(log_pred_density),
        "prop_density =", exp(log_prop_density), "\n")
  }
  exp(log_weight)
}


#' Particle filter for movement in 2D with polygon constraint based on a penalization term
#'
#' @param data A data.frame containing observations with columns: time, Y1, Y2
#' @param sde_params List of parameters: RACVM parameters mu1, mu2,tau, nu, omega 
#' OU parameter : mu1,mu2,sigma,beta
#' @param potential_params Parameters of mixture of gaussian potential. Only necessary
#' when split_around_fixed_point=TRUE
#' @param error_params List of parameters of the measurement error distribution.
#' sigma_obs for normal, scale and df for scaled_t,
#'  sigma_obs,rho, a, df and p for argos.
#' @param error_dist distribution of the measurement noise. Options are
#' normal,scaled_t for student-t as in Weensveen, 
#' and argos for multivariate student as in Brost et .al. 
#
#' @param polygon A matrix defining the boundary polygon (n x 2 coordinates)
#' @param U0 initial state
#' @param lambda penalization parameter 
#' @param num_particles Integer, number of particles to use in the filter
#' @param scheme Either Lie-Trotter or Strang for the numerical approximation of 
#' the underlying dynamic
#' @param split_around_fixed_point Whether to split the SDE adaptatively around fixed point
#' of the drift. Only works for mixture of gaussian potential. 
#' @param ESS_threshold Effective sample size threshold for resampling (default: 0.5)
#' @param proposal_weight Weight in gaussian proposal (between 0 and 1)
#' @param verbose Logical, whether to print progress messages (default: FALSE)
#'
#' @return A list containing:
#' \itemize{
#'   \item particles: Array of particle states (num_particles x num_states x N)
#'   \item weights: Matrix of particle weights (num_particles x N)
#'   \item total_weights: Vector of sum unnormalized weights at each time step
#'   \item loglik: Estimated log-likelihood of the data
#'   \item ancestors: Matrix of ancestor indices for each particle at each time step
#'   \item push_array: Array of inward push values for each particle at each time step
#'   \item ind_fixed_point_mat: Matrix of fixed point indices for each particle at each time step
#'   
#' }
#' @export
#'
particle_filter2D <- function(data,sde_params,potential_params=NULL,
                              error_params,error_dist="normal",
                              polygon,U0,lambda,num_particles,scheme="Lie-Trotter",
                              split_around_fixed_point=FALSE,ESS_threshold=0.95,
                              proposal_weight=0.5,verbose=FALSE) {
  
  
  cat("Running particle filter with params:", as.numeric(sde_params), "\n")
  # Number of observations
  N <- nrow(data)
  
  # Movement parameters
  tau<-sde_params$tau
  nu <- sde_params$nu;omega<-sde_params$omega
  
  # Obsevation link matrix
  I2 <- diag(2)
  M<-cbind(I2, matrix(0,2,2))
  
  # Potential parameters
  if (!is.null(potential_params)) {
    alpha<-potential_params$alpha;B<-potential_params$B
    x_star<-potential_params$x_star
  }
  process_cov <- rbind(
    cbind(matrix(0,2,2), matrix(0,2,2)),
    cbind(matrix(0,2,2), 4*nu^2/pi/tau*I2))
  
  # Initialize particles and weights
  particles <- array(NA, dim = c(num_particles, length(U0), N))
  weights <- matrix(0, nrow = num_particles, ncol = N)
  total_weights <- numeric(N)  # store the total unnormalized weights
  
  #initialize log likelihood
  loglik <- 0
  loglik_vector <- numeric(N-1)
  
  #initialize array to store inward push values
  push_array <- array(0, dim = c(num_particles, 2, N-1))
  
  # track ESS and resampling
  ess_history <- numeric(N)
  resampled_at <- logical(N-1)
  
  # Track index of fixed points to use in splitting
  if (split_around_fixed_point) {
    ind_fixed_point_mat<-matrix(0,nrow=num_particles,ncol=N-1)
  } else { ind_fixed_point_mat<-NULL}
  
  # Track ancestors in resampling
  ancestors <- matrix(0, nrow = num_particles, ncol = N-1)
  
  # Initialize particles with noise
  R0<- rbind(cbind(diag(0.1^2,2),matrix(0,2,2)),
             cbind(matrix(0,2,2),diag(0.5^2,2)))
  
  particles[, , 1]<-matrix(replicate(num_particles,
                                     U0+sqrt(R0)%*%rnorm(length(U0))),
                           ncol=length(U0),nrow=num_particles,byrow=TRUE)
  
  # Initialize equal weights
  weights[, 1] <- 1/num_particles
  
  total_weights[1] <- sum(weights[, 1])
  ess_history[1] <- 1 / sum(weights[,1]^2)
  
  
  if (verbose) cat("Initialization complete. Starting particle filtering...\n")
  
  # Precompute time steps
  deltas <- diff(data$time)
  
  # Precompute link and covariance matrix for solving the SDE when time step 
  # is constant
  L_list<-vector("list",num_particles)
  Q_list<-vector("list",num_particles)
  if (max(deltas)-min(deltas)<1e-6 && !split_around_fixed_point) { 
    
    L_list<-rep(list(RACVM_link(tau,omega,deltas[1])),num_particles)
    Q_list<-rep(list(RACVM_cov(tau,nu,omega,deltas[1])),num_particles)
  }
  
  # Loop over all time steps
  for (j in 1:(N - 1)) {
    
    if (verbose) cat("Time step:", j, "\n")
    
    delta <- deltas[j]
    
    if (verbose) cat("Delta (time step):", delta, "\n")
    
    # Extract current observation
    y <-as.numeric(data[j+1,c("Y1","Y2")])
    if (verbose) cat("Observation:", y, "\n")
    
    # Prediction step
    if (verbose) cat("Prediction step...\n")
    
    
    for (k in 1:num_particles) {
      
      # Get ancestor
      U_prev <- particles[k, ,j]
      
      # Compute inward push
      push<-compute_push_cpp(U_prev[1:2],polygon@coords,lambda)
      
      push_array[k, ,j] <- push
      
      # Find fixed point index for splitting, and update OU covariance and link matrices
      if (split_around_fixed_point) {
        ind_fixed_point <- choose_center_cpp(U_prev[1:2],x_star, 
                                             list(B=B,alpha=alpha))
        ind_fixed_point_mat[k,j]<- ind_fixed_point
        
        # if constant time steps, and fixed point is changed, update L and Q
        if (max(deltas)-min(deltas)<1e-6 && (j==1 || ind_fixed_point!=ind_fixed_point_mat[k,j-1])) {
          
          #cat("Particle", k, "switching to fixed point", ind_fixed_point,"at time",j,"\n")
          
          #extract parameters for the new center
          l<-ind_fixed_point
          B_l<-B[[l]];alpha_l<-alpha[l];x_star_l<-x_star[l,]
          
          
          A<-rbind(cbind(matrix(0,ncol=2,nrow=2),I2),
                   cbind(-2*alpha_l*B_l,-matrix(c(1/tau,-omega,omega,1/tau),
                                                nrow=2,byrow=TRUE)))
          #link and covariance matrix
          L_list[[k]]<-expm(A*delta)
          Q_list[[k]]<-OU_cov_exact(A, process_cov,delta,L_list[[k]])
        }
      }
      else { ind_fixed_point <- NULL
      }
      
      U_next<-propagate_langevin_particle(U_prev,y,M,delta,push,
                                              potential_params,tau,nu,omega,lambda,
                                              error_dist,error_params,
                                              scheme,polygon,ind_fixed_point,
                                              L_list[[k]],
                                              Q_list[[k]],proposal_weight,
                                              verbose=FALSE)
      
      particles[k, , j + 1]=U_next
      
    }
    
    
    if (verbose) cat("Prediction step complete.\n")
    
    # Correction step: Update weights based on observation likelihood
    if (verbose) cat("Correction step...\n")
    
    for (k in 1:num_particles) {
      
      U_pred <- particles[k, , j + 1]  # predicted particle
      U_prev<- particles[k, ,j] # previous particle
      push<-push_array[k, ,j]  #inward push force
      ind_fixed_point<-if (split_around_fixed_point) ind_fixed_point_mat[k,j]
      else NULL  #ind of fixed point for splitting
      
      weights[k, j + 1] <- compute_langevin_weight(U_pred,U_prev,y,M,
                                                       delta,push,
                                                       potential_params,tau,nu,omega,
                                                       error_dist,error_params,
                                                       scheme,
                                                       ind_fixed_point,
                                                       L_list[[k]],
                                                       Q_list[[k]],
                                                       proposal_weight,verbose=FALSE)
      
    }
    
    # Store total weight before normalization
    total_weights[j + 1] <- sum(weights[, j + 1])
    
    
    
    loglik <- loglik + log(total_weights[j + 1]) - log(num_particles)
    loglik_vector[j]<- log(total_weights[j + 1]) - log(num_particles)
    
    # Normalize weights
    weight_sum <- total_weights[j + 1]
    
    if (weight_sum != 0) {
      weights[, j + 1] <- weights[, j + 1] / weight_sum
      ess_history[j + 1] <- 1 / sum(weights[, j + 1]^2)
      if (verbose) cat("Correction step complete.\n")
      
    } else {
      warning("Sum of weights is zero at time step ", j + 1, ". Reinitializing particles...\n")
      
      # Reinitialize particles around the observed position
      particles[, , j + 1] <- replicate(num_particles,c(y,0,0)+sqrt(R0)%*%rnorm(length(U0)))
      weights[,j+1]<-1/num_particles
      ess_history[j + 1] <- 1 / sum(weights[, j + 1]^2)
    }
    
    # Resampling step
    if (j < N - 1) {
      
      ess_ratio <- ess_history[j+1] / num_particles
      
      if (ess_ratio < ESS_threshold) {
        # Systematic resampling 
        resampled_at[j] <- TRUE
        if (verbose) cat("  Resampling (ESS ratio = ", round(ess_ratio, 3), ")\n")
        
        # Vectorized systematic resampling
        u <- (runif(1) + (0:(num_particles-1))) / num_particles
        resample_indices <- findInterval(u, c(0, cumsum(weights[, j+1])))
        
      } else {
        resample_indices <- 1:num_particles
        resampled_at[j] <- FALSE
        if (verbose) cat("  No resampling (ESS ratio = ", round(ess_ratio, 3), ")\n")
      }
      particles[, , j + 1] <- particles[resample_indices, , j + 1]
      ancestors[, j] <- resample_indices
      weights[, j + 1] <- 1/num_particles 
    }
  }
  
  if (verbose) cat("Particle filtering complete.\n")
  
  # Return particles, weights, and total weights
  return(list(particles = particles, weights = weights, 
              total_weights = total_weights,loglik = loglik,
              loglik_vector=loglik_vector,
              push=push_array,
              ind_fixed_point=ind_fixed_point_mat,
              ancestors=ancestors,resampled_at=resampled_at,
              ess_history=ess_history))
}



#' Get estimate of the trajectory (conditional mean) from particle filter 
#' @param PF_results List as returned by \code{particle_filter_2D}
#' @return a matrix of estimated filtered states based on conditional mean
#' @export
get_PF_mean<-function(PF_results) {

  particles <- PF_results$particles
  weights <- PF_results$weights
  
  N<-dim(particles)[3]
  num_states<-dim(particles)[2]
  
  estimated_PF_trajectory <- matrix(0, nrow = N, ncol = num_states) 
  for (t in 1:nrow(estimated_PF_trajectory)) {
    estimated_PF_trajectory[t, ] <- colSums(particles[, , t] * weights[, t])
  }
  
  colnames(estimated_PF_trajectory)<-c("X1","X2","V1","V2")
  return (estimated_PF_trajectory)
}

#' Plot trajectory within polygon boundary
#' Visualizes and compares true and estimated trajectories within a specified domain
#'
#' @param U_true Matrix of true states (x1, x2) from simulation
#' @param U_est Matrix of estimated states (x1, x2) from particle filter
#' @param polygon SpatialPolygons object defining the boundary (from `sp` package)
#' @param true_line_opacity Opacity for true trajectory line (0-1, default: 0.5)
#' @param est_line_opacity Opacity for estimated trajectory line (0-1, default: 0.1)
#' @param zoom_box Optional list with elements `xmin`, `xmax`, `ymin`, `ymax` 
#' to zoom in on a specific area of the plot (default: NULL for no zoom)
#' @param labels Character vector of length 2 for 
#' legend labels of true and estimated paths
#'
#' @return A ggplot object showing:
#' \itemize{
#'   \item Polygon boundary 
#'   \item True trajectory 
#'   \item Estimated trajectory 
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom sp SpatialPolygons
compare_trajectories_in_polygon <- function(U_true, U_est, 
                                            polygon, true_line_opacity = 0.5,
                                            est_line_opacity = 0.1,
                                            zoom_box = NULL,
                                            labels=c("True path","Estimated path")) {
  
  # Extract the coordinates of the polygon
  polygon_coords <- as.data.frame(polygon@coords)
  
  # Create data frames for the true and estimated trajectory points
  U_true_df <- as.data.frame(U_true)
  colnames(U_true_df) <- c("X1", "X2")
  U_true_df$Type <- labels[1]
  
  U_est_df <- as.data.frame(U_est)
  colnames(U_est_df) <- c("X1", "X2")
  U_est_df$Type <- labels[2]
  
  # Combine data for easier legend handling
  combined_df <- rbind(U_true_df, U_est_df)
  combined_df$Type <- factor(combined_df$Type, levels = labels)
  
  colors<-c("black","orange")
  names(colors)=labels
  
  p <- ggplot() +
    geom_polygon(data = polygon_coords, aes(x = V1, y = V2), fill = NA, color = "grey40") +
    geom_point(data = combined_df, aes(x = X1, y = X2, color = Type), size = 1, alpha = 0.2) +
    geom_path(data = U_true_df, aes(x = X1, y = X2, color = Type), size = 1, alpha = true_line_opacity) +
    geom_path(data = U_est_df, aes(x = X1, y = X2, color = Type), size = 1, alpha = est_line_opacity) +
    scale_color_manual(values = colors,
                       name = "Trajectory") +
    labs(title ="",
         x = "X1",
         y = "X2") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Apply zoom if zoom_box is provided
  if (!is.null(zoom_box)) {
    p <- p + 
      coord_cartesian(xlim = c(zoom_box$xmin, zoom_box$xmax),
                      ylim = c(zoom_box$ymin, zoom_box$ymax))
  }
  
  print(p)
  return(p)
}

#' Plot cloud of particles to see where particles are landing 
#' at a specific time step
#' @param t0 time 
#' @param data data frame of observations with columns Y1 and Y2
#' @param particles array of particles as returned by \code{particle_filter2D}
#' @param weights matrix of particle weights as returned by \code{particle_filter2D}
#' @param polygon SpatialPolygons object defining the boundary (from `sp` package`)
#' @export
#' @import ggplot2
#' @importFrom sp SpatialPolygons
#'
particle_cloud<-function(t0,data,particles,weights,polygon=NULL) {
  
  dfp <- data.frame(x = particles[,1,t0], y = particles[,2,t0], w=weights[,t0])
  p<-ggplot() + geom_point(data=dfp, aes(x,y,size=w),alpha=0.6) +
    geom_point(aes(x=data$Y1[t0], y=data$Y2[t0]), colour="red", size=3) +
    coord_equal() + ggtitle(paste("Particle cloud at t=", t0))
  if (!(is.null(polygon))) {
    polygon_coords <- as.data.frame(polygon@coords)
    p<-p+geom_polygon(data = polygon_coords, aes(x = V1, y = V2), fill = "lightblue", color = "blue",alpha=0.1)
  }
  p
}
