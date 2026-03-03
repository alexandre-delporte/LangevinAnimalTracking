
#' Compute the covariance matrix for the RACVM model
#' @param tau time scale parameter
#' @param nu scale parameter
#' @param omega angular velocity parameter
#' @param dt time step
#' @return Covariance matrix Q (4x4 matrix)


RACVM_cov <- function(tau,nu,omega,dt) {
  
  beta<-1/tau
  sigma<-2*nu/sqrt(pi*tau)
  

  C<-beta^2+omega^2
  A<-matrix(c(beta,-omega,omega,beta),nrow=2,byrow=TRUE)
  invA<-1/C*matrix(c(beta,omega,-omega,beta),nrow=2,byrow=TRUE)
  R<-matrix(c(cos(omega*dt),sin(omega*dt),-sin(omega*dt),cos(omega*dt)),byrow=TRUE,nrow=2)
  expAdt<-exp(-beta*dt)*R
  
  # Covariance of next state vector
  var_xi<-sigma^2/C*(dt+(omega^2-3*beta^2)/(2*beta*C)-exp(-2*dt*beta)/(2*beta)+
                      2*exp(-dt*beta)*(beta*cos(omega*dt)-omega*sin(omega*dt))/C)
  var_zeta<-sigma^2/(2*beta)*(1-exp(-2*dt*beta))
  cov1<-sigma^2/(2*C)*(1+exp(-2*dt*beta)-2*exp(-dt*beta)*cos(omega*dt))
  cov2<-sigma^2/C*(exp(-dt*beta)*sin(omega*dt)-omega/(2*beta)*(1-exp(-2*dt*beta)))
  Q<-matrix(c(var_xi,0,cov1,cov2,0,var_xi,-cov2,cov1,cov1,-cov2,var_zeta,0,cov2,cov1,0,var_zeta),
            nrow=4,byrow=TRUE)
  
  return(Q)   
}

#' Computation of the OU covariance matrix with the
#' approximation based on Taylor expansion (Pilipovic, 2025)
#' @param Atilde drift matrix
#' @param Gamma diffusion matrix
#' @param dt time step
#' @return Covariance matrix Q (4x4 matrix)
OU_cov_approx <- function(Atilde, Gamma,dt) {
  A2 <- Atilde %*% Atilde
  A3 <- A2 %*% Atilde
  
  term1 <- dt * Gamma
  term2 <- (dt^2 / 2) * (Atilde %*% Gamma + Gamma %*% t(Atilde))
  term3 <- (dt^3 / 6) * (A2 %*% Gamma + 2 * Atilde %*% Gamma %*% t(Atilde) +
                          Gamma %*% t(A2))
  term4 <- (dt^4 / 24) * (A3 %*% Gamma + 
                           3 * A2 %*% Gamma %*% t(Atilde) + 
                           3 * Atilde %*% Gamma %*% t(A2) + 
                           Gamma %*% t(A3))
  
  Q <- term1 + term2 + term3 + term4
  return(Q)
}

#' Computation of OU covariance matrix based on kronecker sum (Albertsen, 2019)
#' @param A drift matrix
#' @param Gamma diffusion matrix
#' @param h time step
#' @param expAh precomputed matrix exponential of A*h
#' @return Covariance matrix Q (4x4 matrix)
#' 
OU_cov_exact <- function(A, Gamma, h,expAh) {
  
  d <- nrow(A)
  M <- kronecker(A, diag(d)) + kronecker(diag(d), A)   
  vecC <- solve(M, as.vector(Gamma))          
  C <- matrix(vecC, d, d)
  Q <- expAh %*% C %*% t(expAh)-C
  
  return(Q)
}

#' Compute the link matrix for the RACVM model
#' @param tau time scale parameter
#' @param omega angular velocity parameter
#' @param dt time step
#' @return Link matrix L (4x4 matrix)
#' 
RACVM_link <- function(tau,omega,dt) {
  
  beta=1/tau
  
  #relevant matrices
  C=beta^2+omega^2
  A=matrix(c(beta,-omega,omega,beta),nrow=2,byrow=TRUE)
  invA=1/C*matrix(c(beta,omega,-omega,beta),nrow=2,byrow=TRUE)
  R=matrix(c(cos(omega*dt),sin(omega*dt),-sin(omega*dt),cos(omega*dt)),
           byrow=TRUE,nrow=2)
  expAdt=exp(-beta*dt)*R
  
  # link matrix
  L <- matrix(0, nrow = 4, ncol = 4)
  L[1:2, 1:2] <- diag(2)
  L[1:2, 3:4] <- invA%*%(diag(2)-expAdt)
  L[3:4, 1:2] <- matrix(c(0,0,0,0),nrow=2)
  L[3:4, 3:4] <- expAdt
  return(L)   
}



#' Function to project a point on the closest point on the boundary of the polygon
#'@param x Numeric vector of length 2: coordinates of the point
#'@param polygon SpatialPolygon object representing the polygon
#'@param grad logical. If TRUE, also compute the gradient of the projection
#'@return A list with elements:
#' - point: Numeric vector of length 2, coordinates of the closest point on the boundary
#' - gradient: 2x2 matrix representing the gradient of the projection (if grad=TRUE)
#' 
closest_point_on_boundary <- function(x, polygon,grad=TRUE) {
  coords <- polygon@coords
  point_to_segment_closest_point <- function(p, v, w) {
    l2 <- sum((w - v)^2)
    if (l2 == 0) return(list(point = v, t = 0))
    t_val <- max(0, min(1, sum((p - v) * (w - v)) / l2))
    projection <- v + t_val * (w - v)
    return(list(point = projection, t = t_val))
  }
  
  min_dist <- Inf
  closest_point <- NULL
  gradient <- NULL
  
  for (i in 1:(nrow(coords) - 1)) {
    v <- coords[i, ]
    w <- coords[i + 1, ]
    result <- point_to_segment_closest_point(x, v, w)
    projection <- result$point
    t_val <- result$t
    dist <- sqrt(sum((x - projection)^2))
    
    if (dist < min_dist) {
      min_dist <- dist
      closest_point <- projection
      e <- w - v
      # Check if projection is strictly inside the edge (not a vertex)
      if (t_val > 0 && t_val < 1) {
        gradient <- (e %*% t(e)) / sum(e * e) 
      } else {
        gradient <- NULL  # Non-differentiable at vertices
      }
    }
  }
  
  return(list(point = closest_point, gradient = gradient))
}

#' Function to check if a point is inside the polygon
#' @param x Numeric vector of length 2: coordinates of the point
#' @param polygon SpatialPolygon object representing the polygon
#' @return logical: TRUE if the point is inside or on the edge, FALSE otherwise
#' @importFrom sp point.in.polygon
is_point_inside_polygon <- function(x, polygon) {
  coords <- polygon@coords
  status <- point.in.polygon(x[1], x[2], coords[, 1], coords[, 2])
  return(status == 1 || status == 2)  # 1: inside, 2: on the edge
}


#' Compute penalisation term
#' @param x Numeric vector of length 2: coordinates of the point
#' @param polygon SpatialPolygon object representing the polygon
#' @param lambda penalty parameter
#' @return Numeric vector of length 2: penalisation term
compute_push <- function(x, polygon, lambda) {
  
  if (lambda==Inf) {
    return(c(0,0))
  }
  if (is_point_inside_polygon(x, polygon)) {
    return(c(0,0))
  }
  proj <- closest_point_on_boundary(x, polygon)$point
  return ((x - proj) / lambda)
}

#' Compute the gradient of mixture of Gaussian potentials
#'
#' @param x coordinates of the position
#' @param x_star matrix of coordinates of each center in the mixture (default = 0,0)
#' @param params List of parameters:
#'   - alpha: list of weights in the mixture
#'   - B:  list of covariance matrices in the mixture
#' @param exclude vector of indices of mixture components to exclude from the computation
#'
#' @return A numeric vector of length 2: (dphi/dx, dphi/dy)
#' @export
#'
mix_gaussian_grad <- function(x, x_star = matrix(c(0,0), ncol = 2), params,
                              exclude=NULL) {
  
  alpha <- params$alpha
  B_list <- params$B
  m <- length(alpha)
  
  
  grad <- c(0, 0)
  for (j in seq_len(m)) {
    if (!(j %in% exclude)) {
    diff <- x - x_star[j, ]
    quad <- as.numeric(t(diff) %*% B_list[[j]] %*% diff)
    grad <- grad + 2 * alpha[j] * (B_list[[j]] %*% diff) * exp(-quad)
    }
  }
  
  return(as.numeric(grad))
}


#' Compute value of a mixture of Gaussian potentials
#'
#' @param x coordinates of the position
#' @param x_star matrix of coordinates of each center in the mixture
#' @param params List of parameters:
#'   - alpha: list of weights in the mixture
#'   - B:  list of covariance matrices in the mixture
#' @param exclude vector of indices of mixture components to exclude from the computation
#'
#' @return A numeric vector of length 2: (dphi/dx, dphi/dy)
#' @export
mix_gaussian_potential <- function(x, x_star, params,exclude=NULL) {
  alpha <- params$alpha
  B_list <- params$B
  J <- length(alpha)
  
  val <- 0
  for (j in seq_len(J)) {
    if (!(j %in% exclude)) {
      diff <- x - x_star[j, ]
      quad <- as.numeric(t(diff) %*% B_list[[j]] %*% diff)
      val <- val - alpha[j] * exp(-quad)
    }
  }
  return(val)
}


#' Compute the Hessian of a Gaussian mixture potential
#'
#'
#' @param x Numeric vector of length 2: coordinates of the position
#' @param x_star Matrix of size m x 2 giving the coordinates of each center
#' @param params List of parameters:
#'   - alpha: numeric vector of length m, weights of the mixture
#'   - B:     list of length m, each element is a 2x2 positive-definite matrix
#' @param exclude vector of indices of mixture components to exclude 
#' from the computation
#'
#' @return A 2x2 Hessian matrix evaluated at x
#' @export
mix_gaussian_hessian <- function(x, x_star, params,exclude=NULL) {
  
  required_params <- c("alpha", "B")
  if (!all(required_params %in% names(params))) {
    stop("params must include 'alpha' and 'B'")
  }
  
  alpha <- params$alpha
  B_list <- params$B
  J <- length(alpha)
  
  if (nrow(x_star) != J) {
    stop("x_star must have as many rows as length(alpha)")
  }
  if (length(B_list) != J) {
    stop("B must be a list of length m")
  }
  
  H <- matrix(0, nrow = 2, ncol = 2)
  for (j in seq_len(J)) {
    
    if (!(j %in% exclude)) {
      diff <- x - x_star[j, ]
      B <- B_list[[j]]
      quad <- as.numeric(t(diff) %*% B %*% diff)
      e <- exp(-quad)
      u <- B %*% diff
      H <- H + 2 * alpha[j] * (B * e - 2 * B%*%(u %*% t(u))%*%B* e)
    }
  }
  
  return(H)
}


#' Choose the center with the largest gradient norm
#'
#' @param x coordinates of the position (numeric vector of length 2)
#' @param x_star matrix of coordinates of each center in the mixture (J x 2)
#' @param params List of parameters:
#'   - alpha: list of weights in the mixture
#'   - B:  list of precision matrices in the mixture
#'
#' @return index of the chosen center
#' @export
choose_center <- function(x, x_star, params) {
  alpha <- params$alpha
  B_list <- params$B
  J <- length(alpha)
  
  log_grad2 <- function(x, alpha, center, B) {
    diff <- x - center
    quad <- as.numeric(t(diff) %*% B %*% diff)             
    q <- as.numeric(t(diff) %*% t(B) %*% (B %*% diff))     
    return(2*log(alpha) - 2*quad + log(q)) 
  }
  
  L_vals <- vapply(seq_len(J), function(j) {
    log_grad2(x, alpha[j], x_star[j, ], B_list[[j]])
  }, numeric(1))
  
  j_max <- which.max(L_vals)
  
  return(j_max)
}


#' Function to compute RMSE between true and filtered positions
#' @param true matrix of true positions (n x 2)
#' @param filtered matrix of filtered positions (n x 2)
#' @return RMSE value (numeric)
#' @export
rmse=function(true,filtered) {
  
  return (mean((sqrt((true[,1]-filtered[,1])^2+(true[,2]-filtered[,2])^2))))
}

#' Function to compute max error between true and filtered positions
#' @param true matrix of true positions (n x 2)
#' @param filtered matrix of filtered positions (n x 2)
#' @return max error value (numeric)
#' @export
max_error=function(true,filtered) {
  
  return (max((sqrt((true[,1]-filtered[,1])^2+(true[,2]-filtered[,2])^2))))
}

#' Log-density of multivariate normal using Cholesky decomposition
#' @param x Numeric vector of length d: point where to evaluate the density
#' @param mean Numeric vector of length d: mean of the distribution
#' @param cholSigma Upper Cholesky factor of the covariance matrix Sigma (d x d)
#' @return Log-density value (numeric)
#' 
log_dmvnorm_chol <- function(x, mean, cholSigma) {
  # x, mean: numeric vectors of length d
  # cholSigma: upper Cholesky factor of Sigma (d x d)
  d <- length(x)
  diff <- x - mean
  v <- forwardsolve(t(cholSigma), diff)
  
  quad <- crossprod(v)
  logdet <- 2 * sum(log(diag(cholSigma)))
  
  -0.5 * (quad + logdet + d * log(2 * pi))
}


#' Log-density of multivariate t-distribution using Cholesky decomposition
#' @param x Numeric vector of length d: point where to evaluate the density
#' @param mean Numeric vector of length d: mean of the distribution
#' @param chol Upper Cholesky factor of the scale matrix (d x d)
#' @param df Degrees of freedom
#' @return Log-density value (numeric)
log_dmvt_chol <- function(x, mean, chol, df) {
  d  <- length(x)
  z  <- backsolve(chol, x - mean, transpose = TRUE)
  delta <- sum(z^2)
  
  logdet <- 2 * sum(log(diag(chol)))
  
  lgamma((df + d) / 2) -
    lgamma(df / 2) -
    0.5 * (d * log(df * pi) + logdet) -
    0.5 * (df + d) * log1p(delta / df)
}


