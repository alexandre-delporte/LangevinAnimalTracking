



#' Solve deterministic ODE part of the Langevin splitting scheme
#'
#' Computes the ODE flow for one time step of the splitting scheme, 
#' either naively or around a fixed point of a Gaussian mixture potential.
#'
#' @param U Numeric vector of length 4. Current state vector, typically 
#'   \eqn{(x, y, v_x, v_y)} representing position and velocity.
#' @param delta Numeric scalar. Time step size.
#' @param push Numeric vector of length 2. External forcing (push) on the system.
#' @param potential_params List of parameters for the Gaussian mixture potential 
#'   (used if \code{split_around_fixed_point=TRUE}), with elements:
#'   \describe{
#'     \item{\code{alpha}}{Numeric vector of mixture weights.}
#'     \item{\code{B}}{List of precision matrices.}
#'     \item{\code{x_star}}{Matrix of fixed point coordinates (centers).}
#'   }
#' @param ind_fixed_point Index of fixed point of the mixture of gaussian
#' potential
#'
#' @return Numeric vector of length 4. Updated state after ODE step.
solve_ODE<-function(U,delta,push,
                    potential_params=NULL,
                    ind_fixed_point=NULL) {
  
  #get position
  X<-U[1:2]
  
  alpha<-potential_params$alpha;B<-potential_params$B
  x_star<-potential_params$x_star
  
  if (!(is.null(ind_fixed_point))) {
    
    
    
    #choose center
    l<-ind_fixed_point
    
    #extract parameters for the chosen center
    B_l<-B[[l]];alpha_l<-alpha[l];x_star_l<-x_star[l,]
    
    #compute non linear term for ODE solution
    #compute mahalanobis distances
    quad<-t(X-x_star_l)%*%B_l%*%(X-x_star_l)
    e_l<-as.numeric(exp(-quad))
    gv<-push+mix_gaussian_grad_cpp(X,x_star,
                               list(B=B,alpha=alpha),
                               exclude=l)+
      2*alpha_l*(e_l-1)*B_l%*%(X-x_star_l)
    
    #ODE solution
    U_hat<-U-delta*c(0,0,gv)
  }
  
  else {
    
    #ODE solution
    potential_grad<-mix_gaussian_grad_cpp(X,x_star,list(B=B,alpha=alpha),
                                          exclude=integer(0))
    U_hat<-U-delta*c(0,0,push)-delta*c(0,0,potential_grad)
    
  }
  
  return(U_hat)
}

#' Solve stochastic SDE part of the Langevin splitting scheme
#'
#' Computes the SDE flow for one time step of the splitting scheme, 
#' either naively or around a fixed point of a Gaussian mixture potential.
#'
#' @param U Numeric vector of length 4. Current state vector, typically 
#'   \eqn{(x, y, v_x, v_y)} representing position and velocity.
#' @param delta Numeric scalar. Time step size.
#' @param tau Numeric scalar. Correlation timescale parameter.
#' @param nu Numeric scalar. Velocity scale parameter.
#' @param omega Numeric scalar. Rotation parameter.
#' @param potential_params List of parameters for the Gaussian mixture potential 
#'   (used if \code{split_around_fixed_point=TRUE}), with elements:
#'   \describe{
#'     \item{\code{alpha}}{Numeric vector of mixture weights.}
#'     \item{\code{B}}{List of precision matrices.}
#'     \item{\code{x_star}}{Matrix of fixed point coordinates (centers).}
#'   }
#' @param ind_fixed_point Index of fixed point of the mixture of gaussian
#' potential
#' @param L Optional. Precomputed link matrix for the SDE step.
#' @param Q Optional. Precomputed covariance matrix for the SDE step.
#' @return A list with components:
#' \describe{
#'   \item{\code{mean}}{Numeric vector of length 4. Mean of the propagated state.}
#'   \item{\code{Q}}{Covariance matrix of the propagated state.}
#' }
solve_SDE<-function(U,delta,tau,nu,omega,potential_params=NULL,
                    ind_fixed_point=NULL,L=NULL,Q=NULL) {
  
  X<-U[1:2]
  
  # If splitting around fixed point, compute SDE solution accordingly
  if (!(is.null(ind_fixed_point))) {
    
    # If L and Q are already provided, use them directly
    if (!(is.null(L)) & !is.null(Q)) {
      
      x_star<-potential_params$x_star
      l<-ind_fixed_point
      x_star_l<-x_star[l,]
      
      #Mean vector
      mean<-L%*%(U-c(x_star_l,0,0))+c(x_star_l,0,0)
      
    }
    
    #else compute L and Q based on splitting around fixed point
    else {
  
    alpha<-potential_params$alpha;B<-potential_params$B
    x_star<-potential_params$x_star
    process_cov <- rbind(
      cbind(matrix(0,2,2), matrix(0,2,2)),
      cbind(matrix(0,2,2), 4*nu^2/pi/tau*diag(2)))
    
    #choose center
    l<-ind_fixed_point
    
    #extract parameters for the chosen center
    B_l<-B[[l]];alpha_l<-alpha[l];x_star_l<-x_star[l,]
    
    A<-rbind(cbind(matrix(0,ncol=2,nrow=2),diag(2)),
             cbind(-2*alpha_l*B_l,-matrix(c(1/tau,-omega,omega,1/tau),
                                          nrow=2,byrow=TRUE)))
    if (is.null(L)) {
      L<-expm(A*delta)
    }
    if (is.null(Q)) {
      Q<-OU_cov_exact(A, process_cov,delta,L)
    }
    
    mean<-L%*%(U-c(x_star_l,0,0))+c(x_star_l,0,0)
    }

  }
  
  # Else, compute SDE solution naively
  else {
    
    #Link matrix
    if (is.null(L)) {
      L<-RACVM_link(tau,omega,delta)
    }
    
    # Cov matrix
    if (is.null(Q)) {
      Q<-RACVM_cov(tau,nu,omega,delta)
    }
    #Mean vector
    mean<-L%*%U
  }
  
  return(list(mean=mean,Q=Q))
}



#' Simulate penalised Langevin SDE trajectory in two dimensions.
#'
#' @param n_sim number of trajectories to simulate
#' @param sde_params List of movement parameters tau,nu,omega
#' @param potential_params List of parameters for the potential.
#' Only necessary when fixed_point_ind is not null, and only works with 
#' mixture of Gaussian potentials. The list has elements 
#' \code{potential_params$alpha} which is a vector of weights, 
#' \code{potential_params$B} which is a list of covariance matrices, 
#' and \code{potential_params$x_star} which is a matrix of coordinates of each 
#' center in the mixture.
#' @param error_params List of parameters of the measurement error distribution
#' @param error_dist Distribution of the measurement error. 
#' Options are normal,scaled_t and argos
#' @param polygon A matrix defining the boundary polygon (n x 2 coordinates)
#' @param lambda Penalization parameter.
#' @param U0 initial state. Initial position and velocity for RACVM, 
#' and only initial position for BM and OU
#' @param N Integer, number of time steps to simulate
#' @param dt Numeric, time step size
#' @param scheme String, splitting scheme to use. Options are "Lie-Trotter" and "Strang"
#' @param split_around_fixed_point Boolean whether
#' to do the splitting adaptatively around a fixed point of the potential
#' @param ncores Number of cores to use for parallel simulation
#' @param seed Random seed for reproducibility
#'
#' @return A dataframe containing with Nxn_sim rows with columns
#' \itemize{
#'   \item X1, X2: true position
#'   \item V1, V2 : true velocity (for RACVM)
#'   \item Y1, Y2 : Noisy position
#'   \item ID : trajectory identifier
#' }
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom foreach %dopar% foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#'
#'
#' @details The function simulates movement according to a penalised Langevin SDE. When the 
#' animal hits the polygon boundary, a repulsive
#' force is applied. Observations include measurement error.
#'
simulate_2D_trajectory <- function(n_sim,sde_params,
                                   potential_params=NULL,error_params,
                                   error_dist,polygon,lambda,
                                   U0, N, dt,scheme="Lie-Trotter",
                                   split_around_fixed_point=TRUE,
                                   ncores=parallel::detectCores()-1,
                                   seed=42) {
  

  tau<-sde_params$tau;nu <- sde_params$nu
  omega<-sde_params$omega
  
  alpha<-potential_params$alpha;B<-potential_params$B
  x_star<-potential_params$x_star
  
  
  cl <- makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  clusterEvalQ(cl, {
    library(here)
    source(here("src/packages.R"))
    source(here("src/error_models.R"))
    source(here("src/simulate.R"))
    source(here("src/utility.R"))
    Rcpp::sourceCpp(here("src/utility.cpp"))
    
  })
  
  all_trajectories <- foreach(j = 1:n_sim, .combine=rbind) %dopar% {
    
    set.seed(seed*j)
    
    U <- matrix(0, nrow = N, ncol = length(U0))
    Y <- matrix(0, nrow = N, ncol = 2)
    U[1, ] <- U0
    Y[1, ] <- U0[1:2]
    
    for (i in 1:(N - 1)) {
      
      push<-compute_push(U[i,1:2],polygon,lambda)
      
      ind_fixed_point<-if (split_around_fixed_point)
        choose_center_cpp(U[i,1:2], x_star, 
                      list(B=B,alpha=alpha)) 
      else NULL
        
      if (scheme=="Lie-Trotter") {
        U_hat<-solve_ODE(U[i,],dt,push,potential_params,
                                   ind_fixed_point)
        
        #SDE mean and covariance
        OU_solution<-solve_SDE(U_hat,dt,tau,nu,omega,potential_params,
                               ind_fixed_point)
        Q<-OU_solution$Q
        mean<-OU_solution$mean
          
        U[i+1,]<- mvrnorm(1,mean,Q)
        
      } else if (scheme =="Strang") {
          
          U_hat<-solve_ODE(U[i,],dt/2,push,potential_params,
                           ind_fixed_point)
          
          #SDE mean and covariance
          OU_solution<-solve_SDE(U_hat,dt,tau,nu,omega,potential_params,
                                 ind_fixed_point)
          Q<-OU_solution$Q
          mean<-OU_solution$mean
          
          Z<- mvrnorm(1,mean,Q)
          
          new_push<-compute_push(Z[1:2],polygon,lambda)
          
          if (split_around_fixed_point) {
            
            #compute mahalanobis distances
            quad <- sapply(seq_along(alpha), function(idx) {
              t(Z[1:2] - x_star[idx, ]) %*% B[[idx]] %*% (Z[1:2] - x_star[idx, ])
            })
            
            l <- ind_fixed_point
            B_l <- B[[l]]; alpha_l <- alpha[l]; x_star_l <- x_star[l,]
            e_l <- as.numeric(exp(-quad[l]))
            
            grad_Z <- mix_gaussian_grad_cpp(Z[1:2], x_star, list(B=B, alpha=alpha),
                                             exclude = l) +
              2 * alpha_l * (e_l - 1) * B_l %*% (Z[1:2] - x_star_l)
            
            U[i+1,]<-Z-dt/2*c(0,0,new_push)-dt/2*c(0,0,grad_Z)
            
          } else {
            
            potential_grad <- mix_gaussian_grad_cpp(Z[1:2], x_star, list(B=B, alpha=alpha),
                                                  exclude = integer(0))
            U[i+1,]<-Z-dt/2*c(0,0,new_push)-dt/2*c(0,0,potential_grad)
          }
        }
      epsilon <- switch(error_dist,
                        normal   = mvrnorm(1, mu=rep(0,2), Sigma=error_params$sigma_obs^2*diag(2)),
                        scaled_t = c(rscaledt(1, error_params$scale, error_params$df),
                                     rscaledt(1, error_params$scale, error_params$df)),
                        argos    = rmvt_mixture(1, error_params),
                        stop("Unknown error_dist"))
      
      Y[i + 1, ] <- U[i + 1, 1:2] + epsilon
    }
    
    df<-data.frame(U = U, Y = Y, time=seq(0,(N-1)*dt,by=dt),ID=as.factor(rep(j,N)))
    colnames(df)[1:6]<-c("X1","X2","V1","V2","Y1","Y2")
    df
    }
  
  stopCluster(cl)
  return(all_trajectories)
}

#' Visualizes both true and observed animal movement trajectories 
#' within a specified polygon boundary.
#' @param data df from `simulate_2D_trajectory()`. If multiple IDs, individual trajectories are on individual facets
#' @param polygon SpatialPolygons object defining the boundary (from `sp` package)
#' @param true_opacity Opacity for true trajectory line (0-1, default: 0.5)
#' @param obs_opacity Opacity for observed trajectory line (0-1, default: 0.5)
#' @param potential_opacity Opacity for potential grid (0-1, default: 0.4)
#' @param true_size Size of points/lines for true trajectory (default: 0.2)
#' @param obs_size Size of points/lines for observed trajectory (default: 0.1)
#' @param show_obs Logical indicating whether to show observed trajectory (default: TRUE)
#' @param potential_grid spatial grid with columns x, y and H for the underlying 
#' potential driving the movement
#' @param n_sub Integer indicating subsampling rate for plotting 
#' (default: 1, no subsampling)
#' @return A ggplot object showing:
#' \itemize{
#'   \item Polygon boundary (light blue fill)
#'   \item True trajectory (black line and points)
#'   \item Observed trajectory (red line and points with measurement error)
#' }
#' @export
#' @details The plot helps visualize the relationship between true movement and noisy observations.
#'

plot_trajectory_in_polygon <- function(data, polygon, 
                                       true_opacity = 0.2,obs_opacity=0.2,
                                       potential_opacity=0.4,
                                       true_size=0.2,obs_size=0.1,
                                       show_obs=TRUE,
                                       potential_grid=NULL,n_sub=1) {
  #subsample data
  if ("ID" %in% names(data) & length(unique(data$ID))>1) {
    data_sub <- do.call(rbind, lapply(split(data, data$ID), function(df) {
      df[seq(1, nrow(df), by = n_sub), , drop = FALSE]
    }))
    rownames(data_sub) <- NULL
  }
  else {
    data_sub<-data[seq(1, nrow(data), by = n_sub),]
  }
  
  p<-ggplot()
  
  polygon_coords <- as.data.frame(polygon@coords)
  
  p<- p +
    geom_polygon(data = polygon_coords, aes(x = V1, y = V2), 
                 fill = NA, color = "grey40",size=1) +
    geom_point(data = data_sub, aes(x = X1, y = X2), size = true_size,alpha=true_opacity) +
    geom_path(data = data_sub, aes(x = X1, y = X2),size = true_size, alpha = true_opacity)+
    xlab(expression(X[1]))+ylab(expression(X[2]))+ theme_minimal(base_size=18)
  
  if (!(is.null(potential_grid))) {
    
    p<-p+geom_raster(data=potential_grid,aes(x,y,fill = H), interpolate = TRUE,alpha=potential_opacity) +
      geom_contour(data=potential_grid,aes(x,y,z=H),color = "white", alpha = potential_opacity) +
      scale_fill_gradient(low = "grey10", high = "white")+
      coord_equal()
  }
  
  if (show_obs) {
    p <- p +
      geom_point(data = data_sub, aes(x = Y1, y = Y2), 
                 color = "red", size = obs_size, alpha = obs_opacity) +
      geom_segment(data = data_sub, 
                   aes(x = X1, y = X2, xend = Y1, yend = Y2),
                   color = "red", alpha = obs_opacity)
  }
  
    
  if ("ID" %in% names(data) & length(unique(data$ID))>1) {
    p <- p + facet_wrap(~ID)
  }
  

  print(p)
  return (p)
}
