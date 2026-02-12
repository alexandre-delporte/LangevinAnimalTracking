


# Define the polygon 
polygon_coords <- matrix(c(
  0, 0,   
  5, 0,    
  5, 1,   
  4.5, 1,   
  4, 1.5,   
  4, 2.5,   
  4.5, 3,   
  5, 3,     
  5, 4,    
  3.5, 4,   
  3, 1.5,   
  2, 1.5,   
  1.5, 4,   
  1, 4,    
  0, 3.5,   
  0, 2.5,   
  1, 2, 
  1.5, 1,
  0.5, 1,     
  0, 0      # Close the polygon
), ncol = 2, byrow = TRUE)

polygon <- Polygon(10*polygon_coords)


# Define dynamics parameters
sde_params <- list(tau = 1, nu=5,omega = 0.1)


# Define initial position
U0 <- c(35, 15,0,0) 


# Potential parameters
x_star <- matrix(c(25,5,
                   35,15), ncol = 2, byrow = TRUE)

params <- list(
  alpha = c(70, 50),  
  B = list(
    matrix(c(1/3^2, 1/(5*8), 1/(5*8), 1/2^2), ncol = 2),   
    matrix(c(1/6^2, -1/(5*20), -1/(5*20), 1/10^2), ncol = 2) 
  )
)

potential_params<-list(x_star=x_star,B=params$B,alpha=params$alpha)

grid_x <- seq(0, 50, length.out = 100)
grid_y <- seq(0, 40, length.out = 100)
grid <- expand.grid(x = grid_x, y = grid_y)

# Evaluate potential on grid
grid$H <- apply(grid, 1, function(row) {
  mix_gaussian_potential(c(row[1], row[2]), x_star, params)
})


potential_grad<-function(x,x_star=matrix(c(25,5,
                                           35,15), ncol = 2, byrow = TRUE),
                         params =list(
                           alpha = c(70,50),  
                           B = list(
                             matrix(c(1/3^2, 1/(5*8), 1/(5*8), 1/2^2), ncol = 2),   
                             matrix(c(1/6^2, -1/(5*20),-1/(5*20), 1/10^2), ncol = 2) 
                           )
                         )) {
  return (mix_gaussian_grad_cpp(x,x_star,params,integer(0)))
}

potential_hessian<-function(x,x_star=matrix(c(25,5,
                                              35,15), ncol = 2, byrow = TRUE),
                            params =list(
                              alpha = c(70,50),  
                              B = list(
                                matrix(c(1/3^2, 1/(5*8), 1/(5*8), 1/2^2), ncol = 2),   
                                matrix(c(1/6^2, -1/(5*20),-1/(5*20), 1/10^2), ncol = 2) 
                              )
                            )) {
  return (mix_gaussian_hessian(x,x_star,params))
}


