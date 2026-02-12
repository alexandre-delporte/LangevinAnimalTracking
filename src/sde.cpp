#include "sde.h"
#include "utility.h" 
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


//' Compute exact covariance matrix for Ornstein-Uhlenbeck process (C++ version)
//' 
//' @param A Drift matrix (d x d)
//' @param Gamma Process covariance matrix (d x d)
//' @param h Time step
//' @param expAh Matrix exponential exp(A*h) (d x d)
//' @return Covariance matrix Q (d x d)
//' @export
// [[Rcpp::export]]
arma::mat OU_cov_exact_cpp(const arma::mat& A, 
                           const arma::mat& Gamma, 
                           double h,
                           const arma::mat& expAh) {
  
  int d = A.n_rows;
  
  // Compute Kronecker sum: M = A ⊗ I + I ⊗ A
  arma::mat I_d = arma::eye<arma::mat>(d, d);
  arma::mat M = arma::kron(A, I_d) + arma::kron(I_d, A);
  
  // Vectorize Gamma
  arma::vec vecGamma = arma::vectorise(Gamma);
  
  // Solve M * vecC = vecGamma
  arma::vec vecC = arma::solve(M, vecGamma);
  
  // Reshape vecC back to matrix C
  arma::mat C = arma::reshape(vecC, d, d);
  
  // Compute Q = exp(A*h) * C * exp(A*h)^T - C
  arma::mat Q = expAh * C * expAh.t() - C;
  
  return Q;
  }

//' Compute exact covariance matrix with matrix exponential computed internally
//' 
//' @param A Drift matrix (d x d)
//' @param Gamma Process covariance matrix (d x d)
//' @param h Time step
//' @return Covariance matrix Q (d x d)
//' @export
// [[Rcpp::export]]
arma::mat OU_cov_exact_cpp_full(const arma::mat& A, 
                                const arma::mat& Gamma, 
                                double h) {
  
  int d = A.n_rows;
  
  // Compute matrix exponential exp(A*h)
  arma::mat expAh = arma::expmat(A * h);
  
  // Compute Kronecker sum: M = A ⊗ I + I ⊗ A
  arma::mat I_d = arma::eye<arma::mat>(d, d);
  arma::mat M = arma::kron(A, I_d) + arma::kron(I_d, A);
  
  // Vectorize Gamma
  arma::vec vecGamma = arma::vectorise(Gamma);
  
  // Solve M * vecC = vecGamma
  arma::vec vecC = arma::solve(M, vecGamma);
  
  // Reshape vecC back to matrix C
  arma::mat C = arma::reshape(vecC, d, d);
  
  // Compute Q = exp(A*h) * C * exp(A*h)^T - C
  arma::mat Q = expAh * C * expAh.t() - C;
  
  return Q;
}


//' Compute the covariance matrix for the RACVM model (C++ version)
//' 
//' Computes the exact covariance matrix for the Rotational Advective 
//' Correlated Velocity Model (RACVM) used in animal movement modeling.
//' 
//' @param tau Time scale parameter (persistence)
//' @param nu Scale parameter (velocity magnitude)
//' @param omega Angular velocity parameter (rotation)
//' @param dt Time step
//' @return Covariance matrix Q (4x4 matrix)
//' @export
// [[Rcpp::export]]
arma::mat RACVM_cov_cpp(double tau, double nu, double omega, double dt) {
  
  // Derived parameters
  double beta = 1.0 / tau;
  double sigma = 2.0 * nu / std::sqrt(M_PI * tau);
  
  // Commonly used expressions
  double C = beta * beta + omega * omega;
  double sigma2 = sigma * sigma;
  double exp_beta_dt = std::exp(-beta * dt);
  double exp_2beta_dt = std::exp(-2.0 * beta * dt);
  double cos_omega_dt = std::cos(omega * dt);
  double sin_omega_dt = std::sin(omega * dt);
  
  // Variance components
  double var_xi = (sigma2 / C) * (
    dt + 
    (omega * omega - 3.0 * beta * beta) / (2.0 * beta * C) - 
    exp_2beta_dt / (2.0 * beta) + 
    2.0 * exp_beta_dt * (beta * cos_omega_dt - omega * sin_omega_dt) / C
  );
  
  double var_zeta = (sigma2 / (2.0 * beta)) * (1.0 - exp_2beta_dt);
  
  double cov1 = (sigma2 / (2.0 * C)) * (
    1.0 + exp_2beta_dt - 2.0 * exp_beta_dt * cos_omega_dt
  );
  
  double cov2 = (sigma2 / C) * (
    exp_beta_dt * sin_omega_dt - 
    (omega / (2.0 * beta)) * (1.0 - exp_2beta_dt)
  );
  
  // Construct the 4x4 covariance matrix
  arma::mat Q(4, 4);
  
  // Fill in the matrix (row by row)
  Q(0, 0) = var_xi;
  Q(0, 1) = 0.0;
  Q(0, 2) = cov1;
  Q(0, 3) = cov2;
  
  Q(1, 0) = 0.0;
  Q(1, 1) = var_xi;
  Q(1, 2) = -cov2;
  Q(1, 3) = cov1;
  
  Q(2, 0) = cov1;
  Q(2, 1) = -cov2;
  Q(2, 2) = var_zeta;
  Q(2, 3) = 0.0;
  
  Q(3, 0) = cov2;
  Q(3, 1) = cov1;
  Q(3, 2) = 0.0;
  Q(3, 3) = var_zeta;
  
  return Q;
}

//' Compute RACVM link matrix (C++ version)
//' 
//' Computes the state transition matrix for the RACVM model.
//' 
//' @param tau Time scale parameter
//' @param omega Angular velocity parameter
//' @param dt Time step
//' @return Link matrix L (4x4 matrix)
//' @export
// [[Rcpp::export]]
arma::mat RACVM_link_cpp(double tau, double omega, double dt) {
  
  double beta = 1.0 / tau;
  
  // Commonly used values
  double C = beta * beta + omega * omega;
  double exp_beta_dt = std::exp(-beta * dt);
  double cos_omega_dt = std::cos(omega * dt);
  double sin_omega_dt = std::sin(omega * dt);
  
  // Rotation matrix R
  arma::mat R(2, 2);
  R(0, 0) = cos_omega_dt;
  R(0, 1) = sin_omega_dt;
  R(1, 0) = -sin_omega_dt;
  R(1, 1) = cos_omega_dt;
  
  // Matrix A
  arma::mat A(2, 2);
  A(0, 0) = beta;
  A(0, 1) = -omega;
  A(1, 0) = omega;
  A(1, 1) = beta;
  
  // Inverse of A
  arma::mat invA = (1.0 / C) * arma::mat({{beta, omega}, {-omega, beta}});
  
  // exp(-A*dt)
  arma::mat expAdt = exp_beta_dt * R;
  
  // Identity matrix
  arma::mat I2 = arma::eye<arma::mat>(2, 2);
  
  // Construct the 4x4 link matrix
  arma::mat L = arma::zeros<arma::mat>(4, 4);
  
  // Top-left block: I_2
  L.submat(0, 0, 1, 1) = I2;
  
  // Top-right block: invA * (I_2 - expAdt)
  L.submat(0, 2, 1, 3) = invA * (I2 - expAdt);
  
  // Bottom-right block: expAdt
  L.submat(2, 2, 3, 3) = expAdt;
  
  return L;
}


//'
//' Computes the SDE flow for one time step of the splitting scheme, 
//' either naively or around a fixed point of a Gaussian mixture potential.
//'
//' @param U Numeric vector of length 4 (x1, x2, v1, v2)
//' @param delta Time step size
//' @param tau Correlation timescale parameter
//' @param nu Velocity scale parameter
//' @param omega Rotation parameter
//' @param potential_params List with elements a:
//' alpha Numeric vector of mixture weights (empty if not using fixed point)
//' B_list List of precision matrices (empty if not using fixed point)
//' x_star Matrix of fixed point coordinates (empty if not using fixed point)
//' @param ind_fixed_point Index of fixed point (1-based)
//' @param L_provided Optional pre-computed link matrix (4x4, empty matrix if not provided)
//' @param Q_provided Optional pre-computed covariance matrix (4x4, empty matrix if not provided)
//' @return List with 'mean' (length 4) and 'Q' (4x4 matrix)
//' @export
// [[Rcpp::export]]
List solve_SDE_cpp(const arma::vec& U,
                   double delta,
                   double tau,
                   double nu,
                   double omega,
                   const List& potential_params,
                   Nullable<int> ind_fixed_point,
                   Nullable<arma::mat> L_provided,
		   Nullable<arma::mat> Q_provided)
	 {
  
  arma::vec X = U.subvec(0, 1);  // Position
  arma::mat L;
  arma::mat Q;
  arma::vec mean(4);
  
   if (ind_fixed_point.isNotNull()) {
    // Splitting around fixed point

    // Extract parameters for chosen center (convert to 0-based indexing)
    arma::vec alpha = potential_params["alpha"];
    List B_list = potential_params["B"];
    arma::mat x_star = potential_params["x_star"];
    int l = as<int>(ind_fixed_point) - 1;
    arma::mat B_l = B_list[l];
    double alpha_l = alpha[l];
    arma::vec x_star_l = x_star.row(l).t();
    
    // Check if L and Q are provided
    bool L_is_provided = L_provided.isNotNull();
    bool Q_is_provided = Q_provided.isNotNull();

    if (L_is_provided && Q_is_provided) {
      // Use provided matrices
      L = Rcpp::as<arma::mat>(L_provided);
      Q = Rcpp::as<arma::mat>(Q_provided);
    } else {
      // Compute L and Q
      
      // Construct process covariance matrix Gamma
      arma::mat Gamma = arma::zeros<arma::mat>(4, 4);
      double var_process = 4.0 * nu * nu / (M_PI * tau);
      Gamma(2, 2) = var_process;
      Gamma(3, 3) = var_process;
      
      // Construct drift matrix A
      arma::mat A = arma::zeros<arma::mat>(4, 4);
      A.submat(0, 2, 1, 3) = arma::eye<arma::mat>(2, 2);  // Top-right: I_2
      A.submat(2, 0, 3, 1) = -2.0 * alpha_l * B_l;        // Bottom-left: -2*alpha*B
      
      // Bottom-right: velocity damping matrix
      arma::mat damping(2, 2);
      damping(0, 0) = 1.0 / tau;
      damping(0, 1) = -omega;
      damping(1, 0) = omega;
      damping(1, 1) = 1.0 / tau;
      A.submat(2, 2, 3, 3) = -damping;
      
      // Compute matrix exponential
      if (!L_is_provided) {
        L = arma::expmat(A * delta);
      } else {
        L = Rcpp::as<arma::mat>(L_provided);
      }
      
      // Compute covariance matrix
      if (!Q_is_provided) {
        Q = OU_cov_exact_cpp(A, Gamma, delta, L);
      } else {
        Q = Rcpp::as<arma::mat>(Q_provided);
      }
    }
    
    // Compute mean: L * (U - [x_star_l, 0, 0]) + [x_star_l, 0, 0]
    arma::vec center(4);
    center(0) = x_star_l(0);
    center(1) = x_star_l(1);
    center(2) = 0.0;
    center(3) = 0.0;
    
    mean = L * (U - center) + center;
    
  } else {
    // Naive approach (no fixed point)
    
    // Check if L and Q are provided
    bool L_is_provided = L_provided.isNotNull();
    bool Q_is_provided = Q_provided.isNotNull();

    
    // Compute or use provided link matrix
    if (!L_is_provided) {
      L = RACVM_link_cpp(tau, omega, delta);
    } else {
      L = Rcpp::as<arma::mat>(L_provided);;
    }
    
    // Compute or use provided covariance matrix
    if (!Q_is_provided) {
      Q = RACVM_cov_cpp(tau, nu, omega, delta);
    } else {
      Q = Rcpp::as<arma::mat>(Q_provided);;
    }
    
    // Compute mean: L * U
    mean = L * U;
  }
  
  return List::create(
    Named("mean") = mean,
    Named("Q") = Q,
    Named("L") = L  
  );
}
