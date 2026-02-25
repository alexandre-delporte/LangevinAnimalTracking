
#include "utility.h"
#include <RcppArmadillo.h> 
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec mix_gaussian_grad_cpp(const arma::vec& x,
                                     const arma::mat& x_star,
                                     const List& params,
                                     const IntegerVector& exclude) {
    arma::vec alpha = params["alpha"];
    List B_list = params["B"];
    int m = alpha.size();
    arma::vec grad(2, arma::fill::zeros);

    std::unordered_set<int> excl;
    for(int k = 0; k < exclude.size(); ++k) excl.insert(exclude[k]-1);

    for(int j = 0; j < m; ++j){
        if(excl.count(j)) continue;
        arma::vec diff = x - arma::vec(x_star.row(j).t());
        arma::mat B = B_list[j];
        double quad = arma::dot(diff, B*diff);
        grad += 2.0 * alpha[j] * std::exp(-quad) * (B*diff);
    }
    return grad;
}



// [[Rcpp::export]]
double log_dmvnorm_chol_cpp(const arma::vec& x,
                            const arma::vec& mean,
                            const arma::mat& cholSigma) {
  
  int d = x.n_elem;
  
  arma::vec diff = x - mean;

  arma::mat v = arma::solve(
    arma::trimatl(cholSigma.t()),diff);

  double quad = arma::dot(v, v);

  // log determinant term
  double logdet = 2.0 * arma::sum(arma::log(cholSigma.diag()));

  return -0.5 * (quad + logdet + d * std::log(2.0 * M_PI));
}


//' Product of Two Gaussian Distributions
//'
//' Computes the mean and covariance of the product of two Gaussian densities:
//' \deqn{N(mean1, P1^{-1}) \times N(L mean2, P2^{-1}) \propto N(m, \Gamma).}
//'
//' @param P1 Precision matrix (inverse covariance) of dimension \eqn{n \times n}.
//' @param P2 Precision matrix (inverse covariance) of dimension \eqn{q \times q}.
//' @param mean1 Numeric vector of length \eqn{n}. Mean of the first Gaussian.
//' @param mean2 Numeric vector of length \eqn{q}. Mean of the second Gaussian.
//' @param M Link matrix of dimension \eqn{q \times n} that maps the state space.
//' @param proposal_weight Scalar between 0 and 1 to attribute weight to 
//' the first gaussian. Value 0.5 means equal weight for both gaussians.
//'
//' @return A list with components:
//' \describe{
//'   \item{\code{mean}}{Posterior mean vector.}
//'   \item{\code{cov}}{Posterior covariance matrix.}
//'   \item{\code{chol}}{Upper triangular Cholesky factor of covariance matrix,
//'     where \code{cov = t(chol) \%*\% chol} (matching R's \code{chol()} convention).}
//' }
//'
//' @export
// [[Rcpp::export]]
List product_gaussian_cpp(const arma::mat& P1,
                           const arma::mat& P2,
                           const arma::vec& mean1,
                           const arma::vec& mean2,
                           const arma::mat& M,
                           double proposal_weight) {
  
  // Compute precision matrix of the product
 arma::mat Sigma_inv =
	    2*proposal_weight*P1 + 2*(1-proposal_weight)*M.t() * P2 * M;
  
  // Force symmetry (numerical stability)
  Sigma_inv = 0.5 * (Sigma_inv + Sigma_inv.t());
  
  // Compute precision-weighted mean
  arma::vec b =
	    2*proposal_weight * P1 * mean1 + 2*(1 - proposal_weight) * M.t() * P2 * mean2;
  
  // Solve for mean: Sigma_inv * m = b
  arma::vec m = arma::solve(Sigma_inv, b);
  
  // Compute covariance: Gamma = inv(Sigma_inv)
  arma::mat Gamma = arma::inv(Sigma_inv);
  
  // Ensure symmetry of covariance
  Gamma = 0.5 * (Gamma + Gamma.t());
  
  // Compute Cholesky decomposition
  // arma::chol() returns upper triangular U where Gamma = U.t() * U
  arma::mat chol_Gamma = arma::chol(Gamma);
  
  return List::create(
    Named("mean") = m,
    Named("cov") = Gamma,
    Named("chol") = chol_Gamma
  );
}

 
// [[Rcpp::export]]
arma::mat chol_cpp(const arma::mat& Q) {
  // Compute the upper-triangular Cholesky factor
  arma::mat U = arma::chol(Q);  // Armadillo defaults to upper-triangular
  return U;
}

// [[Rcpp::export]]
arma::mat chol2inv_cpp(const arma::mat& R) {
    // R = upper-triangular Cholesky factor such that R'R = Q
    // Solve R' * R * X = I  => X = inv(Q)
    arma::mat I = arma::eye(R.n_rows, R.n_cols);
    arma::mat invQ = arma::solve(arma::trimatl(R.t()), I);
    invQ = arma::solve(arma::trimatu(R), invQ);
    return invQ;
}

// [[Rcpp::export]]
int choose_center_cpp(const arma::vec& x,
                      const arma::mat& x_star,
                      const List& params) {

  arma::vec alpha = as<arma::vec>(params["alpha"]);
  List B_list     = params["B"];

  int J = alpha.n_elem;
  arma::vec L_vals(J);

  for (int j = 0; j < J; ++j) {

    arma::vec center = x_star.row(j).t();
    arma::mat B      = as<arma::mat>(B_list[j]);

    arma::vec diff = x - center;

    // quad = t(diff) %*% B %*% diff
    double quad = arma::as_scalar(diff.t() * B * diff);

    // q = t(diff) %*% t(B) %*% (B %*% diff)
    arma::vec Bdiff = B * diff;
    double q = arma::as_scalar(diff.t() * B.t() * Bdiff);

    L_vals(j) = 2.0 * std::log(alpha(j)) - 2.0 * quad + std::log(q);
  }

  // Return R-style index
  return L_vals.index_max() + 1;
}

// [[Rcpp::export]]
arma::ivec choose_center_matrix_cpp(const arma::mat& X,
                                    const arma::mat& x_star,
                                    const List& params) {
  
  int n_particles = X.n_rows;
  arma::ivec centers(n_particles);
  
  arma::vec alpha = as<arma::vec>(params["alpha"]);
  List B_list     = params["B"];
  int J = alpha.n_elem;
  
  // Process each particle
  for (int k = 0; k < n_particles; ++k) {
    arma::vec x = X.row(k).t();  // Extract position as column vector
    arma::vec L_vals(J);
    
    // Compute L values for all centers
    for (int j = 0; j < J; ++j) {
      arma::vec center = x_star.row(j).t();
      arma::mat B      = as<arma::mat>(B_list[j]);
      
      arma::vec diff = x - center;
      
      // quad = t(diff) %*% B %*% diff
      double quad = arma::as_scalar(diff.t() * B * diff);
      
      // q = t(diff) %*% t(B) %*% (B %*% diff)
      arma::vec Bdiff = B * diff;
      double q = arma::as_scalar(diff.t() * B.t() * Bdiff);
      
      L_vals(j) = 2.0 * std::log(alpha(j)) - 2.0 * quad + std::log(q);
    }
    
    // Return R-style index (1-based)
    centers(k) = L_vals.index_max() + 1;
  }
  
  return centers;
}

// [[Rcpp::export]]
List closest_point_on_boundary_cpp(const arma::vec& x,
                                   const arma::mat& coords,
                                   bool grad) {

  int n = coords.n_rows;
  double min_dist = std::numeric_limits<double>::infinity();

  arma::vec closest_point(2, arma::fill::zeros);
  arma::mat gradient;
  bool has_gradient = false;

  for (int i = 0; i < n - 1; ++i) {

    arma::vec v = coords.row(i).t();
    arma::vec w = coords.row(i + 1).t();
    arma::vec e = w - v;

    double l2 = arma::dot(e, e);

    arma::vec projection;
    double t_val;

    if (l2 == 0.0) {
      projection = v;
      t_val = 0.0;
    } else {
      t_val = arma::dot(x - v, e) / l2;
      t_val = std::max(0.0, std::min(1.0, t_val));
      projection = v + t_val * e;
    }

    double dist = arma::norm(x - projection, 2);

    if (dist < min_dist) {
      min_dist = dist;
      closest_point = projection;

      if (grad && t_val > 0.0 && t_val < 1.0) {
        gradient = (e * e.t()) / l2;
        has_gradient = true;
      } else {
        has_gradient = false;
      }
    }
  }

  if (grad && has_gradient) {
    return List::create(
      Named("point") = closest_point,
      Named("gradient") = gradient
    );
  } else {
    return List::create(
      Named("point") = closest_point,
      Named("gradient") = R_NilValue
    );
  }
}

bool is_point_inside_polygon_cpp(const arma::vec& x,
                                 const arma::mat& coords) {

  // Access sp::point.in.polygon
  Function pip("point.in.polygon", Environment::namespace_env("sp"));

  // Call sp::point.in.polygon(x, y, polyx, polyy)
  IntegerVector status = pip(
    x(0),
    x(1),
    coords.col(0),
    coords.col(1)
  );

  // Same logic as R:
  // status == 1 (inside) or 2 (on edge)
  return (status[0] == 1 || status[0] == 2);
}

// [[Rcpp::export]]
arma::vec compute_push_cpp(const arma::vec& x,
                           const arma::mat& coords,
                           double lambda) {

  arma::vec zero(2, arma::fill::zeros);

  if (!std::isfinite(lambda)) {
    return zero;
  }

  if (is_point_inside_polygon_cpp(x, coords)) {
    return zero;
  }

  List proj = closest_point_on_boundary_cpp(x, coords, false);
  arma::vec p = proj["point"];

  return (x - p) / lambda;
}

// [[Rcpp::export]]
arma::mat compute_push_matrix_cpp(const arma::mat& X,
                                   const arma::mat& coords,
                                   double lambda) {
  
  int n_particles = X.n_rows;
  arma::mat push_matrix(n_particles, 2);
  
  // If lambda is infinite, return matrix of zeros
  if (!std::isfinite(lambda)) {
    push_matrix.zeros();
    return push_matrix;
  }
  
  // Process each particle
  for (int k = 0; k < n_particles; ++k) {
    arma::vec x = X.row(k).t();  // Extract position as column vector
    
    // Check if inside polygon
    if (is_point_inside_polygon_cpp(x, coords)) {
      push_matrix.row(k).zeros();
    } else {
      // Find closest point on boundary
      List proj = closest_point_on_boundary_cpp(x, coords, false);
      arma::vec p = proj["point"];
      
      // Compute push: (x - p) / lambda
      arma::vec push = (x - p) / lambda;
      push_matrix.row(k) = push.t();
    }
  }
  
  return push_matrix;
}


// [[Rcpp::export]]
double dscaledt_cpp(double y,
                    double mean,
                    double scale,
                    double df,
                    bool logp) {
  
  // Standardize
  double z = (y - mean) / scale;
  
  // Evaluate Student-t density using Rf_dt
  double dens = Rf_dt(z, df, logp);  // logp == TRUE gives log density
  
  if (!logp) {
    dens /= scale;  // scale correction for scaled t
  } else {
    dens -= std::log(scale);  // log version
  }
  
  return dens;
}


// [[Rcpp::export]]
double dmvt_mixture_cpp(const arma::vec& x,
                        const arma::vec& mean,
                        const List& params,
                        bool logp) {
  
  double sigma_obs = as<double>(params["sigma_obs"]);
  double rho       = as<double>(params["rho"]);
  double a         = as<double>(params["a"]);
  double df        = as<double>(params["df"]);
  double p         = as<double>(params["p"]);

  // --- Covariance matrices ---
  arma::mat Sigma1 = sigma_obs * sigma_obs * arma::mat{{1, rho*std::sqrt(a)},
                                                       {rho*std::sqrt(a), 1}};
  arma::mat Sigma2 = sigma_obs * sigma_obs * arma::mat{{1, -rho*std::sqrt(a)},
                                                       {-rho*std::sqrt(a), 1}};
  
  // --- Precompute constants ---
  int d = 2; // dimension
  
  arma::vec diff = x - mean;

  // --- Mahalanobis terms ---
  double mahal1 = arma::as_scalar(diff.t() * arma::inv(Sigma1) * diff);
  double mahal2 = arma::as_scalar(diff.t() * arma::inv(Sigma2) * diff);

  // --- Multivariate t densities ---
  double c1 = std::tgamma((df + d)/2.0) / (std::tgamma(df/2.0) * std::pow(df*M_PI, d/2.0) * std::sqrt(arma::det(Sigma1)));
  double c2 = std::tgamma((df + d)/2.0) / (std::tgamma(df/2.0) * std::pow(df*M_PI, d/2.0) * std::sqrt(arma::det(Sigma2)));

  double dens1 = c1 * std::pow(1 + mahal1/df, -(df+d)/2.0);
  double dens2 = c2 * std::pow(1 + mahal2/df, -(df+d)/2.0);

  double mixture = p * dens1 + (1 - p) * dens2;

  if (logp) {
    return std::log(mixture);
  } else {
    return mixture;
  }
}


