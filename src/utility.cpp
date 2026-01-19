#include <RcppArmadillo.h> 
using namespace Rcpp;

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
arma::vec solve_ODE_cpp(const arma::vec& U,
                        double delta,
                        const arma::vec& push,
                        const List& potential_params,
                        Nullable<int> ind_fixed_point = R_NilValue) {
  
  // --- extract position ---
  arma::vec X = U.subvec(0,1); // first 2 elements

  // --- extract potential parameters ---
  arma::vec alpha = potential_params["alpha"];
  List B_list = potential_params["B"];
  arma::mat x_star = potential_params["x_star"];

  arma::vec U_hat = U; 
  
  if (ind_fixed_point.isNotNull()) {
    
    int l = as<int>(ind_fixed_point) - 1; // convert 1-based R index to 0-based C++
    
    // extract B_l and alpha_l
    arma::mat B_l = B_list[l];
    double alpha_l = alpha[l];
    arma::vec x_star_l = x_star.row(l).t();

    // mahalanobis distance
    arma::vec diff = X - x_star_l;
    double quad = arma::dot(diff, B_l * diff);
    double e_l = std::exp(-quad);

    // gradient from mix_gaussian_grad_cpp excluding l
    IntegerVector exclude = IntegerVector::create(l+1); // R-style index for mix_gaussian_grad_cpp
    arma::vec grad = mix_gaussian_grad_cpp(X, x_star, potential_params, exclude);

    // non-linear term
    arma::vec gv = push + grad + 2.0 * alpha_l * (e_l - 1.0) * (B_l * diff);
    
    // Build a vector of zeros with same length as U
    arma::vec tmp(U.n_elem, arma::fill::zeros);

    // Assign push + potential_grad to positions 3 and 4 (0-based indices 2 and 3)
    tmp.subvec(2, 3) = gv;
    U_hat = U - delta * tmp;

  } else {

    // no fixed point
    IntegerVector exclude; // empty
    arma::vec potential_grad = mix_gaussian_grad_cpp(X, x_star, potential_params, exclude);
    
    // Build a vector of zeros with same length as U
    arma::vec tmp(U.n_elem, arma::fill::zeros);

    // Assign push + potential_grad to positions 3 and 4 (0-based indices 2 and 3)
    tmp.subvec(2, 3) = push + potential_grad;
    U_hat = U - delta * tmp;
  }

  return U_hat;
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



// [[Rcpp::export]] 
List product_gaussian_cpp(const arma::mat& P1,
	const arma::mat& P2,
  	const arma::vec& mean1,
   	const arma::vec& mean2, 
   	const arma::mat& M) { 
   
// Sigma_inv = P1 + t(M) * P2 * M 
arma::mat Sigma_inv = P1 + M.t() * P2 * M; 

// Force symmetry 
Sigma_inv = 0.5 * (Sigma_inv + Sigma_inv.t()); 

// b = P1 * mean1 + t(M) * P2 * mean2
arma::vec b = P1 * mean1 + M.t() * P2 * mean2; 

// Cholesky decomposition (upper triangular) 
arma::mat R = arma::chol(Sigma_inv); 

// Solve Sigma_inv * m = b using Cholesky 
arma::vec z = arma::solve(arma::trimatl(R.t()), b);
arma::vec m = arma::solve(arma::trimatu(R), z); 

// Compute Cholesky of covariance: chol(Sigma) 
arma::mat chol_Sigma = arma::solve(arma::trimatl(R.t()), arma::eye(R.n_rows, R.n_cols));

return List::create( Named("mean") = m, Named("chol") = chol_Sigma );
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





