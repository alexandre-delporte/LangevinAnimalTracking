

#include "ode.h"
#include "utility.h" 
#include <RcppArmadillo.h> 
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
arma::vec solve_ODE_cpp(const arma::vec& U,
                        double delta,
                        const arma::vec& push,
                        const List& potential_params,
                        Nullable<int> ind_fixed_point) {
  
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


