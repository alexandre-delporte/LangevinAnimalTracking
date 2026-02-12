#ifndef SDE_H
#define SDE_H

#include <RcppArmadillo.h>

// Compute exact covariance matrix for OU process with given exp(A*h)
arma::mat OU_cov_exact_cpp(const arma::mat& A, 
                           const arma::mat& Gamma, 
                           double h,
                           const arma::mat& expAh);

// Compute exact covariance matrix for OU process, compute exp(A*h) internally
arma::mat OU_cov_exact_cpp_full(const arma::mat& A, 
                                const arma::mat& Gamma, 
                                double h);

// RACVM covariance matrix (4x4)
arma::mat RACVM_cov_cpp(double tau, double nu, double omega, double dt);

// RACVM link matrix (state transition, 4x4)
arma::mat RACVM_link_cpp(double tau, double omega, double dt);

// Solve SDE step for Langevin splitting scheme
Rcpp::List solve_SDE_cpp(const arma::vec& U,
                         double delta,
                         double tau,
                         double nu,
                         double omega,
                         const Rcpp::List& potential_params,
                         Rcpp::Nullable<int> ind_fixed_point = R_NilValue,
                         Rcpp::Nullable<arma::mat> L_provided = R_NilValue,
                         Rcpp::Nullable<arma::mat> Q_provided = R_NilValue);

#endif // SDE_H

