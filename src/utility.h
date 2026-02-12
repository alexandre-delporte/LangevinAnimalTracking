#ifndef UTILITY_H
#define UTILITY_H

#include <RcppArmadillo.h>
#include <unordered_set>

// Note: Do NOT use `using namespace` in header files; instead, fully qualify types

// Gradient of Gaussian mixture excluding some components
arma::vec mix_gaussian_grad_cpp(const arma::vec& x,
                                const arma::mat& x_star,
                                const Rcpp::List& params,
                                const Rcpp::IntegerVector& exclude);

// Log-density of multivariate normal given Cholesky factor
double log_dmvnorm_chol_cpp(const arma::vec& x,
                            const arma::vec& mean,
                            const arma::mat& cholSigma);

// Product of two Gaussian distributions
Rcpp::List product_gaussian_cpp(const arma::mat& P1,
                                const arma::mat& P2,
                                const arma::vec& mean1,
                                const arma::vec& mean2,
                                const arma::mat& M,
                                double proposal_weight);

// Cholesky decomposition and inverse
arma::mat chol_cpp(const arma::mat& Q);
arma::mat chol2inv_cpp(const arma::mat& R);

// Choose center of mixture based on R-style index
int choose_center_cpp(const arma::vec& x,
                      const arma::mat& x_star,
                      const Rcpp::List& params);
                      
// Vectorized version: choose center for multiple positions
arma::ivec choose_center_matrix_cpp(const arma::mat& X,
                                    const arma::mat& x_star,
                                    const Rcpp::List& params);

// Polygon utilities
Rcpp::List closest_point_on_boundary_cpp(const arma::vec& x,
                                         const arma::mat& coords,
                                         bool grad = true);
bool is_point_inside_polygon_cpp(const arma::vec& x,
                                 const arma::mat& coords);

// Compute "push" vector away from polygon boundary
arma::vec compute_push_cpp(const arma::vec& x,
                           const arma::mat& coords,
                           double lambda);
                           
// Vectorized version: compute push for multiple positions
arma::mat compute_push_matrix_cpp(const arma::mat& X,
                                   const arma::mat& coords,
                                   double lambda);

// Scaled t density
double dscaledt_cpp(double y,
                    double mean,
                    double scale,
                    double df,
                    bool logp = false);

// Mixture of multivariate t densities
double dmvt_mixture_cpp(const arma::vec& x,
                        const arma::vec& mean,
                        const Rcpp::List& params,
                        bool logp = false);

#endif // UTILITY_H

