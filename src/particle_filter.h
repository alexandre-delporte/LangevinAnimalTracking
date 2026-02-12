#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include <RcppArmadillo.h>
#include <string>

// ------------------- Particle Propagation -------------------
// Propagate a particle using Langevin dynamics
arma::vec propagate_langevin_particle_cpp(const arma::vec& U,
                                           const arma::vec& y,
                                           const arma::mat& M,
                                           double delta,
                                           const arma::vec& push,
                                           const Rcpp::List& potential_params,
                                           double tau,
                                           double nu,
                                           double omega,
                                           double lambda,
                                           const std::string& error_dist,
                                           const Rcpp::List& error_params,
                                           const std::string& scheme = "Lie-Trotter",
                                           const arma::mat& polygon_coords = arma::mat(),
                                           Rcpp::Nullable<int> ind_fixed_point = R_NilValue,
                                           Rcpp::Nullable<arma::mat> L = R_NilValue,
                                           Rcpp::Nullable<arma::mat> Q = R_NilValue,
                                           double proposal_weight = 0.5,
                                           bool verbose = false);

// Compute the weight of a particle given its prediction
double compute_langevin_weight_cpp(const arma::vec& U_pred,
                                   const arma::vec& U_prev,
                                   const arma::vec& y,
                                   const arma::mat& M,
                                   double delta,
                                   const arma::vec& push,
                                   const Rcpp::List& potential_params,
                                   double tau,
                                   double nu,
                                   double omega,
                                   const std::string& error_dist,
                                   const Rcpp::List& error_params,
                                   const std::string& scheme = "Lie-Trotter",
                                   Rcpp::Nullable<int> ind_fixed_point = R_NilValue,
                                   Rcpp::Nullable<arma::mat> L = R_NilValue,
                                   Rcpp::Nullable<arma::mat> Q = R_NilValue,
                                   double proposal_weight = 0.5,
                                   bool verbose = false);
                                   
// Full particle filter implementation
Rcpp::List particle_filter2D_cpp(const arma::mat& observations,const Rcpp::List& sde_params,
    const Rcpp::List& potential_params,
    const Rcpp::List& error_params,
    const std::string& error_dist,
    const arma::mat& polygon_coords,
    const arma::vec& U0,
    double lambda,
    int num_particles,
    const std::string& scheme,
    bool split_around_fixed_point,
    double ESS_threshold,
    double proposal_weight,
    bool verbose = false);

#endif // PARTICLE_FILTER_H

