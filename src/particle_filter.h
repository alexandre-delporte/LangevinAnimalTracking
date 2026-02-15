#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include <RcppArmadillo.h>
#include <string>

// ------------------- Particle Propagation -------------------
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
    bool verbose = false,
    bool print_timing = false);

#endif // PARTICLE_FILTER_H

