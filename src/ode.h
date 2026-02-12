#ifndef ODE_H
#define ODE_H

#include <RcppArmadillo.h>
#include "utility.h"

// Solve ODE step for Langevin splitting scheme
arma::vec solve_ODE_cpp(const arma::vec& U,
                        double delta,
                        const arma::vec& push,
                        const Rcpp::List& potential_params,
                        Rcpp::Nullable<int> ind_fixed_point = R_NilValue);

#endif // ODE_H

