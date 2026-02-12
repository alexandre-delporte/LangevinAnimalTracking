#include "particle_filter.h"
#include "sde.h"
#include "ode.h"
#include "utility.h"
#include "timing.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Structure to cache intermediate computations to avoid duplication
struct ParticleIntermediate {
    arma::vec U_hat;        // After ODE step
    arma::vec mean;         // Mean from SDE
    arma::mat cholQ;        // Cholesky of Q
    arma::mat invQ;         // Inverse of Q
    arma::mat cholQxx;      // Cholesky of Qxx (Strang only)
    arma::mat invQxx;       // Inverse of Qxx (Strang only)
    List gaussian_proposal; // Gaussian proposal distribution
};

// [[Rcpp::export]]
arma::vec propagate_langevin_particle_cpp(
    const arma::vec& U,
    const arma::vec& y,
    const arma::mat& M,
    double delta,
    const arma::vec& push,
    const List& potential_params,
    double tau,
    double nu,
    double omega,
    double lambda,
    const std::string& error_dist,
    const List& error_params,
    const std::string& scheme,
    const arma::mat& polygon_coords,
    Nullable<int> ind_fixed_point,
    Nullable<arma::mat> L_provided,
    Nullable<arma::mat> Q_provided,
    double proposal_weight,
    bool verbose
) {

  arma::vec X = U.subvec(0,1);
  arma::mat I2 = arma::eye<arma::mat>(2,2);

  arma::vec U_next(4, arma::fill::zeros);

  if (scheme == "Lie-Trotter") {

    // --- ODE ---
    arma::vec U_hat = solve_ODE_cpp(U, delta, push, potential_params, ind_fixed_point);

    // --- SDE ---
    List OU_solution = solve_SDE_cpp(U_hat, delta, tau, nu, omega,
                                     potential_params,
                                     ind_fixed_point, L_provided, Q_provided);
    arma::mat Q = OU_solution["Q"];
    arma::vec mean = OU_solution["mean"];

    arma::mat cholQ = chol_cpp(Q);
    arma::mat invQ  = chol2inv_cpp(cholQ);

    // --- Gaussian proposal ---
    List gaussian_proposal;

    if (error_dist == "scaled_t") {
      double scale = error_params["scale"];
      double df    = error_params["df"];
      gaussian_proposal = product_gaussian_cpp(invQ, ((df-2)/df)/ (scale*scale)*I2,
                                               mean, y, M, proposal_weight);

    } else if (error_dist == "normal") {
      double sigma_obs = error_params["sigma_obs"];
      gaussian_proposal = product_gaussian_cpp(invQ, 1/(sigma_obs*sigma_obs)*I2,
                                               mean, y, M, proposal_weight);

    } else if (error_dist == "argos") {
      double df    = error_params["df"];
      double sigma_obs = error_params["sigma_obs"];
      double rho   = error_params["rho"];
      double a     = error_params["a"];
      double p     = error_params["p"];

      arma::mat Sigma1 = sigma_obs*sigma_obs * arma::mat{{1, rho*sqrt(a)}, {rho*sqrt(a),1}};
      arma::mat Sigma2 = sigma_obs*sigma_obs * arma::mat{{1, -rho*sqrt(a)}, {-rho*sqrt(a),1}};
      arma::mat S1 = df/(df-2) * Sigma1;
      arma::mat S2 = df/(df-2) * Sigma2;

      arma::mat invS1 = arma::inv(S1);
      arma::mat invS2 = arma::inv(S2);

      double u = R::runif(0.0, 1.0);
      arma::mat invS = (u < p) ? invS1 : invS2;

      gaussian_proposal = product_gaussian_cpp(invQ, invS, mean, y, M, proposal_weight);
    }

    if (verbose) {
      Rcout << "proposal mean: " << as<arma::vec>(gaussian_proposal["mean"]).t() << "\n";
    }

    // --- Sample next state ---
    arma::vec z = arma::randn<arma::vec>(4);
    U_next = as<arma::vec>(gaussian_proposal["mean"]) +
             as<arma::mat>(gaussian_proposal["chol"]) * z;

  } else if (scheme == "Strang") {

    // --- Extract potential parameters ---
    arma::vec alpha = potential_params["alpha"];
    List B = potential_params["B"];
    arma::mat x_star = potential_params["x_star"];

    // --- ODE half-step ---
    arma::vec U_hat = solve_ODE_cpp(U, delta/2.0, push, potential_params, ind_fixed_point);

    // --- SDE ---
    List OU_solution = solve_SDE_cpp(U_hat, delta, tau, nu, omega,
                                     potential_params,
                                     ind_fixed_point, L_provided, Q_provided);
    arma::mat Q = OU_solution["Q"];
    arma::vec mean = OU_solution["mean"];

    arma::mat Qxx = Q.submat(0,0,1,1);
    arma::mat cholQxx = chol_cpp(Qxx);
    arma::mat invQxx = chol2inv_cpp(cholQxx);

    arma::mat Q_v_cond_x = Q.submat(2,2,3,3) -
                           Q.submat(2,0,3,1) * invQxx * Q.submat(0,2,1,3);
    arma::mat cholQ_v_cond_x = chol_cpp(Q_v_cond_x);

    // --- Propagate position ---
    List gaussian_proposal;
    arma::vec X_next(2);

    if (error_dist == "scaled_t") {
      double scale = error_params["scale"];
      double df    = error_params["df"];
      gaussian_proposal = product_gaussian_cpp(invQxx, ((df-2)/df)/(scale*scale)*I2,
                                               mean.subvec(0,1), y, I2, proposal_weight);

    } else if (error_dist == "normal") {
      double sigma_obs = error_params["sigma_obs"];
      gaussian_proposal = product_gaussian_cpp(invQxx, 1/(sigma_obs*sigma_obs)*I2,
                                               mean.subvec(0,1), y, I2, proposal_weight);

    } else if (error_dist == "argos") {
      double df    = error_params["df"];
      double sigma_obs = error_params["sigma_obs"];
      double rho   = error_params["rho"];
      double a     = error_params["a"];
      double p     = error_params["p"];

      arma::mat Sigma1 = sigma_obs*sigma_obs * arma::mat{{1, rho*sqrt(a)}, {rho*sqrt(a),1}};
      arma::mat Sigma2 = sigma_obs*sigma_obs * arma::mat{{1, -rho*sqrt(a)}, {-rho*sqrt(a),1}};
      arma::mat S1 = df/(df-2) * Sigma1;
      arma::mat S2 = df/(df-2) * Sigma2;

      arma::mat invS1 = arma::inv(S1);
      arma::mat invS2 = arma::inv(S2);

      double u = R::runif(0.0, 1.0);
      arma::mat invS = (u < p) ? invS1 : invS2;

      gaussian_proposal = product_gaussian_cpp(invQxx, invS,
                                               mean.subvec(0,1), y, I2, proposal_weight);
    }

    X_next = as<arma::vec>(gaussian_proposal["mean"]) +
             as<arma::mat>(gaussian_proposal["chol"]) * arma::randn<arma::vec>(2);

    // --- Compute push ---
    arma::vec push_next = compute_push_cpp(X_next, polygon_coords, lambda);

    // --- Gradient / nonlinear term ---
    arma::vec grad_term(2);
    if (ind_fixed_point.isNotNull()) {
      int l = as<int>(ind_fixed_point) - 1;
      arma::mat B_l = B[l];
      double alpha_l = alpha[l];
      arma::vec x_star_l = x_star.row(l).t();

      double e_l_next = std::exp(-arma::dot(X_next - x_star_l,
                                            B_l * (X_next - x_star_l)));

      arma::vec grad = mix_gaussian_grad_cpp(X_next, x_star, potential_params,
                                            IntegerVector::create(l+1));
      grad_term = push_next + grad + 2.0*alpha_l*(e_l_next-1.0) * (B_l*(X_next - x_star_l));

    } else {
      arma::vec grad = mix_gaussian_grad_cpp(X_next, x_star, potential_params, IntegerVector());
      grad_term = push_next + grad;
    }

    // --- Propagate velocity ---
    arma::vec m_v_cond_x = mean.subvec(2,3) +
                           Q.submat(2,0,3,1) * invQxx * (X_next - mean.subvec(0,1)) -
                           delta/2.0 * grad_term;

    arma::vec V_next = m_v_cond_x + cholQ_v_cond_x * arma::randn<arma::vec>(2);

    U_next = arma::join_vert(X_next, V_next);
  }

  return U_next;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double compute_langevin_weight_cpp(const arma::vec& U_pred,
                                   const arma::vec& U_prev,
                                   const arma::vec& y,
                                   const arma::mat& M,
                                   double delta,
                                   const arma::vec& push,
                                   const List& potential_params,
                                   double tau,
                                   double nu,
                                   double omega,
                                   const std::string& error_dist,
                                   const List& error_params,
                                   const std::string& scheme,
                                   Nullable<int> ind_fixed_point,
                                   Nullable<arma::mat> L,
                                   Nullable<arma::mat> Q,
                                   double proposal_weight,
                                   bool verbose) {

  // --- Identity 2x2 ---
  arma::mat I2 = arma::eye(2,2);

  // --- Step 1: Solve ODE ---
  arma::vec U_hat;
  if (scheme == "Lie-Trotter") {
    U_hat = solve_ODE_cpp(U_prev, delta, push, potential_params, ind_fixed_point);
  } else { // Strang
    U_hat = solve_ODE_cpp(U_prev, delta / 2.0, push, potential_params, ind_fixed_point);
  }

  // --- Step 2: Solve SDE ---
  List OU_solution = solve_SDE_cpp(U_hat, delta, tau, nu, omega, potential_params,
                                   ind_fixed_point, L, Q);

  arma::mat Q_mat = OU_solution["Q"];
  arma::vec mean  = OU_solution["mean"];

  arma::mat cholQ  = chol_cpp(Q_mat);
  arma::mat invQ   = chol2inv_cpp(cholQ);
  arma::mat Qxx, cholQxx, invQxx;

  if (scheme == "Strang") {
    Qxx = Q_mat.submat(0,0,1,1);
    cholQxx = chol_cpp(Qxx);
    invQxx  = chol2inv_cpp(cholQxx);
  }

  double log_weight = 0.0;
  double local_llk = 0.0;
  double log_pred_density = 0.0;
  double log_prop_density = 0.0;

  // --- Step 3: Compute proposal depending on error distribution ---
  if (error_dist == "scaled_t") {
    double scale = as<double>(error_params["scale"]);
    double df    = as<double>(error_params["df"]);

    List gaussian_proposal = product_gaussian_cpp(invQ,
                                                  (df-2)/df/scale/scale * I2,
                                                  mean, y, M, proposal_weight);

    local_llk = dscaledt_cpp(y[0], U_pred[0], scale, df, true) +
                dscaledt_cpp(y[1], U_pred[1], scale, df, true);

    log_pred_density = log_dmvnorm_chol_cpp(U_pred, mean, cholQ);
    log_prop_density = log_dmvnorm_chol_cpp(U_pred,
                                            gaussian_proposal["mean"],
                                            gaussian_proposal["chol"]);

    log_weight = local_llk + log_pred_density - log_prop_density;
  }
  else if (error_dist == "normal") {
    double sigma_obs = as<double>(error_params["sigma_obs"]);

    List gaussian_proposal;
    if (scheme == "Lie-Trotter") {
      gaussian_proposal = product_gaussian_cpp(invQ, (1.0/(sigma_obs*sigma_obs)) * I2,
                                               mean, y, M, proposal_weight);
      local_llk = log_dmvnorm_chol_cpp(y, U_pred.subvec(0,1), sigma_obs*I2);
      log_pred_density = log_dmvnorm_chol_cpp(U_pred, mean, cholQ);
      log_prop_density = log_dmvnorm_chol_cpp(U_pred,
                                              gaussian_proposal["mean"],
                                              gaussian_proposal["chol"]);
    } else {
      gaussian_proposal = product_gaussian_cpp(invQxx, (1.0/(sigma_obs*sigma_obs)) * I2,
                                               mean.subvec(0,1), y, I2, proposal_weight);
      local_llk = log_dmvnorm_chol_cpp(y, U_pred.subvec(0,1), sigma_obs*I2);
      log_pred_density = log_dmvnorm_chol_cpp(U_pred.subvec(0,1), mean.subvec(0,1), cholQxx);
      log_prop_density = log_dmvnorm_chol_cpp(U_pred.subvec(0,1),
                                              gaussian_proposal["mean"],
                                              gaussian_proposal["chol"]);
    }
    log_weight = local_llk + log_pred_density - log_prop_density;
  }
  else if (error_dist == "argos") {
    double df        = as<double>(error_params["df"]);
    double sigma_obs = as<double>(error_params["sigma_obs"]);
    double rho       = as<double>(error_params["rho"]);
    double a         = as<double>(error_params["a"]);
    double p         = as<double>(error_params["p"]);

    arma::mat Sigma1 = sigma_obs*sigma_obs * arma::mat{{1, rho*std::sqrt(a)},
                                                       {rho*std::sqrt(a),1}};
    arma::mat Sigma2 = sigma_obs*sigma_obs * arma::mat{{1, -rho*std::sqrt(a)},
                                                       {-rho*std::sqrt(a),1}};

    arma::mat S1 = df/(df-2.0) * Sigma1;
    arma::mat S2 = df/(df-2.0) * Sigma2;

    List gaussian_proposal1 = product_gaussian_cpp(invQ, chol2inv_cpp(chol_cpp(S1)),
                                                   mean, y, M, proposal_weight);
    List gaussian_proposal2 = product_gaussian_cpp(invQ, chol2inv_cpp(chol_cpp(S2)),
                                                   mean, y, M, proposal_weight);

    local_llk = dmvt_mixture_cpp(y, U_pred.subvec(0,1), error_params, true);
    log_pred_density = log_dmvnorm_chol_cpp(U_pred, mean, cholQ);

    double l1 = std::log(p) +
                log_dmvnorm_chol_cpp(U_pred, gaussian_proposal1["mean"], gaussian_proposal1["chol"]);
    double l2 = std::log(1-p) +
                log_dmvnorm_chol_cpp(U_pred, gaussian_proposal2["mean"], gaussian_proposal2["chol"]);

    double m = std::max(l1,l2);
    log_prop_density = m + std::log(std::exp(l1-m) + std::exp(l2-m));
    log_weight = local_llk + log_pred_density - log_prop_density;
  }

  if (verbose) {
    Rcout << "local_llk = " << std::exp(local_llk)
          << ", pred_density = " << std::exp(log_pred_density)
          << ", prop_density = " << std::exp(log_prop_density) << std::endl;
  }

  return std::exp(log_weight);
}

// Propagate particle and cache intermediate computations for weight calculation
// Properly handles both Lie-Trotter and Strang schemes
static inline arma::vec propagate_particle_with_cache(
    const arma::vec& U_prev,
    const arma::vec& y,
    const arma::mat& M,
    double delta,
    const arma::vec& push,
    const List& potential_params,
    double tau,
    double nu,
    double omega,
    double lambda,
    const std::string& error_dist,
    const List& error_params,
    const std::string& scheme,
    const arma::mat& polygon_coords,
    int ind_fixed_point,
    bool use_precomputed_LQ,
    const arma::mat& L_precomputed,
    const arma::mat& Q_precomputed,
    double proposal_weight,
    ParticleIntermediate& cache  // OUTPUT: cached values
) {
  
  static const arma::mat I2 = arma::eye<arma::mat>(2, 2);
  arma::vec U_next(4, arma::fill::zeros);
  
  if (scheme == "Lie-Trotter") {
    
    // --- ODE FULL step ---
    global_timer.start();
    if (ind_fixed_point > 0) {
      cache.U_hat = solve_ODE_cpp(U_prev, delta, push, potential_params, wrap(ind_fixed_point));
    } else {
      cache.U_hat = solve_ODE_cpp(U_prev, delta, push, potential_params, R_NilValue);
    }
    global_timer.record("prop_ode");
    
    // --- SDE step ---
    List OU_solution;
    if (use_precomputed_LQ) {
      if (ind_fixed_point > 0) {
        OU_solution = solve_SDE_cpp(cache.U_hat, delta, tau, nu, omega, potential_params,
                                    wrap(ind_fixed_point), wrap(L_precomputed), wrap(Q_precomputed));
      } else {
        OU_solution = solve_SDE_cpp(cache.U_hat, delta, tau, nu, omega, potential_params,
                                    R_NilValue, wrap(L_precomputed), wrap(Q_precomputed));
      }
    } else {
      if (ind_fixed_point > 0) {
        OU_solution = solve_SDE_cpp(cache.U_hat, delta, tau, nu, omega, potential_params,
                                    wrap(ind_fixed_point), R_NilValue, R_NilValue);
      } else {
        OU_solution = solve_SDE_cpp(cache.U_hat, delta, tau, nu, omega, potential_params,
                                    R_NilValue, R_NilValue, R_NilValue);
      }
    }
    global_timer.record("prop_sde");
    
    arma::mat Q = as<arma::mat>(OU_solution["Q"]);
    cache.mean = as<arma::vec>(OU_solution["mean"]);
    
    // Cache Cholesky decompositions
    cache.cholQ = chol_cpp(Q);
    cache.invQ = chol2inv_cpp(cache.cholQ);
    global_timer.record("prop_cholesky");
    
    // Cache 4D Gaussian proposal
    if (error_dist == "scaled_t") {
      double scale = error_params["scale"];
      double df = error_params["df"];
      cache.gaussian_proposal = product_gaussian_cpp(cache.invQ, ((df-2)/df)/(scale*scale)*I2,
                                                      cache.mean, y, M, proposal_weight);
    } else if (error_dist == "normal") {
      double sigma_obs = error_params["sigma_obs"];
      cache.gaussian_proposal = product_gaussian_cpp(cache.invQ, 1/(sigma_obs*sigma_obs)*I2,
                                                      cache.mean, y, M, proposal_weight);
    } else if (error_dist == "argos") {
      double df = error_params["df"];
      double sigma_obs = error_params["sigma_obs"];
      double rho = error_params["rho"];
      double a = error_params["a"];
      double p = error_params["p"];
      
      arma::mat Sigma1 = sigma_obs*sigma_obs * arma::mat{{1, rho*sqrt(a)}, {rho*sqrt(a), 1}};
      arma::mat Sigma2 = sigma_obs*sigma_obs * arma::mat{{1, -rho*sqrt(a)}, {-rho*sqrt(a), 1}};
      arma::mat S1 = df/(df-2) * Sigma1;
      arma::mat S2 = df/(df-2) * Sigma2;
      
      arma::mat invS1 = arma::inv(S1);
      arma::mat invS2 = arma::inv(S2);
      
      double u = R::runif(0.0, 1.0);
      arma::mat invS = (u < p) ? invS1 : invS2;
      
      cache.gaussian_proposal = product_gaussian_cpp(cache.invQ, invS, cache.mean, y, M, proposal_weight);
    }
    global_timer.record("prop_proposal");
    
    // Sample from cached 4D proposal
    arma::vec z = arma::randn<arma::vec>(4);
    U_next = as<arma::vec>(cache.gaussian_proposal["mean"]) +
             as<arma::mat>(cache.gaussian_proposal["chol"]) * z;
    global_timer.record("prop_sample");
    
  } else if (scheme == "Strang") {
    
    // --- Extract potential parameters ---
    arma::vec alpha = potential_params["alpha"];
    List B = potential_params["B"];
    arma::mat x_star = potential_params["x_star"];
    
    // --- ODE HALF-step ---
    global_timer.start();
    if (ind_fixed_point > 0) {
      cache.U_hat = solve_ODE_cpp(U_prev, delta/2.0, push, potential_params, wrap(ind_fixed_point));
    } else {
      cache.U_hat = solve_ODE_cpp(U_prev, delta/2.0, push, potential_params, R_NilValue);
    }
    global_timer.record("prop_ode");
    
    // --- SDE step ---
    List OU_solution;
    if (use_precomputed_LQ) {
      if (ind_fixed_point > 0) {
        OU_solution = solve_SDE_cpp(cache.U_hat, delta, tau, nu, omega, potential_params,
                                    wrap(ind_fixed_point), wrap(L_precomputed), wrap(Q_precomputed));
      } else {
        OU_solution = solve_SDE_cpp(cache.U_hat, delta, tau, nu, omega, potential_params,
                                    R_NilValue, wrap(L_precomputed), wrap(Q_precomputed));
      }
    } else {
      if (ind_fixed_point > 0) {
        OU_solution = solve_SDE_cpp(cache.U_hat, delta, tau, nu, omega, potential_params,
                                    wrap(ind_fixed_point), R_NilValue, R_NilValue);
      } else {
        OU_solution = solve_SDE_cpp(cache.U_hat, delta, tau, nu, omega, potential_params,
                                    R_NilValue, R_NilValue, R_NilValue);
      }
    }
    global_timer.record("prop_sde");
    
    arma::mat Q = as<arma::mat>(OU_solution["Q"]);
    cache.mean = as<arma::vec>(OU_solution["mean"]);
    
    // Cache Cholesky for position block
    arma::mat Qxx = Q.submat(0, 0, 1, 1);
    cache.cholQxx = chol_cpp(Qxx);
    cache.invQxx = chol2inv_cpp(cache.cholQxx);
    
    // Compute conditional covariance for velocity
    arma::mat Q_v_cond_x = Q.submat(2,2,3,3) -
                           Q.submat(2,0,3,1) * cache.invQxx * Q.submat(0,2,1,3);
    arma::mat cholQ_v_cond_x = chol_cpp(Q_v_cond_x);
    global_timer.record("prop_cholesky");
    
    // --- Propagate position only (2D proposal) ---
    if (error_dist == "scaled_t") {
      double scale = error_params["scale"];
      double df = error_params["df"];
      cache.gaussian_proposal = product_gaussian_cpp(cache.invQxx, ((df-2)/df)/(scale*scale)*I2,
                                                      cache.mean.subvec(0,1), y, I2, proposal_weight);
    } else if (error_dist == "normal") {
      double sigma_obs = error_params["sigma_obs"];
      cache.gaussian_proposal = product_gaussian_cpp(cache.invQxx, 1/(sigma_obs*sigma_obs)*I2,
                                                      cache.mean.subvec(0,1), y, I2, proposal_weight);
    } else if (error_dist == "argos") {
      double df = error_params["df"];
      double sigma_obs = error_params["sigma_obs"];
      double rho = error_params["rho"];
      double a = error_params["a"];
      double p = error_params["p"];
      
      arma::mat Sigma1 = sigma_obs*sigma_obs * arma::mat{{1, rho*sqrt(a)}, {rho*sqrt(a), 1}};
      arma::mat Sigma2 = sigma_obs*sigma_obs * arma::mat{{1, -rho*sqrt(a)}, {-rho*sqrt(a), 1}};
      arma::mat S1 = df/(df-2) * Sigma1;
      arma::mat S2 = df/(df-2) * Sigma2;
      
      arma::mat invS1 = arma::inv(S1);
      arma::mat invS2 = arma::inv(S2);
      
      double u = R::runif(0.0, 1.0);
      arma::mat invS = (u < p) ? invS1 : invS2;
      
      cache.gaussian_proposal = product_gaussian_cpp(cache.invQxx, invS,
                                                      cache.mean.subvec(0,1), y, I2, proposal_weight);
    }
    global_timer.record("prop_proposal");
    
    // Sample position
    arma::vec X_next = as<arma::vec>(cache.gaussian_proposal["mean"]) +
                       as<arma::mat>(cache.gaussian_proposal["chol"]) * arma::randn<arma::vec>(2);
    
    // --- Compute push at new position ---
    arma::vec push_next = compute_push_cpp(X_next, polygon_coords, lambda);
    
    // --- Recompute gradient / nonlinear term at new position ---
    arma::vec grad_term(2);
    if (ind_fixed_point > 0) {
      int l = ind_fixed_point - 1;
      arma::mat B_l = B[l];
      double alpha_l = alpha[l];
      arma::vec x_star_l = x_star.row(l).t();
      
      double e_l_next = std::exp(-arma::dot(X_next - x_star_l,
                                            B_l * (X_next - x_star_l)));
      
      arma::vec grad = mix_gaussian_grad_cpp(X_next, x_star, potential_params,
                                            IntegerVector::create(l+1));
      grad_term = push_next + grad + 2.0*alpha_l*(e_l_next-1.0) * (B_l*(X_next - x_star_l));
      
    } else {
      arma::vec grad = mix_gaussian_grad_cpp(X_next, x_star, potential_params, IntegerVector());
      grad_term = push_next + grad;
    }
    
    // --- Propagate velocity conditionally ---
    arma::vec m_v_cond_x = cache.mean.subvec(2,3) +
                           Q.submat(2,0,3,1) * cache.invQxx * (X_next - cache.mean.subvec(0,1)) -
                           delta/2.0 * grad_term;
    
    arma::vec V_next = m_v_cond_x + cholQ_v_cond_x * arma::randn<arma::vec>(2);
    global_timer.record("prop_sample");
    
    U_next = arma::join_vert(X_next, V_next);
  }
  
  return U_next;
}

// Compute weight using cached intermediate values (avoids recomputation)
// Properly handles both Lie-Trotter and Strang schemes
static inline double compute_weight_from_cache(
    const arma::vec& U_pred,
    const arma::vec& y,
    const arma::mat& M,
    const std::string& error_dist,
    const std::string& scheme,
    double proposal_weight,
    // Precomputed error parameters
    double sigma_obs,
    double scale,
    double df,
    const arma::mat& invS_argos1,
    const arma::mat& invS_argos2,
    double log_p_argos,
    double log_1mp_argos,
    // Cached intermediate values
    const ParticleIntermediate& cache
) {
  
  static const arma::mat I2 = arma::eye<arma::mat>(2, 2);
  
  // Skip ODE, SDE - use cached values!
  // Already have: cache.U_hat, cache.mean, cache.cholQ/cholQxx, cache.invQ/invQxx, cache.gaussian_proposal
  
  double log_weight = 0.0;
  double local_llk = 0.0;
  double log_pred_density = 0.0;
  double log_prop_density = 0.0;
  
  // Compute densities using cached values
  if (scheme == "Lie-Trotter") {
    
    // Lie-Trotter: 4D state, use cached 4D proposal
    if (error_dist == "scaled_t") {
      local_llk = dscaledt_cpp(y[0], U_pred[0], scale, df, true) +
                  dscaledt_cpp(y[1], U_pred[1], scale, df, true);
      
      log_pred_density = log_dmvnorm_chol_cpp(U_pred, cache.mean, cache.cholQ);
      log_prop_density = log_dmvnorm_chol_cpp(U_pred,
                                              cache.gaussian_proposal["mean"],
                                              cache.gaussian_proposal["chol"]);
      global_timer.record("weight_densities");
      
    } else if (error_dist == "normal") {
      local_llk = log_dmvnorm_chol_cpp(y, U_pred.subvec(0, 1), sigma_obs * I2);
      log_pred_density = log_dmvnorm_chol_cpp(U_pred, cache.mean, cache.cholQ);
      log_prop_density = log_dmvnorm_chol_cpp(U_pred,
                                              cache.gaussian_proposal["mean"],
                                              cache.gaussian_proposal["chol"]);
      global_timer.record("weight_densities");
      
    } else if (error_dist == "argos") {
      // For argos, we still need to compute proposals (different random draw in weight computation)
      List gaussian_proposal1 = product_gaussian_cpp(cache.invQ, invS_argos1,
                                                     cache.mean, y, M, proposal_weight);
      List gaussian_proposal2 = product_gaussian_cpp(cache.invQ, invS_argos2,
                                                     cache.mean, y, M, proposal_weight);
      global_timer.record("weight_proposal");
      
      List error_params_local = List::create(
        Named("sigma_obs") = sigma_obs,
        Named("df") = df,
        Named("rho") = 0.0,
        Named("a") = 0.0,
        Named("p") = std::exp(log_p_argos)
      );
      local_llk = dmvt_mixture_cpp(y, U_pred.subvec(0, 1), error_params_local, true);
      log_pred_density = log_dmvnorm_chol_cpp(U_pred, cache.mean, cache.cholQ);
      
      double l1 = log_p_argos +
                  log_dmvnorm_chol_cpp(U_pred, gaussian_proposal1["mean"], gaussian_proposal1["chol"]);
      double l2 = log_1mp_argos +
                  log_dmvnorm_chol_cpp(U_pred, gaussian_proposal2["mean"], gaussian_proposal2["chol"]);
      
      double m = std::max(l1, l2);
      log_prop_density = m + std::log(std::exp(l1 - m) + std::exp(l2 - m));
      global_timer.record("weight_densities");
    }
    
    log_weight = local_llk + log_pred_density - log_prop_density;
    
  } else if (scheme == "Strang") {
    
    // Strang: 2D position, use cached 2D proposal
    if (error_dist == "scaled_t") {
      local_llk = dscaledt_cpp(y[0], U_pred[0], scale, df, true) +
                  dscaledt_cpp(y[1], U_pred[1], scale, df, true);
      
      log_pred_density = log_dmvnorm_chol_cpp(U_pred.subvec(0, 1), cache.mean.subvec(0, 1), cache.cholQxx);
      log_prop_density = log_dmvnorm_chol_cpp(U_pred.subvec(0, 1),
                                              cache.gaussian_proposal["mean"],
                                              cache.gaussian_proposal["chol"]);
      global_timer.record("weight_densities");
      
    } else if (error_dist == "normal") {
      local_llk = log_dmvnorm_chol_cpp(y, U_pred.subvec(0, 1), sigma_obs * I2);
      log_pred_density = log_dmvnorm_chol_cpp(U_pred.subvec(0, 1), cache.mean.subvec(0, 1), cache.cholQxx);
      log_prop_density = log_dmvnorm_chol_cpp(U_pred.subvec(0, 1),
                                              cache.gaussian_proposal["mean"],
                                              cache.gaussian_proposal["chol"]);
      global_timer.record("weight_densities");
      
    } else if (error_dist == "argos") {
      // For argos, compute 2D proposals
      List gaussian_proposal1 = product_gaussian_cpp(cache.invQxx, invS_argos1,
                                                     cache.mean.subvec(0,1), y, I2, proposal_weight);
      List gaussian_proposal2 = product_gaussian_cpp(cache.invQxx, invS_argos2,
                                                     cache.mean.subvec(0,1), y, I2, proposal_weight);
      global_timer.record("weight_proposal");
      
      List error_params_local = List::create(
        Named("sigma_obs") = sigma_obs,
        Named("df") = df,
        Named("rho") = 0.0,
        Named("a") = 0.0,
        Named("p") = std::exp(log_p_argos)
      );
      local_llk = dmvt_mixture_cpp(y, U_pred.subvec(0, 1), error_params_local, true);
      log_pred_density = log_dmvnorm_chol_cpp(U_pred.subvec(0, 1), cache.mean.subvec(0, 1), cache.cholQxx);
      
      double l1 = log_p_argos +
                  log_dmvnorm_chol_cpp(U_pred.subvec(0, 1), gaussian_proposal1["mean"], gaussian_proposal1["chol"]);
      double l2 = log_1mp_argos +
                  log_dmvnorm_chol_cpp(U_pred.subvec(0, 1), gaussian_proposal2["mean"], gaussian_proposal2["chol"]);
      
      double m = std::max(l1, l2);
      log_prop_density = m + std::log(std::exp(l1 - m) + std::exp(l2 - m));
      global_timer.record("weight_densities");
    }
    
    log_weight = local_llk + log_pred_density - log_prop_density;
  }
  
  return std::exp(log_weight);
}


//' Run particle filter algorithm for the Langevin diffusion with specified measurement errors
//' and splitting schemes
//' 
//' @param observations Matrix with dimensions N x 3: time, Y1, Y2
//' @param sde_params List with elements: tau (persistence), nu (velocity scale), omega (rotation)
//' @param potential_params List with elements:
//'   \itemize{
//'     \item alpha: Numeric vector of mixture weights
//'     \item B: List of precision matrices
//'     \item x_star: Matrix of fixed point coordinates
//'   }
//' @param error_params List of parameters for the measurement error distribution (contents depend on error_dist):
//'   \itemize{
//'     \item For "normal": sigma_obs (observation error standard deviation)
//'     \item For "scaled_t": scale and df (degrees of freedom)
//'     \item For "argos": sigma_obs, df, rho, a, p
//'   }
//' @param error_dist String specifying error distribution: "normal", "scaled_t", or "argos"
//' @param polygon_coords Matrix of polygon boundary coordinates (N x 2)
//' @param U0 Numeric vector of initial state (length 4: position and velocity)
//' @param lambda Numeric, penalization parameter (positive number, Inf for no penalization)
//' @param num_particles Integer, number of particles to use in the filter
//' @param scheme String, splitting scheme: "Lie-Trotter" or "Strang"
//' @param split_around_fixed_point Logical, whether to split around fixed points
//' @param ESS_threshold Numeric between 0 and 1, threshold on effective sample size for resampling
//' @param proposal_weight Numeric between 0 and 1, weight attributed to the previous state in the Gaussian proposal
//' @param verbose Logical, whether to print progress messages
//' 
//' @return List with elements:
//'   \itemize{
//'     \item particles: Array of particles (num_particles x 4 x N)
//'     \item weights: Matrix of normalized weights (num_particles x N)
//'     \item total_weights: Vector of sum of unnormalized weights at each time step
//'     \item loglik: Estimated log-likelihood
//'     \item loglik_vector: Estimated log-likelihood at each time step
//'     \item push: Array of inward pushes for each particle at each time step
//'     \item ancestors: Matrix of ancestor indices
//'     \item resampled_at: Logical vector indicating steps where resampling was performed
//'     \item ess_history: Vector of ESS at each time step
//'     \item timing: List of timing information for code profiling
//'     \item ind_fixed_point: (if split_around_fixed_point=TRUE) Matrix of fixed point indices
//'   }
//' @export
// [[Rcpp::export]]
List particle_filter2D_cpp(
    const arma::mat& observations,           
    const List& sde_params,                  
    const List& potential_params,           
    const List& error_params,                
    const std::string& error_dist,           
    const arma::mat& polygon_coords,         
    const arma::vec& U0,                     
    double lambda,                          
    int num_particles,
    const std::string& scheme,              
    bool split_around_fixed_point,
    double ESS_threshold,
    double proposal_weight,
    bool verbose
) {
  
  // Initialize timing
  global_timer.reset();
  global_timer.start();
  
  // Extract dimensions
  int N = observations.n_rows;
  int num_states = 4;
  
  // Extract SDE parameters
  double tau = as<double>(sde_params["tau"]);
  double nu = as<double>(sde_params["nu"]);
  double omega = as<double>(sde_params["omega"]);
  
  // Observation link matrix
  arma::mat I2 = arma::eye<arma::mat>(2, 2);
  arma::mat M = arma::join_horiz(I2, arma::zeros<arma::mat>(2, 2));
  
  // Extract potential parameters if provided
  arma::vec alpha = potential_params["alpha"];
  List B_list = potential_params["B"];
  arma::mat x_star = potential_params["x_star"];
  
  // Initialize output arrays
  arma::cube particles(num_particles, num_states, N);
  arma::mat weights(num_particles, N);
  arma::vec total_weights(N);
  arma::vec loglik_vector(N - 1);
  arma::cube push_array(num_particles, 2, N - 1);
  arma::imat ind_fixed_point_mat;
  arma::imat ancestors(num_particles, N - 1);
  arma::vec ess_history(N);
  LogicalVector resampled_at(N - 1);
  
  if (split_around_fixed_point) {
    ind_fixed_point_mat = arma::imat(num_particles, N - 1);
  }
  
  double loglik = 0.0;
  
  // Initialize particles with noise around U0
  arma::mat R0 = arma::diagmat(arma::vec{0.01, 0.01, 0.25, 0.25});
  for (int k = 0; k < num_particles; ++k) {
    arma::vec noise = arma::sqrtmat_sympd(R0) * arma::randn<arma::vec>(num_states);
    particles.slice(0).row(k) = (U0 + noise).t();
  }
  
  // Initialize weights
  weights.col(0).fill(1.0 / num_particles);
  total_weights(0) = 1.0;
  ess_history(0) = num_particles;
  
  if (verbose) {
    Rcout << "Initialization complete. Starting particle filtering...\n";
    Rcout << "Running particle filter with tau=" << tau 
          << ", nu=" << nu << ", omega=" << omega << "\n";
  }
  
  global_timer.record("initialization");
  
  // Precompute time steps
  arma::vec times = observations.col(0);
  arma::vec deltas = arma::diff(times);
  
  // Precompute L and Q if time step is constant and not splitting
  bool use_precomputed_LQ = false;
  arma::mat L_precomputed;
  arma::mat Q_precomputed;
  
  if (arma::max(deltas) - arma::min(deltas) < 1e-6 && !split_around_fixed_point) {
    use_precomputed_LQ = true;
    double delta = deltas(0);
    L_precomputed = RACVM_link_cpp(tau, omega, delta);
    Q_precomputed = RACVM_cov_cpp(tau, nu, omega, delta);
    if (verbose) {
      Rcout << "Using precomputed L and Q matrices (constant time step = " 
            << delta << ").\n";
    }
  }
  
  global_timer.record("precompute_LQ");
  
  // Precompute error distribution parameters to avoid repeated extraction
  double sigma_obs = 0.0;
  double scale = 0.0, df = 0.0;
  double rho = 0.0, a = 0.0, p_argos = 0.0;
  arma::mat invS_argos1, invS_argos2;
  double log_p_argos = 0.0, log_1mp_argos = 0.0;
  
  if (error_dist == "normal") {
    sigma_obs = as<double>(error_params["sigma_obs"]);
  } else if (error_dist == "scaled_t") {
    scale = as<double>(error_params["scale"]);
    df = as<double>(error_params["df"]);
  } else if (error_dist == "argos") {
    df = as<double>(error_params["df"]);
    sigma_obs = as<double>(error_params["sigma_obs"]);
    rho = as<double>(error_params["rho"]);
    a = as<double>(error_params["a"]);
    p_argos = as<double>(error_params["p"]);
    
    // Precompute inverse covariance matrices
    arma::mat Sigma1 = sigma_obs*sigma_obs * arma::mat{{1, rho*std::sqrt(a)},
                                                       {rho*std::sqrt(a), 1}};
    arma::mat Sigma2 = sigma_obs*sigma_obs * arma::mat{{1, -rho*std::sqrt(a)},
                                                       {-rho*std::sqrt(a), 1}};
    arma::mat S1 = df/(df-2.0) * Sigma1;
    arma::mat S2 = df/(df-2.0) * Sigma2;
    
    invS_argos1 = arma::inv(S1);
    invS_argos2 = arma::inv(S2);
    log_p_argos = std::log(p_argos);
    log_1mp_argos = std::log(1.0 - p_argos);
  }
  
  // Main particle filter loop
  for (int j = 0; j < N - 1; ++j) {
    
    if (verbose) {
      Rcout << "Time step: " << j + 1 << " / " << N - 1 << "\n";
    }
    
    double delta = deltas(j);
    arma::vec y = observations.row(j + 1).subvec(1, 2).t();
    
    // --- PREDICTION STEP ---
    if (verbose) Rcout << "  Prediction step...\n";
    
    // Extract current particle positions
    arma::mat X_positions = particles.slice(j).cols(0, 1);
    
    // Compute push for all particles (vectorized)
    arma::mat push_matrix = compute_push_matrix_cpp(X_positions, polygon_coords, lambda);
    push_array.slice(j) = push_matrix;
    
    global_timer.record("compute_push");
    
    // Choose centers for all particles if splitting
    arma::ivec ind_fixed_point_vec;
    if (split_around_fixed_point) {
      ind_fixed_point_vec = choose_center_matrix_cpp(
        X_positions,
        x_star,
        potential_params
      );
      ind_fixed_point_mat.col(j) = ind_fixed_point_vec;
      global_timer.record("choose_centers");
    }
    
    // Cache intermediate values for each particle to avoid recomputation
    std::vector<ParticleIntermediate> particle_cache(num_particles);
    
    // --- PREDICTION STEP ---
    // Propagate each particle and cache intermediate computations
    for (int k = 0; k < num_particles; ++k) {
      
      arma::vec U_prev = particles.slice(j).row(k).t();
      arma::vec push = push_matrix.row(k).t();
      
      // Get fixed point index for this particle
      int ind_fp = (split_around_fixed_point && ind_fixed_point_vec(k) > 0) 
                ? ind_fixed_point_vec(k) : 0;
      
      // Propagate particle WITH CACHING
      arma::vec U_next = propagate_particle_with_cache(
        U_prev, y, M, delta, push,
        potential_params, tau, nu, omega, lambda,
        error_dist, error_params, scheme, polygon_coords,
        ind_fp, use_precomputed_LQ, L_precomputed, Q_precomputed,
        proposal_weight,
        particle_cache[k]);  // OUTPUT: cache intermediate values
      
      particles.slice(j + 1).row(k) = U_next.t();
    }
    
    global_timer.record("propagate_particles");
    
    if (verbose) Rcout << "  Prediction step complete.\n";
  
   // --- CORRECTION STEP ---
    if (verbose) Rcout << "  Correction step...\n";
    
    // Compute weights using CACHED intermediate values
    for (int k = 0; k < num_particles; ++k) {
      
      arma::vec U_pred = particles.slice(j + 1).row(k).t();
      
      // Use cached values - NO recomputation of ODE, SDE, Cholesky!
      weights(k, j + 1) = compute_weight_from_cache(
        U_pred, y, M,
        error_dist, scheme,
        proposal_weight, sigma_obs, scale, df,
        invS_argos1, invS_argos2, log_p_argos, log_1mp_argos,
        particle_cache[k]);  // INPUT: use cached values
     }
    
    global_timer.record("compute_weights");
    
    // Store total weight
    total_weights(j + 1) = arma::sum(weights.col(j + 1));
    
    // Update log-likelihood
    loglik += std::log(total_weights(j + 1)) - std::log(num_particles);
    loglik_vector(j) = std::log(total_weights(j + 1)) - std::log(num_particles);
    
    // Normalize weights
    double weight_sum = total_weights(j + 1);
    
    if (weight_sum > 0) {
      weights.col(j + 1) /= weight_sum;
      ess_history(j + 1) = 1.0 / arma::sum(arma::square(weights.col(j + 1)));
      if (verbose) Rcout << "  Correction step complete. ESS = " << ess_history(j + 1) << "\n";
    } else {
      Rcpp::warning("Sum of weights is zero at time step %d. Reinitializing particles...", j + 1);
      
      // Reinitialize around observation
      for (int k = 0; k < num_particles; ++k) {
        arma::vec noise = arma::sqrtmat_sympd(R0) * arma::randn<arma::vec>(num_states);
        arma::vec new_particle = arma::join_vert(y, arma::zeros<arma::vec>(2)) + noise;
        particles.slice(j + 1).row(k) = new_particle.t();
      }
      weights.col(j + 1).fill(1.0 / num_particles);
      ess_history(j + 1) = num_particles;
    }
    
    global_timer.record("normalize_weights");
    
    // --- RESAMPLING STEP ---
    if (j < N - 2) {
      double ess_ratio = ess_history(j + 1) / num_particles;
      
      if (ess_ratio < ESS_threshold) {
        // Resample using systematic resampling (much faster than naive search)
        resampled_at(j) = true;
        if (verbose) Rcout << "  Resampling (ESS ratio = " << ess_ratio << ")\n";
        
        // Systematic resampling - O(N) instead of O(N²)
        arma::vec cum_weights = arma::cumsum(weights.col(j + 1));
        double u0 = R::runif(0.0, 1.0 / num_particles);
        
        arma::mat particles_resampled = particles.slice(j + 1);
        int idx = 0;
        
        for (int k = 0; k < num_particles; ++k) {
          double u = u0 + k / double(num_particles);
          
          // Move idx forward until cum_weights[idx] >= u
          while (idx < num_particles - 1 && cum_weights(idx) < u) {
            idx++;
          }
          
          particles.slice(j + 1).row(k) = particles_resampled.row(idx);
          ancestors(k, j) = idx + 1;  // R-style 1-based indexing
        }
        
        // Reset weights
        weights.col(j + 1).fill(1.0 / num_particles);
        
        global_timer.record("resampling");
        
      } else {
        resampled_at(j) = false;
        for (int k = 0; k < num_particles; ++k) {
          ancestors(k, j) = k + 1;  // No resampling, keep same indices
        }
        if (verbose) Rcout << "  No resampling (ESS ratio = " << ess_ratio << ")\n";
      }
    }
  }  // End main loop
  
  if (verbose) Rcout << "Particle filtering complete.\n";
  
  // Print timing results
  global_timer.print_timings();
  
  // Prepare output
  List result = List::create(
    Named("particles") = particles,
    Named("weights") = weights,
    Named("total_weights") = total_weights,
    Named("loglik") = loglik,
    Named("loglik_vector") = loglik_vector,
    Named("push") = push_array,
    Named("ancestors") = ancestors,
    Named("resampled_at") = resampled_at,
    Named("ess_history") = ess_history,
    Named("timing") = global_timer.get_timings()
  );
  
  if (split_around_fixed_point) {
    result["ind_fixed_point"] = ind_fixed_point_mat;
  }
  
  return result;
}



