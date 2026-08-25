#include "particle_filter.h"
#include "sde.h"
#include "ode.h"
#include "utility.h"
#include "timing.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Structure to cache intermediate computations to avoid computing same matrices two times
struct ParticleIntermediate {
    arma::vec U_hat;        // After ODE step
    arma::vec mean;         // Mean from SDE
    arma::mat cholQ;        // Cholesky of Q
    arma::mat invQ;         // Inverse of Q
    arma::mat cholQxx;      // Cholesky of Qxx (Strang only)
    arma::mat invQxx;       // Inverse of Qxx (Strang only)
    List gaussian_proposal; // Gaussian proposal distribution
};

// Propagate particle and cache intermediate computations for weight calculation
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
             as<arma::mat>(cache.gaussian_proposal["chol"]).t() * z;
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
                       as<arma::mat>(cache.gaussian_proposal["chol"]).t() * arma::randn<arma::vec>(2);
    
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
    
    arma::vec V_next = m_v_cond_x + cholQ_v_cond_x.t() * arma::randn<arma::vec>(2);
    global_timer.record("prop_sample");
    
    U_next = arma::join_vert(X_next, V_next);
  }
  
  return U_next;
}

// Compute weight using cached intermediate values (avoids recomputation)
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
    double a,
    double rho,
    // Cached intermediate values
    const ParticleIntermediate& cache
) {
  
  static const arma::mat I2 = arma::eye<arma::mat>(2, 2);
  
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
      // For argos, still needed to compute proposals (different random draw in weight computation)
      List gaussian_proposal1 = product_gaussian_cpp(cache.invQ, invS_argos1,
                                                     cache.mean, y, M, proposal_weight);
      List gaussian_proposal2 = product_gaussian_cpp(cache.invQ, invS_argos2,
                                                     cache.mean, y, M, proposal_weight);
      global_timer.record("weight_proposal");
      
      List error_params_local = List::create(
        Named("sigma_obs") = sigma_obs,
        Named("df") = df,
        Named("rho") = rho,
        Named("a") = a,
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
        Named("rho") = rho,
        Named("a") = a,
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
//' @param print_timing Logical, whether to print profiling timing results (default: FALSE)
//' @param obs_error_params Optional list of length N, where each element is an error_params
//'   list for the corresponding observation. When provided (non-NULL), overrides error_params
//'   at each time step. Useful for ARGOS data where location class varies per observation.
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
    bool verbose,
    bool print_timing,
    Nullable<List> obs_error_params
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
  
  // Unpack obs_error_params if provided (for time-varying error params, e.g. ARGOS classes)
  bool has_obs_error_params = obs_error_params.isNotNull();
  List obs_ep_list;
  if (has_obs_error_params) {
    obs_ep_list = as<List>(obs_error_params);
    if ((int)obs_ep_list.size() != N) {
      Rcpp::stop("obs_error_params must have length equal to the number of observations (%d)", N);
    }
  }
  
  // Main particle filter loop
  for (int j = 0; j < N - 1; ++j) {
    
    if (verbose) {
      Rcout << "Time step: " << j + 1 << " / " << N - 1 << "\n";
    }
    
    double delta = deltas(j);
    arma::vec y = observations.row(j + 1).subvec(1, 2).t();
    
    // --- Get error parameters for this time step ---
    // Use per-observation params if supplied, otherwise fall back to global precomputed values
    double step_sigma_obs   = sigma_obs;
    double step_scale       = scale;
    double step_df          = df;
    double step_rho         = rho;
    double step_a           = a;
    double step_log_p       = log_p_argos;
    double step_log_1mp     = log_1mp_argos;
    arma::mat step_invS1    = invS_argos1;
    arma::mat step_invS2    = invS_argos2;
    List step_error_params  = error_params;   // used by propagate_particle_with_cache
    
    if (has_obs_error_params && error_dist == "argos") {
      // obs_error_params is indexed by observation (j+1 is the obs used at step j)
      step_error_params = as<List>(obs_ep_list[j + 1]);
      step_df          = as<double>(step_error_params["df"]);
      step_sigma_obs   = as<double>(step_error_params["sigma_obs"]);
      step_rho         = as<double>(step_error_params["rho"]);
      step_a           = as<double>(step_error_params["a"]);
      double step_p    = as<double>(step_error_params["p"]);
      step_log_p       = std::log(step_p);
      step_log_1mp     = std::log(1.0 - step_p);
      
      arma::mat Sigma1 = step_sigma_obs*step_sigma_obs *
                         arma::mat{{1, step_rho*std::sqrt(step_a)},
                                   {step_rho*std::sqrt(step_a), 1}};
      arma::mat Sigma2 = step_sigma_obs*step_sigma_obs *
                         arma::mat{{1, -step_rho*std::sqrt(step_a)},
                                   {-step_rho*std::sqrt(step_a), 1}};
      step_invS1 = arma::inv(step_df/(step_df-2.0) * Sigma1);
      step_invS2 = arma::inv(step_df/(step_df-2.0) * Sigma2);
    }
    
    // --- PREDICTION STEP ---
    if (verbose) Rcout << "  Prediction step...\n";
    
    // Extract current particle positions
    arma::mat X_positions = particles.slice(j).cols(0, 1);
    
    // Compute push for all particles (vectorized)
    arma::mat push_matrix = compute_push_matrix_cpp(X_positions, polygon_coords, lambda);
    push_array.slice(j) = push_matrix;
    
    global_timer.record("compute_push");
    
    // Choose centers for all particles if splitting around a fixed point
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
    
 
    // Propagate each particle and cache intermediate computations
    for (int k = 0; k < num_particles; ++k) {
      
      arma::vec U_prev = particles.slice(j).row(k).t();
      arma::vec push = push_matrix.row(k).t();
      
      // Get fixed point index for this particle
      int ind_fp = (split_around_fixed_point && ind_fixed_point_vec(k) > 0) 
                ? ind_fixed_point_vec(k) : 0;
      
      // Propagate particle 
      arma::vec U_next = propagate_particle_with_cache(
        U_prev, y, M, delta, push,
        potential_params, tau, nu, omega, lambda,
        error_dist, step_error_params, scheme, polygon_coords,
        ind_fp, use_precomputed_LQ, L_precomputed, Q_precomputed,
        proposal_weight,
        particle_cache[k]);
      
      particles.slice(j + 1).row(k) = U_next.t();
    }
    
    global_timer.record("propagate_particles");
    
    if (verbose) Rcout << "  Prediction step complete.\n";
  
   // --- CORRECTION STEP ---
    if (verbose) Rcout << "  Correction step...\n";
    
    // Compute weights using cached intermediate values.
    // weights(k, j+1) = weights(k, j) * incremental(k): the entering weight
    // weights(k, j) is already normalized (sums to 1), and is only uniform
    // (1/num_particles) when step j resampled. When resampling was skipped
    // (adaptive ESS_threshold < 1), weights(k, j) carries the actual
    // non-uniform prior weight, which must be propagated forward 
    for (int k = 0; k < num_particles; ++k) {

      arma::vec U_pred = particles.slice(j + 1).row(k).t();

      // Use cached values to avoid recomputation of ODE, SDE, Cholesky!
      double incremental = compute_weight_from_cache(
        U_pred, y, M,
        error_dist, scheme,
        proposal_weight, step_sigma_obs, step_scale, step_df,
        step_invS1, step_invS2, step_log_p, step_log_1mp, step_a, step_rho,
        particle_cache[k]); 

      weights(k, j + 1) = weights(k, j) * incremental;
     }

    global_timer.record("compute_weights");

    // Store total weight
    total_weights(j + 1) = arma::sum(weights.col(j + 1));

    // Update log-likelihood.
    loglik += std::log(total_weights(j + 1));
    loglik_vector(j) = std::log(total_weights(j + 1));
    
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
        // Resample using systematic resampling
        resampled_at(j) = true;
        if (verbose) Rcout << "  Resampling (ESS ratio = " << ess_ratio << ")\n";
        
        // Systematic resampling 
        arma::vec cum_weights = arma::cumsum(weights.col(j + 1));
        double u0 = R::runif(0.0, 1.0 / num_particles);
        
        arma::mat particles_resampled = particles.slice(j + 1);
        int idx = 0;
        
        for (int k = 0; k < num_particles; ++k) {
          double u = u0 + k / double(num_particles);
          
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
          ancestors(k, j) = k + 1;  // No resampling
        }
        if (verbose) Rcout << "  No resampling (ESS ratio = " << ess_ratio << ")\n";
      }
    }
  }  
  
  if (verbose) Rcout << "Particle filtering complete.\n";
  
  // Print timing results if requested
  if (print_timing) {
    global_timer.print_timings();
  }
  
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


// ==================== Conditional Particle Filter with Ancestor Sampling (CPF-AS) ====================
// Implements the CPF-AS kernel (Svensson et al. 2015; Lindsten, Jordan & Schon 2014) as an
// MCMC (Particle Gibbs) alternative to forward-filtering backward-sampling for drawing
// trajectories from the smoothing distribution.

// Log transition density p(U_target | U_prev) under the numerical splitting scheme.
// Mirrors the backward-weight computation used in forward_filtering_backward_sampling() (R):
// for Strang, the final gradient term always uses the naive (non-fixed-point) mixture
// gradient, matching that function's documented limitation.
static inline double transition_log_density(
    const arma::vec& U_prev,
    const arma::vec& U_target,
    double delta,
    const arma::vec& push,
    const List& potential_params,
    double tau,
    double nu,
    double omega,
    double lambda,
    const std::string& scheme,
    const arma::mat& polygon_coords,
    int ind_fixed_point,
    bool use_precomputed_LQ,
    const arma::mat& L_precomputed,
    const arma::mat& Q_precomputed
) {
  Nullable<int> ind_fp_arg = (ind_fixed_point > 0) ? Nullable<int>(wrap(ind_fixed_point)) : Nullable<int>(R_NilValue);

  if (scheme == "Lie-Trotter") {
    arma::vec U_hat = solve_ODE_cpp(U_prev, delta, push, potential_params, ind_fp_arg);

    List OU_solution = use_precomputed_LQ
      ? solve_SDE_cpp(U_hat, delta, tau, nu, omega, potential_params, ind_fp_arg, wrap(L_precomputed), wrap(Q_precomputed))
      : solve_SDE_cpp(U_hat, delta, tau, nu, omega, potential_params, ind_fp_arg, R_NilValue, R_NilValue);

    arma::mat Q = as<arma::mat>(OU_solution["Q"]);
    arma::vec mean = as<arma::vec>(OU_solution["mean"]);

    return log_dmvnorm_chol_cpp(U_target, mean, chol_cpp(Q));

  } else { // Strang
    arma::mat x_star = potential_params["x_star"];

    arma::vec U_hat = solve_ODE_cpp(U_prev, delta / 2.0, push, potential_params, ind_fp_arg);

    List OU_solution = use_precomputed_LQ
      ? solve_SDE_cpp(U_hat, delta, tau, nu, omega, potential_params, ind_fp_arg, wrap(L_precomputed), wrap(Q_precomputed))
      : solve_SDE_cpp(U_hat, delta, tau, nu, omega, potential_params, ind_fp_arg, R_NilValue, R_NilValue);

    arma::mat Q = as<arma::mat>(OU_solution["Q"]);
    arma::vec mean = as<arma::vec>(OU_solution["mean"]);

    arma::vec X_next = U_target.subvec(0, 1);
    arma::vec V_next = U_target.subvec(2, 3);

    arma::vec push_next = compute_push_cpp(X_next, polygon_coords, lambda);
    arma::vec grad_next = mix_gaussian_grad_cpp(X_next, x_star, potential_params, IntegerVector());

    arma::vec V_tilde = V_next + (delta / 2.0) * (push_next + grad_next);
    arma::vec U_tilde_next = arma::join_vert(X_next, V_tilde);

    return log_dmvnorm_chol_cpp(U_tilde_next, mean, chol_cpp(Q));
  }
}

// Run a single CPF-AS sweep conditional on a reference trajectory, returning a new sampled trajectory.
static arma::mat run_cpf_as_sweep(
    const arma::mat& observations,
    const List& sde_params,
    const List& potential_params,
    const List& error_params,
    const std::string& error_dist,
    const arma::mat& polygon_coords,
    const arma::vec& U0,
    double lambda,
    int K,
    const std::string& scheme,
    bool split_around_fixed_point,
    double proposal_weight,
    double ESS_threshold,
    const arma::mat& ref_traj,
    bool has_obs_error_params,
    const List& obs_ep_list,
    arma::vec& ess_history_out,
    LogicalVector& resampled_at_out
) {
  int N = observations.n_rows;
  int num_states = 4;

  double tau = as<double>(sde_params["tau"]);
  double nu = as<double>(sde_params["nu"]);
  double omega = as<double>(sde_params["omega"]);

  arma::mat I2 = arma::eye<arma::mat>(2, 2);
  arma::mat M = arma::join_horiz(I2, arma::zeros<arma::mat>(2, 2));

  arma::mat x_star = potential_params["x_star"];

  arma::vec times = observations.col(0);
  arma::vec deltas = arma::diff(times);

  bool use_precomputed_LQ = false;
  arma::mat L_precomputed, Q_precomputed;
  if (arma::max(deltas) - arma::min(deltas) < 1e-6 && !split_around_fixed_point) {
    use_precomputed_LQ = true;
    double delta0 = deltas(0);
    L_precomputed = RACVM_link_cpp(tau, omega, delta0);
    Q_precomputed = RACVM_cov_cpp(tau, nu, omega, delta0);
  }

  double sigma_obs = 0.0, scale = 0.0, df = 0.0, rho = 0.0, a = 0.0, p_argos = 0.0;
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

    arma::mat Sigma1 = sigma_obs*sigma_obs * arma::mat{{1, rho*std::sqrt(a)}, {rho*std::sqrt(a), 1}};
    arma::mat Sigma2 = sigma_obs*sigma_obs * arma::mat{{1, -rho*std::sqrt(a)}, {-rho*std::sqrt(a), 1}};
    invS_argos1 = arma::inv(df/(df-2.0) * Sigma1);
    invS_argos2 = arma::inv(df/(df-2.0) * Sigma2);
    log_p_argos = std::log(p_argos);
    log_1mp_argos = std::log(1.0 - p_argos);
  }

  arma::cube particles(K, num_states, N);
  arma::mat weights(K, N);
  arma::imat ancestors(K, N - 1);

  arma::mat R0 = arma::diagmat(arma::vec{0.01, 0.01, 0.25, 0.25});
  arma::mat sqrtR0 = arma::sqrtmat_sympd(R0);
  for (int k = 0; k < K - 1; ++k) {
    arma::vec noise = sqrtR0 * arma::randn<arma::vec>(num_states);
    particles.slice(0).row(k) = (U0 + noise).t();
  }
  particles.slice(0).row(K - 1) = ref_traj.row(0);
  weights.col(0).fill(1.0 / K);
  ess_history_out(0) = K;

  for (int j = 0; j < N - 1; ++j) {

    double delta = deltas(j);
    arma::vec y = observations.row(j + 1).subvec(1, 2).t();

    double step_sigma_obs = sigma_obs, step_scale = scale, step_df = df, step_rho = rho, step_a = a;
    double step_log_p = log_p_argos, step_log_1mp = log_1mp_argos;
    arma::mat step_invS1 = invS_argos1, step_invS2 = invS_argos2;
    List step_error_params = error_params;

    if (has_obs_error_params && error_dist == "argos") {
      step_error_params = as<List>(obs_ep_list[j + 1]);
      step_df = as<double>(step_error_params["df"]);
      step_sigma_obs = as<double>(step_error_params["sigma_obs"]);
      step_rho = as<double>(step_error_params["rho"]);
      step_a = as<double>(step_error_params["a"]);
      double step_p = as<double>(step_error_params["p"]);
      step_log_p = std::log(step_p);
      step_log_1mp = std::log(1.0 - step_p);

      arma::mat Sigma1 = step_sigma_obs*step_sigma_obs * arma::mat{{1, step_rho*std::sqrt(step_a)}, {step_rho*std::sqrt(step_a), 1}};
      arma::mat Sigma2 = step_sigma_obs*step_sigma_obs * arma::mat{{1, -step_rho*std::sqrt(step_a)}, {-step_rho*std::sqrt(step_a), 1}};
      step_invS1 = arma::inv(step_df/(step_df-2.0) * Sigma1);
      step_invS2 = arma::inv(step_df/(step_df-2.0) * Sigma2);
    }

    arma::mat X_positions = particles.slice(j).cols(0, 1);
    arma::mat push_matrix = compute_push_matrix_cpp(X_positions, polygon_coords, lambda);

    arma::ivec ind_fp_vec;
    if (split_around_fixed_point) {
      ind_fp_vec = choose_center_matrix_cpp(X_positions, x_star, potential_params);
    }

    arma::vec w_prev = weights.col(j);
    arma::vec cw_prev = arma::cumsum(w_prev);

    double ess = 1.0 / arma::sum(arma::square(w_prev));
    bool do_resample = (ess / K) < ESS_threshold;
    resampled_at_out(j) = do_resample;

    // --- Resample ancestors for particles 0..K-2 (systematic, adaptive on ESS) ---
    arma::ivec anc(K);
    if (do_resample) {
      double u0 = R::runif(0.0, 1.0 / (K - 1));
      int idx = 0;
      for (int k = 0; k < K - 1; ++k) {
        double u = u0 + k / double(K - 1);
        while (idx < K - 1 && cw_prev(idx) < u) idx++;
        anc(k) = idx;
      }
    } else {
      for (int k = 0; k < K - 1; ++k) anc(k) = k;
    }

    std::vector<ParticleIntermediate> cache(K);

    // --- Propagate particles 0..K-2 via the proposal ---
    for (int k = 0; k < K - 1; ++k) {
      int anc_idx = anc(k);
      arma::vec U_prev = particles.slice(j).row(anc_idx).t();
      arma::vec push = push_matrix.row(anc_idx).t();
      int ind_fp = (split_around_fixed_point && ind_fp_vec(anc_idx) > 0) ? ind_fp_vec(anc_idx) : 0;

      arma::vec U_next = propagate_particle_with_cache(
        U_prev, y, M, delta, push, potential_params,
        tau, nu, omega, lambda, error_dist, step_error_params, scheme, polygon_coords,
        ind_fp, use_precomputed_LQ, L_precomputed, Q_precomputed, proposal_weight, cache[k]);

      particles.slice(j + 1).row(k) = U_next.t();
    }

    // --- Pin particle K-1 to the reference trajectory ---
    particles.slice(j + 1).row(K - 1) = ref_traj.row(j + 1);
    arma::vec U_ref_next = ref_traj.row(j + 1).t();

    // --- Ancestor sampling for the reference particle ---
    arma::vec log_as_weights(K);
    for (int m = 0; m < K; ++m) {
      arma::vec U_prev_m = particles.slice(j).row(m).t();
      arma::vec push_m = push_matrix.row(m).t();
      int ind_fp_m = (split_around_fixed_point && ind_fp_vec(m) > 0) ? ind_fp_vec(m) : 0;

      double log_trans = transition_log_density(
        U_prev_m, U_ref_next, delta, push_m, potential_params,
        tau, nu, omega, lambda, scheme, polygon_coords,
        ind_fp_m, use_precomputed_LQ, L_precomputed, Q_precomputed);

      log_as_weights(m) = std::log(w_prev(m)) + log_trans;
    }
    double mx = log_as_weights.max();
    arma::vec as_w = arma::exp(log_as_weights - mx);
    as_w /= arma::sum(as_w);
    arma::vec cw_as = arma::cumsum(as_w);
    double u_as = R::runif(0.0, 1.0);
    int a_ref = 0;
    while (a_ref < K - 1 && cw_as(a_ref) < u_as) a_ref++;
    anc(K - 1) = a_ref;

    // Cache needed to compute the reference particle's importance weight
    {
      arma::vec U_prev = particles.slice(j).row(a_ref).t();
      arma::vec push = push_matrix.row(a_ref).t();
      int ind_fp = (split_around_fixed_point && ind_fp_vec(a_ref) > 0) ? ind_fp_vec(a_ref) : 0;
      ParticleIntermediate ref_cache;
      propagate_particle_with_cache(
        U_prev, y, M, delta, push, potential_params,
        tau, nu, omega, lambda, error_dist, step_error_params, scheme, polygon_coords,
        ind_fp, use_precomputed_LQ, L_precomputed, Q_precomputed, proposal_weight, ref_cache);
      cache[K - 1] = ref_cache;
    }

    // --- Compute importance weights for all K particles ---
    // When particles 0..K-2 were not resampled, their previous (non-uniform) weight must
    // carry forward multiplicatively; when resampled, the incremental factor alone suffices
    // since the resampling step already incorporated the previous weight into the selection
    // probability. The reference particle (K-1) is always "resampled" via ancestor sampling,
    // so its weight is always incremental only.
    for (int k = 0; k < K; ++k) {
      arma::vec U_pred = particles.slice(j + 1).row(k).t();
      double incremental = compute_weight_from_cache(
        U_pred, y, M, error_dist, scheme, proposal_weight,
        step_sigma_obs, step_scale, step_df,
        step_invS1, step_invS2, step_log_p, step_log_1mp, step_a, step_rho,
        cache[k]);
      weights(k, j + 1) = (k < K - 1 && !do_resample) ? w_prev(anc(k)) * incremental : incremental;
    }

    ancestors.col(j) = anc;

    double wsum = arma::sum(weights.col(j + 1));
    weights.col(j + 1) /= wsum;
    ess_history_out(j + 1) = 1.0 / arma::sum(arma::square(weights.col(j + 1)));
  }

  // --- Final draw and genealogy tracing ---
  arma::vec wN = weights.col(N - 1);
  arma::vec cwN = arma::cumsum(wN);
  double u_final = R::runif(0.0, 1.0);
  int idx = 0;
  while (idx < K - 1 && cwN(idx) < u_final) idx++;

  arma::mat traj(N, num_states);
  traj.row(N - 1) = particles.slice(N - 1).row(idx);
  for (int t = N - 2; t >= 0; --t) {
    idx = ancestors(idx, t);
    traj.row(t) = particles.slice(t).row(idx);
  }

  return traj;
}

//' Conditional particle filter with ancestor sampling (CPF-AS)
//'
//' Implements the CPF-AS kernel (Svensson et al. 2015; Lindsten, Jordan and Schon 2014)
//' as an MCMC (Particle Gibbs) alternative to forward-filtering backward-sampling for
//' drawing trajectories from the smoothing distribution. One or more CPF-AS sweeps are
//' run, each conditioning on the trajectory sampled by the previous sweep (or on
//' \code{reference_trajectory} for the first sweep), and the final sampled trajectory is
//' returned. Calling this repeatedly across SGD iterations, each time passing back in the
//' previously returned trajectory as \code{reference_trajectory}, forms a single persistent
//' MCMC chain over trajectories (Markovian stochastic approximation).
//'
//' @param observations Matrix with dimensions N x 3: time, Y1, Y2
//' @param sde_params List with elements: tau, nu, omega
//' @param potential_params List with elements: alpha, B, x_star
//' @param error_params List of parameters for the measurement error distribution
//' @param error_dist String specifying error distribution: "normal", "scaled_t", or "argos"
//' @param polygon_coords Matrix of polygon boundary coordinates (N x 2)
//' @param U0 Numeric vector of initial state (length 4)
//' @param lambda Numeric, penalization parameter
//' @param num_particles Integer, number of particles K (K-1 propagated + 1 reference)
//' @param scheme String, splitting scheme: "Lie-Trotter" or "Strang"
//' @param split_around_fixed_point Logical, whether to split around fixed points
//' @param proposal_weight Numeric between 0 and 1, weight for the Gaussian proposal
//' @param reference_trajectory Matrix (N x 4) giving the initial conditioning trajectory
//' @param n_sweeps Integer, number of CPF-AS sweeps to run (default 1)
//' @param ESS_threshold Numeric between 0 and 1, threshold on effective sample size (as a
//'   fraction of num_particles) below which particles 0..K-2 are resampled (systematic
//'   resampling) at a given step. The reference particle's ancestor is always resampled via
//'   ancestor sampling regardless of this threshold. Use 1 to resample at every step.
//' @param obs_error_params Optional list of length N of per-observation error_params (ARGOS)
//'
//' @return List with elements:
//'   \itemize{
//'     \item trajectory: matrix (N x 4) sampled from the (approximate) smoothing distribution
//'     \item ess_history: vector of ESS at each time step of the last sweep
//'     \item resampled_at: logical vector indicating steps where particles 0..K-2 were
//'       resampled in the last sweep
//'   }
//' @export
// [[Rcpp::export]]
List conditional_particle_filter_cpp(
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
    double proposal_weight,
    const arma::mat& reference_trajectory,
    int n_sweeps,
    double ESS_threshold,
    Nullable<List> obs_error_params
) {
  int N = observations.n_rows;

  bool has_obs_error_params = obs_error_params.isNotNull();
  List obs_ep_list;
  if (has_obs_error_params) {
    obs_ep_list = as<List>(obs_error_params);
    if ((int)obs_ep_list.size() != N) {
      Rcpp::stop("obs_error_params must have length equal to the number of observations (%d)", N);
    }
  }

  arma::mat ref = reference_trajectory;
  arma::vec ess_history(N);
  LogicalVector resampled_at(N - 1);

  for (int s = 0; s < n_sweeps; ++s) {
    ref = run_cpf_as_sweep(
      observations, sde_params, potential_params, error_params, error_dist,
      polygon_coords, U0, lambda, num_particles, scheme, split_around_fixed_point,
      proposal_weight, ESS_threshold, ref, has_obs_error_params, obs_ep_list,
      ess_history, resampled_at);
  }

  return List::create(
    Named("trajectory") = ref,
    Named("ess_history") = ess_history,
    Named("resampled_at") = resampled_at
  );
}


