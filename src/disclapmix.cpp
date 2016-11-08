#include <iostream>
#include <vector>
#include <ctime>

#include "disclapmix.h"

using namespace Eigen;

void check_input(const MatrixXi& x, MatrixXi& y, 
  VectorXd& beta, VectorXd& tau, MatrixXd& v_mat,
  size_t verbose_level,
  size_t maxit_centers_change, size_t maxit_disclap_em, size_t maxit_irls, 
  double eps_v_mat,
  bool use_deviance, double eps_deviance, 
  bool use_beta, double eps_beta) {

  if (verbose_level < 0) {
    throw "verbose_level must be >= 0";
  }

  if (maxit_centers_change <= 0) {
    throw "maxit_centers_change must be > 0";
  }
  
  if (maxit_disclap_em <= 0) {
    throw "maxit_disclap_em must be > 0";
  }
  
  if (maxit_irls <= 0) {
    throw "maxit_irls must be > 0";
  }
  
  if (eps_v_mat <= 0) {
    throw "eps_v_mat must be > 0";
  }
  
  // maxit_disclap_em already checked (ensures correct name in output)
  check_common_input(maxit_disclap_em, use_deviance, eps_deviance, use_beta, eps_beta);
  
  check_x_y(x, y);
  check_beta(beta, x, y);
  check_tau(tau, y);
  check_v_mat(v_mat, x, y);
}

MatrixXd calculate_non_norm_v_mat(const MatrixXi& x, MatrixXi& y, const MatrixXd& disclap_ps, const VectorXd& tau) {
  check_x_y(x, y);

  size_t individuals = x.rows();
  size_t clusters = y.rows();
  size_t loci = x.cols();
  
  MatrixXd wij(individuals, clusters);
  
  for (size_t j = 0; j < clusters; j++) {
    for (size_t i = 0; i < individuals; i++) {
      double prod = tau[j];
      
      for (size_t k = 0; k < loci; k++) {
        double p = disclap_ps(j, k);
        int d = std::abs(x(i, k) - y(j, k));
        prod *= std::pow(p, d)*((1-p)/(1+p));
      }
      
      wij(i, j) = prod;
    }
  }

  return wij;
}

MatrixXd calculate_norm_v_mat(const MatrixXd& wij) {
/*  if (wij == NULL) {
    throw "wij is NULL";
  }
  */
  size_t individuals = wij.rows();
  size_t clusters = wij.cols();
  
  MatrixXd vic(individuals, clusters);
  VectorXd sums = wij.rowwise().sum();
  
  for (size_t i = 0; i < individuals; i++) {
    for (size_t j = 0; j < clusters; j++) {
      vic(i, j) = wij(i, j) / sums[i];
    }
  }
  
  return(vic);
}

void update_weight_vec(std::vector<MatrixXd>& weight_vec, MatrixXd& v) {
/*  if (weight_vec == NULL) {
    throw "w is NULL";
  }
  
  if (v == NULL) {
    throw "w is NULL";
  }
  */
  size_t individuals = weight_vec.size();
  size_t clusters = weight_vec.at(0).rows();
  size_t loci = weight_vec.at(0).cols();
  
  //OUTPUT << "=====" << std::endl; 
  //OUTPUT << *v << std::endl;
  //OUTPUT << "=====" << std::endl;
  
  for (size_t i = 0; i < individuals; i++) {
    //OUTPUT << weight_vec[i] << std::endl;
    //OUTPUT << "----" << std::endl;
    
    for (size_t j = 0; j < clusters; j++) {
      for (size_t k = 0; k < loci; k++) {
        weight_vec[i](j, k) = v(i, j);
      }
    }
  }
  /*
  OUTPUT << "=====" << std::endl;
   
  for (size_t i = 0; i < individuals; i++) {
    OUTPUT << weight_vec[i] << std::endl;
    OUTPUT << "----" << std::endl;
  }
  
  OUTPUT << "=====" << std::endl;
  */
}

bool change_centers(const MatrixXi& x, MatrixXi& y, VectorXi& x_mins, VectorXi& x_maxs, MatrixXd& v) {
  check_x_y(x, y);
  
  // FIXME: Check if some centers are duplicate?!
  
/*  if (x_mins == NULL) {
    throw "x_mins is NULL";
  }
  
  if (x_maxs == NULL) {
    throw "x_maxs is NULL";
  }
  
  if (v == NULL) {
    throw "w is NULL";
  }
 */
  size_t individuals = x.rows();
  size_t clusters = y.rows();
  size_t loci = x.cols();
  
  MatrixXi new_y(clusters, loci);
  
  if (x_mins.size() != loci) {
    throw "length(x_mins) != loci";
  }

  if (x_maxs.size() != loci) {
    throw "length(x_maxs) != loci";
  }
  
  for (size_t k = 0; k < loci; k++) {
    if (x_mins[k] >= x_maxs[k]) {
      throw "x_mins[k] >= x_maxs[k]";
    }
  }
  
  for (size_t k = 0; k < loci; k++) {    
    for (size_t j = 0; j < clusters; j++) {
      if (x_mins[k] == x_maxs[k]) {
        new_y(j, k) = y(j, k);
        continue;
      }
      
      int y_best = -1;
      bool y_best_found = false;
      double logl = std::numeric_limits<double>::infinity();
      
      for (int y_cand = x_mins[k]; y_cand <= x_maxs[k]; y_cand++) {
        double y_cand_log_l = 0.0;
        
        //OUTPUT << "y_cand = " << y_cand << std::endl;
        
        for (size_t i = 0; i < individuals; i++) {
          double v_ij = v(i, j);
          y_cand_log_l += v_ij * std::abs(x(i, k) - y_cand);
          //OUTPUT << "  v_ij * |x - y| = " << v_ij << " * | " << std::abs(x(i, k) - y_cand) << "|" << std::endl;
        }
        
        if (y_cand_log_l < logl) {
          logl = y_cand_log_l;
          y_best = y_cand;
          y_best_found = true;
        }
        
        //OUTPUT << "logl = " << logl << ", y_cand_log_l = " << y_cand_log_l << std::endl;
      }
      
      if (!y_best_found) {
        //OUTPUT << "x min = " << x_mins[k] << std::endl;
        //OUTPUT << "x max = " << x_maxs[k] << std::endl;
        //FLUSH_OUTPUT;
        
        throw "Unexpected error when finding new center";
      }
      
      new_y(j, k) = y_best;
    }
  }
  
  bool changed = false;

  for (size_t j = 0; j < clusters; j++) {  
    for (size_t k = 0; k < loci; k++) {    
      if (y(j, k) != new_y(j, k)) {
        changed = true;
        break;
      }
    }
  }
  
  if (changed) {
    y = new_y;
    return true;
  }
  
  return false;
}

bool disclapmix(const MatrixXi& x, MatrixXi& y, 
                VectorXd& beta, VectorXd& tau, MatrixXd& v_mat, 
  MatrixXd& vcov, 
  double& deviance,
  std::vector<IterationsInfo>& iter_info,
  size_t& model_observations, 
  size_t& model_parameters, 
  size_t& iter_irls,
  size_t& iter_em,
  size_t& iter_centers,
  size_t verbose_level, 
  size_t maxit_centers_change, size_t maxit_disclap_em, size_t maxit_irls,
  double eps_v_mat,
  bool use_deviance, double eps_deviance, 
  bool use_beta, double eps_beta) {
  
  clock_t clock_begin = clock();
  
  check_input(x, y, beta, tau, v_mat, verbose_level,
    maxit_centers_change, maxit_disclap_em, maxit_irls,
    eps_v_mat, use_deviance, eps_deviance, use_beta, eps_beta);

  size_t individuals = x.rows();
  size_t clusters = y.rows();
  size_t loci = x.cols();
  
  model_observations = individuals*loci;
  //                  [central haplotypes]  [beta . disclap pars]  [tau]
  model_parameters = (clusters*loci) +     (clusters + loci - 1) + (clusters - 1);
  
  if (iter_info.size() != 0) {
    throw "iter_info must be initialized and empty";
  }
  
  IterationsInfo iterinfo = IterationsInfo();
  iterinfo.y = y;
  iter_info.push_back(iterinfo);
  
  VectorXd lin_pred(individuals * clusters * loci);
  MatrixXd disclap_ps;
  double logL_full;
  double logL_marginal;

  VectorXi x_mins = x.colwise().minCoeff();
  VectorXi x_maxs = x.colwise().maxCoeff();
  x_mins = x_mins.array() - 1;
  x_maxs = x_maxs.array() + 1;

  if (verbose_level >= 2) {
    OUTPUT << "Input data characteristics:" << std::endl;
    OUTPUT << "  Individuals          = " << individuals << std::endl;
    OUTPUT << "  Clusters             = " << clusters << std::endl;
    OUTPUT << "  Loci                 = " << loci << std::endl;
    OUTPUT << std::endl;
    
    OUTPUT << "Running with the following parameters:" << std::endl;
    OUTPUT << "  verbose_level        = " << verbose_level << std::endl;
    OUTPUT << "  maxit_centers_change = " << maxit_centers_change << std::endl;
    OUTPUT << "  maxit_disclap_em     = " << maxit_disclap_em << std::endl;
    OUTPUT << "  maxit_irls           = " << maxit_irls << std::endl;
    OUTPUT << "  eps_v_mat            = " << eps_v_mat << std::endl;
    OUTPUT << "  use_deviance         = " << use_deviance << std::endl;
    
    if (use_deviance) {
      OUTPUT << "  eps_deviance         = " << eps_deviance << std::endl;
    }
    
    OUTPUT << "  use_beta             = " << use_beta << std::endl;
    
    if (use_beta) {
      OUTPUT << "  eps_beta             = " << eps_beta << std::endl;
    }
    
    OUTPUT << std::endl;
  }  
    
  std::vector<MatrixXi> response_vec = init_response_vec(x, y);
  std::vector<MatrixXd> weight_vec = init_weight_vec(individuals, clusters, loci);
  
  bool converged = false;
  bool irls_verbose = verbose_level >= 4;
  
  size_t iter_total_irls = 0;
  size_t iter_total_em = 0;
  size_t iter_total_centers = 1; // we start with some centers
    
  for (size_t iter_outer = 0; iter_outer < maxit_centers_change; iter_outer++) {    
    if (verbose_level == 1) {
      OUTPUT << "[" << (iter_outer + 1) << "]: ";
    }

    if (verbose_level >= 2) {
      OUTPUT << "================================================================================" << std::endl;
      OUTPUT << "Outer (change central haplotypes, y) iteration " << (iter_outer + 1) << std::endl;
      OUTPUT << "================================================================================" << std::endl;
    }

    if (verbose_level >= 2) {
      OUTPUT << "  Initial tau  = " << tau.transpose().format(disclapEigenFmt) << std::endl;
      OUTPUT << "  Initial beta = " << beta.transpose().format(disclapEigenFmt) << std::endl;
    }
          
    for (size_t iter_inner = 0; iter_inner < maxit_disclap_em; iter_inner++) {
      if (verbose_level == 1) {
        OUTPUT << "(" << (iter_inner + 1) << ")";
      }
      
      if (verbose_level >= 3) {
        OUTPUT << "  ------------------------------------------------------------------------------" << std::endl;
        OUTPUT << "  Inner (EM) iteration " << (iter_inner + 1) << std::endl;
        OUTPUT << "  ------------------------------------------------------------------------------" << std::endl;
      }      
      
      size_t irls_iterations = 0;
            
      converged = irls(response_vec, weight_vec, beta, lin_pred, vcov, deviance, irls_iterations, 
        irls_verbose, maxit_irls, use_deviance, eps_deviance, use_beta, eps_beta, false);
      
      iter_total_irls += irls_iterations;
      iter_total_em += 1;
      
      if (verbose_level >= 3) {
        if (converged) {
          OUTPUT << "    IRLS did converge in " << irls_iterations << " iterations" << std::endl;
        } else {
          OUTPUT << "    !!! IRLS did NOT converge (in " << irls_iterations << " iterations" << ") !!!" << std::endl;
        }
        
        OUTPUT << "    Deviance          = " << deviance << std::endl;
        OUTPUT << "    Coefficients      = " << beta.transpose().format(disclapEigenFmt) << std::endl;
        //OUTPUT << "    Coefficients ext  = " << get_beta_extended(beta).transpose() << std::endl;
      }

      disclap_ps = beta_to_discrete_laplace_parameters(x, y, beta); 
      const MatrixXd wij = calculate_non_norm_v_mat(x, y, disclap_ps, tau);    
      MatrixXd new_v_mat = calculate_norm_v_mat(wij);
      update_weight_vec(weight_vec, new_v_mat);
 
      double impr = (new_v_mat.array() - v_mat.array()).cwiseAbs().maxCoeff() / v_mat.maxCoeff();

      v_mat = new_v_mat;
      VectorXd new_tau = v_mat.colwise().sum() / individuals;
      tau = new_tau;

      logL_full = loglikelihood_full(x, y, v_mat, disclap_ps, tau);
      logL_marginal = loglikelihood_marginal(x, y, disclap_ps, tau);

      iter_info.back().beta_vec.push_back(beta);
      iter_info.back().deviance_vec.push_back(deviance);
      iter_info.back().logL_full_vec.push_back(logL_full);
      iter_info.back().logL_marginal_vec.push_back(logL_marginal);
            
      if (isnan(logL_full)) {
        throw "logL full is NA, maybe IWLS did not converge; try another init_y";
      }
            
      if (verbose_level >= 3) {
        if (impr < eps_v_mat) {
          OUTPUT << "    v_mat improvement = " << impr << " < " << eps_v_mat << " [CONVERGENCE]" << std::endl;
        } else {
          OUTPUT << "    v_mat improvement = " << impr << " > " << eps_v_mat << std::endl;
        }

        OUTPUT << "    Tau               = " << tau.transpose().format(disclapEigenFmt) << std::endl;
        OUTPUT << "    log(L) full       = " << logL_full << std::endl;
        OUTPUT << "    log(L) marginal   = " << logL_marginal << std::endl;

        /*
        OUTPUT << "improvement = " << impr << std::endl;
        OUTPUT << "    v_ij matrix:" << std::endl;
        OUTPUT << v_mat << std::endl;
        OUTPUT << "    Resulting tau = " << tau.transpose() << std::endl;
        */
      }
      
      converged = converged && impr < eps_v_mat;
      
      if (verbose_level == 2) {
        OUTPUT << "    Inner (EM) iteration " << (iter_inner + 1) << ": ";
        
        if (converged) {
          OUTPUT << "EM did converge" << std::endl;
        } else {
          OUTPUT << "EM did *NOT* converge" << std::endl;
        }
      } else if (verbose_level >= 3) {
        if (converged) {
          OUTPUT << "  EM did converge (on both IRLS and v_mat improvement)" << std::endl;
        } else {
          OUTPUT << "  EM did *NOT* converge in this iteration (either due to IRLS or v_mat improvement)" << std::endl;
        }
      }
      
      if (converged) {
        break;
      }
    }
    
    if (verbose_level >= 2) {
      OUTPUT << "  Post tau  = " << tau.transpose().format(disclapEigenFmt) << std::endl;
      OUTPUT << "  Post beta = " << beta.transpose().format(disclapEigenFmt) << std::endl;
    }
    
    // FIXME: Reset v_{ij} when moving centers?
    MatrixXi y_old = y;
    bool y_changed = change_centers(x, y, x_mins, x_maxs, v_mat);
    
    if (y_changed) {
      iter_total_centers += 1;
      
      IterationsInfo iterinfo_newy = IterationsInfo();
      iterinfo_newy.y = y;
      iter_info.push_back(iterinfo_newy);
  
      
      if (verbose_level >= 3) {
        OUTPUT << "Centers changed from" << std::endl << y_old.format(disclapEigenFmt) << std::endl << "to" << std::endl << y.format(disclapEigenFmt) << std::endl;
      } else if (verbose_level >= 2) {
        OUTPUT << "Centers changed." << std::endl;
      }
    } else {
      if (verbose_level >= 2) {
        OUTPUT << "Centers did not change." << std::endl;
      }
      
      break;
    }
    
    if (verbose_level == 1) {
      OUTPUT << std::endl;
    }
  }
  
  if (verbose_level >= 1) {
    clock_t clock_end = clock();
    double elapsed_secs = double(clock_end - clock_begin) / CLOCKS_PER_SEC;
    
    OUTPUT << std::endl;
    OUTPUT << "================================================================================" << std::endl;
    OUTPUT << "Final results" << std::endl;
    OUTPUT << "================================================================================" << std::endl;
    OUTPUT << "Total number of IRLS iterations             = " << iter_total_irls << std::endl;
    OUTPUT << "Total number of EM iterations               = " << iter_total_em << std::endl;
    OUTPUT << "Total number of changing central haplotypes = " << iter_total_centers << std::endl;
    OUTPUT << "Elapsed CPU time                            = " << elapsed_secs << " sec" << std::endl;
    OUTPUT << "Deviance                                    = " << deviance << std::endl;
    OUTPUT << "log(L) full                                 = " << logL_full << std::endl;
    OUTPUT << "log(L) marginal                             = " << logL_marginal << std::endl;
    
    OUTPUT << "Coefficients                                = " << beta.transpose().format(disclapEigenFmt) << std::endl;
    OUTPUT << "                                              ";
    beta_description(x, y, beta);
    OUTPUT << std::endl;

    OUTPUT << "Coefficients variance-covariance matrix assuming that tau and y are known" << std::endl;
    OUTPUT << vcov.format(disclapEigenFmt) << std::endl;

    //OUTPUT << "Coefficients ext                            = " << get_beta_extended(beta).transpose() << std::endl;
    OUTPUT << "Tau                                         = " << tau.transpose().format(disclapEigenFmt) << std::endl;

    MatrixXd discrete_laplace_parameters = beta_to_discrete_laplace_parameters(x, y, beta);
    OUTPUT << "y:" << std::endl;
    OUTPUT << y.format(disclapEigenFmt) << std::endl;

    OUTPUT << "Discrete Laplace parameters (rows: clusters / columns: loci):" << std::endl;
    OUTPUT << discrete_laplace_parameters.format(disclapEigenFmt) << std::endl;
    
    OUTPUT << std::endl;
    
    if (converged) {
      OUTPUT << "Last EM did converge (on both IRLS and v_mat improvement)" << std::endl;
    } else {
      OUTPUT << "Last EM did *NOT* converge (either due to IRLS or v_mat improvement)" << std::endl;
    }
  }

  iter_irls = iter_total_irls;
  iter_em = iter_total_em;
  iter_centers = iter_total_centers;

  /*
  if (verbose_level >= 2) {
    OUTPUT << "y:" << std::endl << y << std::endl;
    OUTPUT << "v_mat:" << std::endl << v_mat << std::endl;
  }
  */
  
  return converged;
}

/*
// Try step-halving
      if (deviance > old_deviance) {
        OUTPUT << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        OUTPUT << "!! Deviance increased!";
        OUTPUT << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        
        OUTPUT << "Start beta       = " << start_beta.transpose() << std::endl;
          
        OUTPUT << "beta that failed = " << beta.transpose() << std::endl;
        OUTPUT << "  -> deviance = " << deviance << std::endl;
        OUTPUT << std::endl;
                  
        for (size_t step_halving_i = 0; step_halving_i < 10; ++step_halving_i) {
          double fraction = (1.0d / pow(2.0d, (double)(step_halving_i + 1)));
          
          beta = (1.0d - fraction)*(start_beta.array()) + fraction*(beta.array());
          //beta = beta.array() + 10 * VectorXd::Random(beta.size()).array();
          
          OUTPUT << "fraction = " << fraction << std::endl;
          OUTPUT << "beta new(" << step_halving_i << ")      = " << beta.transpose() << std::endl;
    
          converged = irls(response_vec, weight_vec, beta, lin_pred, deviance, irls_iterations, 
            irls_verbose, maxit_irls, use_deviance, eps_deviance, use_beta, eps_beta, false);
 
          OUTPUT << "  -> deviance = " << deviance << std::endl;
          OUTPUT << std::endl;
          
          if (deviance < old_deviance) {
            break;          
          }
        }
        
        if (deviance > old_deviance) {
          throw "FAIL";
        }
      }
*/

