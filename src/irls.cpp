#include <iostream>
#include <Eigen/LU>

#include "irls.h"
#include "input_validation.h"
#include "disclap_family.h"
#include "beta.h"

using namespace Eigen;

void check_input(const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w, VectorXd& beta, size_t maxit, bool use_deviance, double eps_deviance, bool use_beta, double eps_beta) {
  check_common_input(maxit, use_deviance, eps_deviance, use_beta, eps_beta);
  check_d_w(d, w);
  check_beta(beta, d, w);
}

void update_lin_pred(VectorXd& lin_pred, const VectorXd& beta_ext, size_t individuals, size_t clusters, size_t loci) {
  size_t idx = 0;
  
  for (size_t k = 0; k < loci; k++) {
    double b_k = beta_ext[clusters + k];
    
    for (size_t j = 0; j < clusters; j++) {
      double b_j = beta_ext[j];
      
      double bj_bk = b_j + b_k;
      
      for (size_t i = 0; i < individuals; i++) {
        lin_pred[idx] = bj_bk;
        idx += 1;      
      }
    }
  }  
}

bool deviance_converged(const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w, 
  bool& deviance_calculated,
  VectorXd& lin_pred, const VectorXd& beta, double& deviance, double& old_deviance, double eps_deviance, 
  size_t individuals, size_t clusters, size_t loci, bool verbose) {

  VectorXd beta_ext = get_beta_extended(beta);
    
  update_lin_pred(lin_pred, beta_ext, individuals, clusters, loci);
  
  double dev = 0;
  
  for (size_t i = 0; i < individuals; i++) {
    MatrixXi d_i = d[i];
    MatrixXd w_i = w[i];

    for (size_t j = 0; j < clusters; j++) {
      double b_j = beta_ext[j];

      for (size_t k = 0; k < loci; k++) {
        double b_k = beta_ext[clusters + k];
        
        double eta = b_j + b_k;
        double mu = linkinv(eta);
        int d_ele = d_i(j, k);
        
        double dev_contrib = 0.0;

        if (d_ele == 0) {  
          // y == 0: dev = 2*log((1+p)/(1-p))
          double p = (mu < 1e-6) ? 0.5 * mu : (sqrt(1.0 + mu * mu) - 1.0) / mu;
          dev_contrib = 2 * log((1.0 + p) / (1.0 - p));
        } else {
          // y != 0
          double d_ele_double = (double)d_ele;
          dev_contrib = 2 * (loglikeh(d_ele_double, d_ele_double) - loglikeh(mu, d_ele_double));
        }
        
        dev_contrib *= w_i(j, k);
        dev += dev_contrib;
      }
    }
  }  

  old_deviance = deviance;
  deviance = dev;
  deviance_calculated = true;
  
  double c = std::abs(dev - old_deviance)/(0.1 + std::abs(dev));
  
  if (verbose) {
    OUTPUT << "      Deviance convergence investigation:" << std::endl;
    OUTPUT << "        deviance     = " << dev << std::endl;
  }

  if (c < eps_deviance) {
    if (verbose) {
      OUTPUT << "        criteria     = " << c << " < " << eps_deviance << " [CONVERGENCE]" << std::endl;
    }
    
    return true;
  }

  if (verbose) {
    OUTPUT << "        criteria     = " << c << " >= " << eps_deviance << std::endl;
  }
  
  return false;
}

bool beta_converged(const VectorXd& beta, const VectorXd& old_beta, double eps_beta, bool verbose) {
  VectorXd den = 0.1 + beta.cwiseAbs().array();
  VectorXd diff = beta - old_beta;
  
  double c = diff.cwiseAbs().cwiseQuotient(den).maxCoeff();

  if (verbose) {
    OUTPUT << "      Beta convergence investigation:" << std::endl;
    OUTPUT << "        beta         = " << beta.transpose() << std::endl;
  }
    
  if (c < eps_beta) {
    if (verbose) {
      OUTPUT << "        criteria     = " << c << " < " << eps_beta << " [CONVERGENCE]" << std::endl;
    }
    
    return true;
  }

  if (verbose) {
    OUTPUT << "        criteria     = " << c << " >= " << eps_beta << std::endl;
  }
  
  return false;
}

void fill_generics(size_t clusters, size_t loci, VectorXd& beta_ext, 
  MatrixXd& lin_pred_generic, MatrixXd& mu_generic) {
  
  for (size_t j = 0; j < clusters; j++) {
    double beta_j = beta_ext[j];
    
    for (size_t k = 0; k < loci; k++) {
      double eta = beta_j + beta_ext[clusters + k];
      lin_pred_generic(j, k) = eta;
      mu_generic(j, k) = linkinv(eta);
    }
  }
}

void fill_H_ab(size_t individuals, size_t clusters, size_t loci,   
  const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w, 
  const MatrixXd& mu_generic,
  MatrixXd& H, MatrixXd& ab) {
   
  for (size_t i = 0; i < individuals; i++) {
    MatrixXi d_i = d[i];
    MatrixXd w_i = w[i];
        
    for (size_t j = 0; j < clusters; j++) {
      for (size_t k = 0; k < loci; k++) {
        double mu_generic_jk = mu_generic(j, k);
        double w_ijk = w_i(j, k);
        
        // H
        double psi_jk = w_ijk * varfunc(mu_generic_jk);
        H(j, k) += psi_jk;

        // ab
        double d_mu_res = d_i(j, k) - mu_generic_jk;
        ab(j, k) += w_ijk * d_mu_res;
      }
    }
  }
}

/*
void fill_H_ab(size_t individuals, size_t clusters, size_t loci,   
  const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w, 
  const MatrixXd& mu_generic,
  MatrixXd& H, MatrixXd& ab) {
  
  // varfunc(double mu) = mu * sqrt(1.0 + mu*mu)
  MatrixXd mu_1_sq = mu_generic.array().square() + 1;
  MatrixXd mu_sqrt_1_sq = mu_1_sq.cwiseSqrt();
  MatrixXd mu_varfunc = mu_generic.cwiseProduct(mu_sqrt_1_sq);
   
  for (size_t i = 0; i < individuals; i++) {
    MatrixXi d_i = d[i];
    MatrixXd w_i = w[i];
    MatrixXd d_mu = d_i.cast<double>().array() - mu_generic.array();
    
    H = H.array() + w_i.cwiseProduct(mu_varfunc).array();
    ab = ab.array() + w_i.cwiseProduct(d_mu).array();
  }
}
*/

void update_beta(const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w, VectorXd& beta, 
  MatrixXd& vcov, 
  size_t individuals, size_t clusters, size_t loci, bool verbose) {
  
  VectorXd beta_ext = get_beta_extended(beta);
  MatrixXd lin_pred_generic(clusters, loci);
  MatrixXd mu_generic(clusters, loci);
  fill_generics(clusters, loci, beta_ext, lin_pred_generic, mu_generic);
  
  MatrixXd H = MatrixXd::Zero(clusters, loci);
  MatrixXd ab = MatrixXd::Zero(clusters, loci);
  fill_H_ab(individuals, clusters, loci, d, w, mu_generic, H, ab);
  
  VectorXd Drvec_raw = H.colwise().sum();
  VectorXd Drvec(loci - 1);
  for (size_t k = 0; k < (loci - 1); k++) {
    Drvec[k] = Drvec_raw[k];
  }  
  MatrixXd Dr = Drvec.asDiagonal();
  
  VectorXd Dcinvvec(clusters);
  for (size_t j = 0; j < clusters; j++) {
    Dcinvvec[j] = 1 / H.row(j).sum();
  }
  MatrixXd Dcinv = Dcinvvec.asDiagonal();

  H.conservativeResize(NoChange, loci - 1);
  
  MatrixXd E = Dr - H.transpose() * Dcinv * H;
  MatrixXd Einv = E.lu().solve(MatrixXd::Identity(loci - 1, loci - 1));
  MatrixXd F = Dcinv * H;
  MatrixXd Ft = F.transpose(); 
  
  /*
  OUTPUT << "H:" << std::endl << H << std::endl;
  OUTPUT << "Dr:" << std::endl << Dr << std::endl;
  OUTPUT << "Dcinv:" << std::endl << Dcinv << std::endl;
  OUTPUT << "E:" << std::endl << E << std::endl;
  OUTPUT << "Einv:" << std::endl << Einv << std::endl;
  OUTPUT << "F:" << std::endl << F << std::endl;
  OUTPUT << "Ft:" << std::endl << Ft << std::endl;
  */
  
  MatrixXd P = MatrixXd::Zero(clusters + loci - 1, clusters + loci - 1);
  P.topLeftCorner(clusters, clusters) = Dcinv + (F * Einv) * Ft;
  P.topRightCorner(clusters, loci - 1) = -F * Einv;
  P.bottomLeftCorner(loci - 1, clusters) = -Einv * Ft;
  P.bottomRightCorner(loci - 1, loci - 1) = Einv;  
  
  vcov = P;

  VectorXd a = ab.rowwise().sum();  
  VectorXd b = ab.colwise().sum();  
  VectorXd gamma(clusters + loci - 1);
  
  for (size_t j = 0; j < clusters; j++) {
    gamma[j] = a[j];
  }
  
  for (size_t k = 0; k < (loci - 1); k++) {
    gamma[clusters + k] = b[k];
  }
  
  VectorXd beta_correction = P * gamma;
  
  beta = beta + beta_correction;
}

/*
// iter >= ensures old_deviance has been calculated
if (deviance_calculated && iter >= 2 && !isinf(old_deviance) && old_deviance < dev) {
  OUTPUT << "old_deviance = " << old_deviance << std::endl;
  OUTPUT << "dev          = " << dev << std::endl;
  throw "Deviance increased! Normally step-halving would be tried, but it has yet to be implemented...";
}
*/

/*
Returns true if converged, false otherwise
*/
bool irls(const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w, 
  VectorXd& beta, VectorXd& lin_pred, 
  MatrixXd& vcov,
  double& deviance, size_t& iterations,
  bool verbose = true, size_t maxit = 25, 
  bool use_deviance = true, double eps_deviance = 1e-6, 
  bool use_beta = true, double eps_beta = 1e-6,
  bool force_calculate_deviance = false) {

  check_input(d, w, beta, maxit, use_deviance, eps_deviance, use_beta, eps_beta);

  assert(maxit > 0);  
  assert(!use_deviance || (use_deviance && eps_deviance > 0));
  assert(!use_beta || (use_beta && eps_beta > 0));
  
  size_t individuals = d.size();
  assert(individuals > 0);
  assert(w.size() == individuals);

  size_t clusters = d.at(0).rows();
  size_t loci = d.at(0).cols();
  
  bool converged = false;
  
  bool deviance_calculated = false;
  double dev = std::numeric_limits<double>::infinity();
  double old_deviance = std::numeric_limits<double>::infinity();
  
  VectorXd start_beta = VectorXd(beta);
  
  VectorXd old_beta(beta.size());
  for (size_t idx = 0; idx < old_beta.size(); idx++) {
    old_beta[idx] = std::numeric_limits<double>::infinity();
  }

  if (verbose) {
    OUTPUT << "    Initial coefficients = " << beta.transpose() << std::endl;
  }
  
  iterations = 0;
  
  for (size_t iter = 0; iter < maxit; iter++) {
    iterations = iterations + 1;

    if (verbose) {
      OUTPUT << "    IRLS iteration " << (iter + 1) << std::endl;
    }
    
    if (use_deviance && use_beta) {
      bool beta_conv = beta_converged(beta, old_beta, eps_beta, verbose);
      
      if (beta_conv) {
        bool dev_conv = deviance_converged(d, w, deviance_calculated, lin_pred, beta, dev, old_deviance, eps_deviance, individuals, clusters, loci, verbose);
        
        if (dev_conv) {
          converged = true;
          break;
        }
      }
    } else if (use_deviance) {
      bool dev_conv = deviance_converged(d, w, deviance_calculated, lin_pred, beta, dev, old_deviance, eps_deviance, individuals, clusters, loci, verbose);
      
      if (dev_conv) {
        converged = true;
        break;
      }
    } else if (use_beta) {
      bool beta_conv = beta_converged(beta, old_beta, eps_beta, verbose);
      
      if (beta_conv) {
        converged = true;
        break;      
      }
    } else {
      throw "Unexpected error (!use_deviance && !use_beta)";
    }
    
    old_beta = beta;
    update_beta(d, w, beta, vcov, individuals, clusters, loci, verbose);
    
    #if defined(DISCLAPMIX_USED_IN_R)
    R_CheckUserInterrupt();
    #endif
  }
  
  if (!deviance_calculated) {
    if (force_calculate_deviance) {
      deviance_converged(d, w, deviance_calculated, lin_pred, beta, dev, old_deviance, eps_deviance, individuals, clusters, loci, false);
    }
  }

  deviance = dev;
  
  return converged;
}

