#include <math.h>

#include "likelihood.h"

double loglikelihood_full(
  const MatrixXi& x, 
  const MatrixXi& y,
  const MatrixXd& v_mat, 
  const MatrixXd& disclap_parameters, 
  const VectorXd& tau) {
  
  check_x_y(x, y);

  size_t individuals = x.rows();
  size_t clusters = y.rows();
  size_t loci = x.cols();
  
  double logL = 0.0d;
  
  double loci_inv = 1.0d / (double)loci;
  
  for (size_t i = 0; i < individuals; i++) {
    for (size_t j = 0; j < clusters; j++) {
      double tau_j_norm = pow(tau[j], loci_inv);
      double v_ij = v_mat(i, j);
      
      for (size_t k = 0; k < loci; k++) {
        double d = std::abs(x(i, k) - y(j, k));
        double p_jk = disclap_parameters(j, k);
        double density = pow(p_jk, d) * (1.0d - p_jk) / (1.0d + p_jk);
        
        double logLterm = v_ij * log( tau_j_norm * density );
        
        // FIXME: Why is term NaN?
        // 0 * log(0) = 0 as x*log(x) -> 0 as x -> 0:
        if (v_ij == 0 && density == 0) {
          continue; // skip term (instead of adding 0)
        }
        /*
        if (isnan(logLterm)) {
          OUTPUT << "logLterm       = " << logLterm << std::endl;
          OUTPUT << "    v_ij       = " << v_ij << std::endl;
          OUTPUT << "    tau_j_norm = " << tau_j_norm << std::endl;
          OUTPUT << "    density    = " << density << std::endl;
        }
        */
        
        logL += logLterm;
      }
    }
  }
  
  return logL;
}

double loglikelihood_marginal(
  const MatrixXi& x, 
  const MatrixXi& y,
  const MatrixXd& disclap_parameters, 
  const VectorXd& tau) {
  
  VectorXd haplotype_probabilities = predict(x, y, disclap_parameters, tau);

  double logL = haplotype_probabilities.array().log().array().sum();
  
  return logL;
}

/*
double information_criterion(double logL, size_t parameters, double m) {
  double p = (double)parameters;

  return m*p - 2*logL;
}

double AIC(double logL, size_t parameters) {
  return information_criterion(logL, parameters, 2.0d);
}

double AICc(double AIC, size_t parameters, size_t observations) {
  double k = (double)parameters;
  double n = (double)observations;
  
  return AIC + 2*k*(k+1) / (n-k-1);
}

double BIC(double logL, size_t parameters, size_t observations) {
  return information_criterion(logL, parameters, (double)observations);
}
*/

