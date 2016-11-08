#include <Eigen/SVD>

#include "beta.h"
#include "input_validation.h"

using namespace Eigen;

/*
beta = [clus_1, clus_2, ..., clus_c, loc_1, loc_2, ..., loc_r]
*/

VectorXd get_beta_extended(const VectorXd& beta) {
  size_t m = beta.size();
  VectorXd beta_ext(m + 1);
  beta_ext[m] = 0;
  for (size_t idx = 0; idx < m; idx++) {
    beta_ext[idx] = beta[idx];
  }
  
  return beta_ext;
}

void beta_description(const MatrixXi& x, const MatrixXi& y, const VectorXd& beta) {
  check_x_y(x, y);
  
  size_t clusters = y.rows();
  size_t loci = x.cols();
  
  for (size_t j = 0; j < clusters; j++) {
    OUTPUT << "clus_" << (j+1) << " ";
  }
  
  for (size_t k = 0; k < (loci - 1); k++) {
    if (k > 0) {
      OUTPUT << " ";
    }
    
    OUTPUT << "loc_" << (k+1);    
  }
}

MatrixXd beta_to_discrete_laplace_parameters(const MatrixXi& x, const MatrixXi& y, const VectorXd& beta) {
  check_x_y(x, y);
  
  size_t clusters = y.rows();
  size_t loci = x.cols();
  
  VectorXd beta_ext = get_beta_extended(beta);
  MatrixXd pars(clusters, loci);
  
  for (size_t j = 0; j < clusters; j++) {
    double b_j = beta_ext[j];
    
    for (size_t k = 0; k < loci; k++) {
      double b_k = beta_ext[clusters + k];
      double p = std::exp(b_j + b_k);
      pars(j, k) = p;      
    }
  }
  
  return pars;
}

