#include "init.h"
#include "input_validation.h"

using namespace Eigen;

VectorXd init_beta(const MatrixXi& x, const MatrixXi& y) {
  check_x_y(x, y);
  
  size_t clusters = y.rows();
  size_t loci = x.cols();
  
  VectorXd beta(clusters + loci - 1);

  for (size_t j = 0; j < clusters; j++) { 
    beta[j] = -1;  
  }
  
  for (size_t k = 0; k < (loci - 1); k++) { 
    beta[clusters + k] = 0.5;  
  }
  
  /*
  for (size_t j = 0; j < clusters; j++) { 
    beta[j] = 0.1;
  }
  
  for (size_t k = 0; k < (loci - 1); k++) { 
    beta[clusters + k] = 0.1;  
  }
  */
  
  return beta;
}

VectorXd init_tau(const MatrixXi& y) {
  check_y(y);
  
  size_t clusters = y.rows();
  
  VectorXd tau(clusters);
  for (size_t j = 0; j < clusters; j++) { 
    tau[j] = 1.0 / (double)clusters;  
  }  
  
  return tau;
}

std::vector<MatrixXi> init_response_vec(const MatrixXi& x, const MatrixXi& y) {
  check_x_y(x, y);

  size_t individuals = x.rows();
  size_t clusters = y.rows();
  size_t loci = x.cols();
  
  std::vector<MatrixXi> d(individuals);
  //std::vector<MatrixXi> d;
  //d.reserve(individuals); // reserve does not initialise values

  for (size_t i = 0; i < individuals; i++) {
    d[i] = MatrixXi(clusters, loci);
    
    for (size_t j = 0; j < clusters; j++) {
      for (size_t k = 0; k < loci; k++) {
        d[i](j, k) = std::abs(x(i, k) - y(j, k));
      }
    }    
  }
  
  return d;
}

std::vector<MatrixXd> init_weight_vec(size_t individuals, size_t clusters, size_t loci) {
  std::vector<MatrixXd> w(individuals);
  //std::vector<MatrixXd> w;
  //w.reserve(individuals); // reserve does not initialise values

  for (size_t i = 0; i < individuals; i++) {
    w[i] = MatrixXd(clusters, loci);
    
    for (size_t j = 0; j < clusters; j++) {
      for (size_t k = 0; k < loci; k++) {
        w[i](j, k) = 1.0 / clusters;        
      }
    }    
  }
  
  return w;
}

MatrixXd init_v_mat(const MatrixXi& x, const MatrixXi& y) {
  check_x_y(x, y);
  
  size_t individuals = x.rows();
  size_t clusters = y.rows();
  
  MatrixXd v(individuals, clusters);

  for (size_t i = 0; i < individuals; i++) {
    for (size_t j = 0; j < clusters; j++) {
      v(i, j) = 1.0 / (double)clusters;
    }
  }
  
  return v;
}

