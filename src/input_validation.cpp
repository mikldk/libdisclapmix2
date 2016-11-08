#include <vector>
#include <iostream>

#include "input_validation.h"

using namespace Eigen;

void check_common_input(size_t maxit, bool use_deviance, double eps_deviance, bool use_beta, double eps_beta) {
  if (maxit <= 0) {
    throw "maxit must be > 0";
  }
  
  if (use_deviance == false && use_beta == false) {
    throw "Either use_deviance or use_beta (or both) must be true";
  }
  
  if (use_deviance && eps_deviance <= 0) {
    throw "eps_deviance must be > 0";
  }
  
  if (use_beta && eps_beta <= 0) {
    throw "eps_beta must be > 0";
  }  
}

void check_y(const MatrixXi& y) {
/*  if (y == NULL) {
    throw "y cannot be NULL";
  }
  */
  size_t clusters = y.rows();
  size_t loci = y.cols();
  
  if (clusters <= 0) {
    throw "At least one cluster is needed";
  }
  
  if (loci <= 0) {
    throw "At least one locus is needed";
  }
}

void check_x_y(const MatrixXi& x, const MatrixXi& y) {
/*  if (x == NULL) {
    throw "x cannot be NULL";
  }
  
  if (y == NULL) {
    throw "y cannot be NULL";
  }
  */
  size_t individuals = x.rows();
  size_t clusters = y.rows();
  size_t loci = x.cols();
  
  if (individuals <= 0) {
    throw "At least one individual is needed";
  }
  
  if (clusters <= 0) {
    throw "At least one cluster is needed";
  }
  
  if (loci <= 0) {
    throw "At least one locus is needed";
  }
  
  if (loci != y.cols()) {    
    throw "x and y must have the same number of columns";
  } 
}


void check_d_w(const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w) {
/* if (d == NULL) {
    throw "d cannot be NULL";
  }
  
  if (w == NULL) {
    throw "w cannot be NULL";
  }
  */
  size_t individuals = d.size();
  
  if (individuals != w.size()) {    
    throw "d and w must not differ in size";
  } else if (individuals <= 0) {
    throw "At least one individual needed";
  }
  
  size_t clusters = d.at(0).rows();
  size_t loci = d.at(0).cols();
  
  if (clusters <= 0) {
    throw "clusters must be at least 1";
  } else if (loci <= 0) {
    throw "loci must be at least 1";
  }
}

void check_beta(const VectorXd& beta, size_t clusters, size_t loci) {
/*  if (beta == NULL) {
    throw "beta cannot be NULL (initialise one by using init_beta function)";
  }
  */
  if (beta.size() != (clusters + loci - 1)) {
    throw "beta must have size = clusters + loci - 1";
  }
}

void check_beta(const VectorXd& beta, const MatrixXi& x, const MatrixXi& y) {
  check_x_y(x, y);
  size_t clusters = y.rows();
  size_t loci = x.cols();
   
  check_beta(beta, clusters, loci);
}

void check_beta(const VectorXd& beta, const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w) {
  check_d_w(d, w);
  size_t clusters = d.at(0).rows();
  size_t loci = d.at(0).cols();
   
  check_beta(beta, clusters, loci);
}

void check_tau(const VectorXd& tau, const MatrixXi& y) {
  check_y(y);
  size_t clusters = y.rows();
   
/*  if (tau == NULL) {
    throw "tau cannot be NULL (initialise one by using init_tau function)";
  }
  */
  if (tau.size() != clusters) {
    throw "tau must have size = clusters";
  }
}

void check_v_mat(const MatrixXd& v_mat, size_t individuals, size_t clusters) {
/*  if (v_mat == NULL) {
    throw "v_mat cannot be NULL (initialise one by using init_v_mat function)";
  }
  */
  if (v_mat.rows() != individuals) {
    throw "v_mat must have rows() = individuals";
  }
  
  if (v_mat.cols() != clusters) {
    throw "v_mat must have cols() = clusters";
  }
}

void check_v_mat(const MatrixXd& v_mat, const MatrixXi& x, const MatrixXi& y) {
  check_x_y(x, y);
  size_t individuals = x.rows();
  size_t clusters = y.rows();
   
  check_v_mat(v_mat, individuals, clusters);
}

