#ifndef IRLS_H
#define IRLS_H

#include "common.h"

bool irls(const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w, 
  VectorXd& beta, VectorXd& lin_pred, MatrixXd& vcov,
  double& deviance, size_t& iterations, 
  bool verbose, size_t maxit,
  bool use_deviance, double eps_deviance, bool use_beta, double eps_beta,
  bool force_calculate_deviance);

#endif

