#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include "common.h"
#include "input_validation.h"
#include "predict.h"

double loglikelihood_full(
  const MatrixXi& x, 
  const MatrixXi& y,
  const MatrixXd& v_mat, 
  const MatrixXd& disclap_parameters, 
  const VectorXd& tau);

double loglikelihood_marginal(
  const MatrixXi& x, 
  const MatrixXi& y,
  const MatrixXd& disclap_parameters, 
  const VectorXd& tau);

/*
double AIC(double logL, size_t parameters);
double AICc(double AIC, size_t parameters, size_t observations);
double BIC(double logL, size_t parameters, size_t observations);
*/

#endif

