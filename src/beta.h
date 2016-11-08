#ifndef BETA_H
#define BETA_H

#include "common.h"

VectorXd get_beta_extended(const VectorXd& beta);
void beta_description(const MatrixXi& x, const MatrixXi& y, const VectorXd& beta);
MatrixXd beta_to_discrete_laplace_parameters(const MatrixXi& x, const MatrixXi& y, const VectorXd& beta);

#endif

