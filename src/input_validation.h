#ifndef INPUT_VALIDATION_H
#define INPUT_VALIDATION_H

#include "common.h"

void check_common_input(size_t maxit, bool use_deviance, double eps_deviance, bool use_beta, double eps_beta);
void check_y(const MatrixXi& y);
void check_x_y(const MatrixXi& x, const MatrixXi& y);
void check_d_w(const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w);

void check_beta(const VectorXd& beta, const MatrixXi& x, const MatrixXi& y);
void check_beta(const VectorXd& beta, const std::vector<MatrixXi>& d, const std::vector<MatrixXd>& w);

void check_tau(const VectorXd& tau, const MatrixXi& y);

void check_v_mat(const MatrixXd& v_mat, const MatrixXi& x, const MatrixXi& y);

#endif

