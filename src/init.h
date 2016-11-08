#ifndef INIT_H
#define INIT_H

#include "common.h"

VectorXd init_beta(const MatrixXi& x, const MatrixXi& y);
VectorXd init_tau(const MatrixXi& y);
std::vector<MatrixXi> init_response_vec(const MatrixXi& x, const MatrixXi& y);
std::vector<MatrixXd> init_weight_vec(size_t individuals, size_t clusters, size_t loci);
MatrixXd init_v_mat(const MatrixXi& x, const MatrixXi& y);

#endif

