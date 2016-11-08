#ifndef PREDICT_H
#define PREDICT_H

#include "common.h"
#include "input_validation.h"

VectorXd predict(const MatrixXi& new_data, const MatrixXi& y, const MatrixXd& disclap_parameters, const VectorXd& tau);

#endif

