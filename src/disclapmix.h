#ifndef DISCLAPMIX_H
#define DISCLAPMIX_H

// cat(paste0(paste("#include \"", sort(list.files("/home/mikl/work-aau/software/libdisclapmix/code/libdisclapmix", "hpp")), "\"", sep = ""), collapse = "\n"))
#include "beta.h"
#include "common.h"
#include "disclap_family.h"
//#include "disclapmix.h"
#include "init.h"
#include "input_validation.h"
#include "irls.h"
#include "likelihood.h"
#include "predict.h"

class IterationsInfo {
  public:
    MatrixXi y;
    std::vector<VectorXd> beta_vec;
    std::vector<double> deviance_vec;
    std::vector<double> logL_full_vec;
    std::vector<double> logL_marginal_vec;
};

bool disclapmix(const MatrixXi& x, MatrixXi& y, 
  VectorXd& beta, VectorXd& tau, MatrixXd& v_mat, MatrixXd& vcov,
  double& deviance, 
  std::vector<IterationsInfo>& iter_info,
  size_t& model_observations, 
  size_t& model_parameters, 
  size_t& iter_irls,
  size_t& iter_em,
  size_t& iter_centers,
  size_t verbose_level, size_t maxit_centers_change, size_t maxit_disclap_em, size_t maxit_irls, 
  double eps_v_mat,
  bool use_deviance, double eps_deviance, 
  bool use_beta, double eps_beta);
  
#endif

