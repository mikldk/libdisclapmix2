#include <math.h>

#include "predict.h"

VectorXd predict(const MatrixXi& new_data, const MatrixXi& y, const MatrixXd& disclap_parameters, const VectorXd& tau) {
  check_x_y(new_data, y);

  size_t individuals = new_data.rows();
  size_t clusters = y.rows();
  size_t loci = new_data.cols();
  
  VectorXd happrobs(individuals);
  
  for (int i = 0; i < individuals; i++) {
    VectorXi h = new_data.row(i);
    double hprob = 0.0;
    
    for (int j = 0; j < clusters; j++) {
      VectorXi yhap = y.row(j);
      double component_prob = tau[j];
      
      for (int k = 0; k < loci; k++) {
        double p_jk = disclap_parameters(j, k);
        
        component_prob *= pow(p_jk, abs(h(k) - yhap(k)))*((1 - p_jk)/(1 + p_jk));
      }
      
      hprob += component_prob;
    }
    
    happrobs(i) = hprob;
  }
  
  return happrobs;
}

