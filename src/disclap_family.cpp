#include <math.h>

#include "disclap_family.h"

double loglikeh(double mu, double y) {
  /*
    R:
    library(Rmpfr)
    mus_org <- 10^(-(1:14))
    mus <- mpfr(mus_org, precBits = 100)
    ps_cor <- (sqrt(1.0 + mus*mus) - 1.0) / mus
    ps_app <- 0.5*mus
    cbind(mus_org, abs(ps_cor - ps_app))
  */
  double p = (mu < 1e-6) ? 0.5*mu : (sqrt(1.0 + mu*mu) - 1.0) / mu;
  double logl = log(1.0-p) - log(1.0+p) + y*log(p);  
  return logl;
}

double linkinv(double eta) {
  // returns mu = g^-1(eta), where this is g^-1(eta) = 2*e^eta / (1 - e^(2*eta))
  double expeta = exp(eta);
  return 2*expeta / (1.0 - expeta*expeta);
}

double mu_eta(double eta) {
  double exp_eta = exp(eta);
  double exp_eta_sq = exp_eta * exp_eta;
  double exp_eta_sq_one = exp_eta_sq - 1.0;  
  return (2*exp_eta*(1.0 + exp_eta_sq)) / (exp_eta_sq_one * exp_eta_sq_one);
}

double varfunc(double mu) {
  return mu * sqrt(1.0 + mu*mu);
}

