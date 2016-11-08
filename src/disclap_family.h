#ifndef DISCLAP_FAMILY_H
#define DISCLAP_FAMILY_H

#include "common.h"

double loglikeh(double mu, double y);
double linkinv(double eta);
double mu_eta(double eta);
double varfunc(double mu);

#endif

