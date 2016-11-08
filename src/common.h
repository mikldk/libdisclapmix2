#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

#if defined(DISCLAPMIX_USED_IN_R)
  #include <R.h>
  #include <Rcpp.h>
  #define OUTPUT Rcpp::Rcout
  #define FLUSH_OUTPUT R_FlushConsole()
#else
  #define OUTPUT std::cout
  #define FLUSH_OUTPUT std::cout << std::flush
#endif

//const IOFormat disclapEigenFmt(FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
const int disclapEigenFmtPrecision = 6;
const IOFormat disclapEigenFmt(disclapEigenFmtPrecision, 0, ", ", ";\n", "(", ")", "[", "]");

#endif

