#include <Rcpp.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <math.h>
#include <stdio.h>
using namespace std;

// [[Rcpp::export]]
RcppExport SEXP lgammad(SEXP x_, SEXP a_, SEXP b_) {
  double x = Rcpp::as<double>(x_);
  double b = Rcpp::as<double>(b_);
  double a = Rcpp::as<double>(a_);
  double ldens;
  ldens = a*log(b) - lgamma(a) + (a-1)*log(x) - x*b;
  return Rcpp::wrap(ldens);
}