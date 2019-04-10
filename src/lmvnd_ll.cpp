#include <Rcpp.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include "lmvnd_ll.h"
// load header files and libraries
using namespace std;

// [[Rcpp::export]]
double lmvnd(vector<double> Sigma, vector<double> X, vector<double> mean, int imax) {
  //vector<double> mean = Rcpp::as<vector<double> >(mean_);
  //vector<double> Sigma = Rcpp::as<vector<double> >(Sigma_);
  //vector<double> X = Rcpp::as<vector<double> >(X_);
  double log_det = 0;
  double diff = 0;
  double log_exp = 0;
  double ldens;
  for (int i=0; i<imax; i++) {
    log_det = log_det + -.5*log(Sigma[i]);
    diff = X[i] - mean[i];
    log_exp = log_exp + -.5*pow(diff,2) / Sigma[i];
  }
  ldens = log_det + log_exp;
  return ldens;
}