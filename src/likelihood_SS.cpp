#include <Rcpp.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include "lmvnd_ll.h"
//#include "lgd.h"
using namespace std;

// [[Rcpp::export]]
RcppExport SEXP SSLL(SEXP Cobs_, SEXP C_, SEXP Sigma_, SEXP Cpred_, SEXP sd_, SEXP init_mean_, SEXP init_sd_, int N) {
  //vector<double> mean = Rcpp::as<vector<double> >(mean_);
  //vector<double> Sigma = Rcpp::as<vector<double> >(Sigma_);
  //vector<double> X = Rcpp::as<vector<double> >(X_);
  double log_lik = 0;
  Rcpp::NumericMatrix Cobs = Rcpp::as<Rcpp::NumericMatrix >(Cobs_);
  Rcpp::NumericMatrix C = Rcpp::as<Rcpp::NumericMatrix >(C_);
  Rcpp::NumericVector Sigma = Rcpp::as<Rcpp::NumericVector >(Sigma_);
  Rcpp::NumericVector sd = Rcpp::as<Rcpp::NumericVector >(sd_);
  Rcpp::NumericVector init_mean = Rcpp::as<Rcpp::NumericVector>(init_mean_);
  Rcpp::NumericVector init_sd = Rcpp::as<Rcpp::NumericVector>(init_sd_);
  Rcpp::NumericMatrix Cpred = Rcpp::as<Rcpp::NumericMatrix >(Cpred_);
  Rcpp::NumericVector Cobsrow = Cobs(0, Rcpp::_ );
  Rcpp::NumericVector Crow = C(0,Rcpp::_);
  Rcpp::NumericVector Cprow = Cpred(0,Rcpp::_);
  log_lik = LMVND_LL::lmvnd(Rcpp::as<vector<double> >(Sigma), Rcpp::as<vector<double> > (Cobsrow), Rcpp::as<vector<double> > (Crow));
  log_lik = log_lik + LMVND_LL::lmvnd(Rcpp::as<vector<double> >(init_sd), Rcpp::as<vector<double> >(Crow), Rcpp::as<vector<double> >(init_mean));
  for (int i=1; i<N; i++) {
    Cobsrow = Cobs(i,Rcpp::_);
    Crow = C(i,Rcpp::_);
    Cprow = Cpred(i, Rcpp::_);
    log_lik = log_lik + LMVND_LL::lmvnd(Rcpp::as<vector<double> >(Sigma), Rcpp::as<vector<double> > (Cobsrow), Rcpp::as<vector<double> > (Crow));
    log_lik = log_lik + LMVND_LL::lmvnd(Rcpp::as<vector<double> >(sd), Rcpp::as<vector<double> > (Crow), Rcpp::as<vector<double> > (Cprow));
  }
  //double a = Rcpp::as<double>(a_);
  //double b = Rcpp::as<double>(b_);
  //vector<double> X = Rcpp::as<vector<double> >(sd_);
  
  //for (int j=0; j<5; j++) {
  //  double lgx = X[j];
  //  log_lik = log_lik + LMVND_LL::lgammad(lgx, a, b);
  //}
  return Rcpp::wrap(log_lik);
}