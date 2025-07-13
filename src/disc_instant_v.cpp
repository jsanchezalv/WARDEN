#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector disc_instant_rcpp(double lcldr, NumericVector lclcurtime, NumericVector lclval) {
  int n = lclval.size();
  NumericVector result(n);
  double dr = 1.0 + lcldr;
  
  for (int i = 0; i < n; ++i) {
    result[i] = lclval[i] * std::pow(dr, -lclcurtime[i]);
  }
  
  return result;
}
