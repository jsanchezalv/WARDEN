#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector disc_instant_rcpp(double lcldr, NumericVector lclcurtime, NumericVector lclval) {
  int n = lclval.size();
  NumericVector result(n);
  
  if (lclcurtime.size() == 1) lclcurtime = rep(lclcurtime, n);
  
  for (int i = 0; i < n; ++i) {
    result[i] = lclval[i] * std::pow(1.0 + lcldr, -lclcurtime[i]);
  }
  
  return result;
}
