#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector disc_ongoing_rcpp(double lcldr,
                                NumericVector lclprvtime,
                                NumericVector lclcurtime,
                                NumericVector lclval) {
  
  int n = lclval.size();
  NumericVector result(n);
  
  double instantdr = std::log(1.0 + lcldr);
  
  for (int i = 0; i < n; ++i) {
    double t0 = lclprvtime.size() == 1 ? lclprvtime[0] : lclprvtime[i];
    double t1 = lclcurtime.size() == 1 ? lclcurtime[0] : lclcurtime[i];
    double val = lclval[i];
    
    if (lcldr == 0) {
      result[i] = val * (t1 - t0);
    } else {
      result[i] = (val / (-instantdr)) * (std::exp(-instantdr * t1) - std::exp(-instantdr * t0));
    }
  }
  
  return result;
}
