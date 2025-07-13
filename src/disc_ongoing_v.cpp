#include <Rcpp.h>
using namespace Rcpp;

//' Calculate discounted costs and qalys between events for vectors
//'
//' @param lcldr The discount rate 
//' @param lclprvtime The time of the previous event in the simulation
//' @param lclcurtime The time of the current event in the simulation
//' @param lclval The value to be discounted
//'
//' @return Double based on continuous time discounting
//'
//' @examples 
//' disc_ongoing_v(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)
//' 
//'
//' @export
// [[Rcpp::export]]
NumericVector disc_ongoing_v(double lcldr,
                                NumericVector lclprvtime,
                                NumericVector lclcurtime,
                                NumericVector lclval) {
  
  int n = lclcurtime.size();
  NumericVector result(n);
  
  double instantdr = std::log1p(lcldr);
  
  if (lcldr == 0) {
    for (int i = 0; i < n; ++i) {
      double t0 = lclprvtime.size() == 1 ? lclprvtime[0] : lclprvtime[i];
      double t1 = lclcurtime.size() == 1 ? lclcurtime[0] : lclcurtime[i];
      double val = lclval.size() == 1 ? lclval[0] : lclval[i];
      result[i] = val * (t1 - t0);
    }
  } else {
    for (int i = 0; i < n; ++i) {
      double t0 = lclprvtime.size() == 1 ? lclprvtime[0] : lclprvtime[i];
      double t1 = lclcurtime.size() == 1 ? lclcurtime[0] : lclcurtime[i];
      double val = lclval.size() == 1 ? lclval[0] : lclval[i];
      result[i] = (val / (-instantdr)) * (std::exp(-instantdr * t1) - std::exp(-instantdr * t0));
    }
  }
  
  return result;
}
