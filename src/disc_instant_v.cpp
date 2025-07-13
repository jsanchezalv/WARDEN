#include <Rcpp.h>
using namespace Rcpp;

//' Calculate instantaneous discounted costs or qalys for vectors
//'
//' @param lcldr The discount rate
//' @param lclcurtime The time of the current event in the simulation
//' @param lclval The value to be discounted
//'
//' @return Double based on discrete time discounting
//'
//' @examples
//' disc_instant_v(lcldr=0.035, lclcurtime=3, lclval=2500)
//' 
//'
//' @export
// [[Rcpp::export]]
NumericVector disc_instant_v(double lcldr, NumericVector lclcurtime, NumericVector lclval) {
  int n = lclcurtime.size();
  NumericVector result(n);
  double dr = 1.0 + lcldr;
  
  for (int i = 0; i < n; ++i) {
    double val = lclval.size() == 1 ? lclval[0] : lclval[i];
    result[i] = val * std::pow(dr, -lclcurtime[i]);
  }
  
  return result;
}
