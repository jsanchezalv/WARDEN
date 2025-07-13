#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector luck_adj_rcpp(NumericVector prevsurv,
                            NumericVector cursurv,
                            NumericVector luck,
                            bool condq = true) {
  
  int n = luck.size();
  NumericVector result(n);
  
  for (int i = 0; i < n; i++) {
    double p = prevsurv.size() == 1 ? prevsurv[0] : prevsurv[i];
    double c = cursurv.size() == 1 ? cursurv[0] : cursurv[i];
    double l = luck[i];
    
    double adj;
    if (p != 0) {
      if (condq) {
        adj = 1.0 - ((1.0 - l) * (p / c));
      } else {
        adj = 1.0 - ((1.0 - l) * (c / p));
      }
    } else {
      adj = l;
    }
    
    // Clamp to [1e-9, 0.999999999]
    adj = std::min(std::max(adj, 1e-9), 0.999999999);
    result[i] = adj;
  }
  
  // Preserve names
  result.attr("names") = luck.attr("names");
  
  return result;
}
