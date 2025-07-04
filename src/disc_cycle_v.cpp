#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector disc_cycle_rcpp(
    double lcldr,
    NumericVector lclprvtime,
    NumericVector cyclelength,
    NumericVector lclcurtime,
    NumericVector lclval,
    NumericVector starttime,
    Nullable<NumericVector> max_cycles = R_NilValue
) {
  int n = lclval.size();
  NumericVector result(n);
  NumericVector max_c(n, NA_REAL);
  
  if (max_cycles.isNotNull()) {
    max_c = as<NumericVector>(max_cycles);
    if (max_c.size() == 1) {
      std::fill(max_c.begin(), max_c.end(), max_c[0]);
    }
  }
  
  if (lclprvtime.size() == 1) lclprvtime = rep(lclprvtime, n);
  if (cyclelength.size() == 1) cyclelength = rep(cyclelength, n);
  if (lclcurtime.size() == 1) lclcurtime = rep(lclcurtime, n);
  if (starttime.size() == 1) starttime = rep(starttime, n);
  
  for (int i = 0; i < n; ++i) {
    double prv = lclprvtime[i];
    double cur = lclcurtime[i];
    double start = starttime[i];
    double clen = cyclelength[i];
    double val = lclval[i];
    
    double total_cycles = std::ceil((cur - start) / clen);
    double prev_cycles = std::ceil((prv - start) / clen);
    if (total_cycles < 0) total_cycles = 0;
    if (prev_cycles < 0) prev_cycles = 0;
    
    double n_cycles = total_cycles - prev_cycles;
    
    // Add 1 cycle if current time is exactly a cycle boundary
    if (prv == cur && fmod(cur - start, clen) == 0) {
      n_cycles += 1;
    }
    
    double max_eff = n_cycles;
    if (!NumericVector::is_na(max_c[i])) {
      double rem = std::max(0.0, max_c[i] - prev_cycles);
      max_eff = std::min(n_cycles, rem);
    }
    
    double s = std::pow(1.0 + lcldr, clen) - 1.0;
    double d = (prv == 0 ? 0 : prv / clen) - 1;
    
    if (lcldr == 0) {
      result[i] = val * max_eff;
    } else {
      result[i] = val * (1 - std::pow(1 + s, -max_eff)) / (s * std::pow(1 + s, d));
    }
  }
  
  return result;
}
