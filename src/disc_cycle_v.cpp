#include <Rcpp.h>
using namespace Rcpp;

//' Cycle discounting for vectors
//'
//' @param lcldr The discount rate
//' @param lclprvtime The time of the previous event in the simulation
//' @param cyclelength The cycle length
//' @param lclcurtime The time of the current event in the simulation
//' @param lclval The  value to be discounted
//' @param starttime The start time for accrual of cycle costs (if not 0)
//' @param max_cycles The maximum number of cycles
//'
//' @return Double vector based on cycle discounting
//' 
//' @details
//' This function per cycle discounting, i.e., considers that the cost/qaly is accrued
//' per cycles, and performs it automatically without needing to create new events.
//' It can accommodate changes in cycle length/value/starttime (e.g., in the case of 
//' induction and maintenance doses) within the same item.
//'
//' @examples 
//' disc_cycle_v(lcldr=0.03, lclprvtime=0, cyclelength=1/12, lclcurtime=2, lclval=500,starttime=0)
//' disc_cycle_v(
//'  lcldr=0.000001,
//'  lclprvtime=0,
//'  cyclelength=1/12,
//'  lclcurtime=2,
//'  lclval=500,
//'  starttime=0,
//'  max_cycles = 4)
//' 
//' #Here we have a change in cycle length, max number of cylces and starttime at time 2
//'  #(e.g., induction to maintenance)
//' #In the model, one would do this by redifining cycle_l, max_cycles and starttime
//'  #of the corresponding item at a given event time. 
//' disc_cycle_v(lcldr=0,
//'  lclprvtime=c(0,1,2,2.5),
//'  cyclelength=c(1/12, 1/12,1/2,1/2),
//'  lclcurtime=c(1,2,2.5,4), lclval=c(500,500,500,500),
//'  starttime=c(0,0,2,2), max_cycles = c(24,24,2,2)
//'   )
//' @export
// [[Rcpp::export]]
NumericVector disc_cycle_v(
    double lcldr,
    NumericVector lclprvtime,
    NumericVector cyclelength,
    NumericVector lclcurtime,
    NumericVector lclval,
    NumericVector starttime,
    Nullable<NumericVector> max_cycles = R_NilValue
) {
  int n = lclcurtime.size();
  NumericVector result(n);
  NumericVector max_c(n, NA_REAL);
  
  if (max_cycles.isNotNull()) {
    max_c = as<NumericVector>(max_cycles);
    if (max_c.size() == 1) {
      std::fill(max_c.begin(), max_c.end(), max_c[0]);
    }
  }
  
  if (lclprvtime.size() == 1 && n > 1) lclprvtime = rep(lclprvtime, n);
  if (cyclelength.size() == 1 && n > 1) cyclelength = rep(cyclelength, n);
  if (lclval.size() == 1 && n > 1) lclval = rep(lclval, n);
  if (starttime.size() == 1 && n > 1) starttime = rep(starttime, n);
  
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
    
    double s = std::pow(1.0 + lcldr, clen) - 1;
    double d = (prv == 0 ? 0 : prv / clen) - 1;
    
    if (lcldr == 0) {
      result[i] = val * max_eff;
    } else {
      result[i] = val * (1 - std::pow(1 + s, -max_eff)) / (s * std::pow(1 + s, d));
    }
  }
  
  return result;
}
