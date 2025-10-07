#include <Rcpp.h>
using namespace Rcpp;

//' Draw time to event (tte) from a Poisson or Poisson-Gamma (PG) Mixture/Negative Binomial (NB) Process using C++
//' 
//' @param n The number of observations to be drawn
//' @param rate rate of the event (events per unit time)
//' @param obs_time period over which events are observable
//' @param theta Optional. If provided, Poisson-Gamma (NB). Represents gamma shape.
//' @param t_reps Optional. Number of TBEs to be generated to capture events within the observation window.
//' @param seed Optional integer seed for reproducibility.
//' @param return_ind_rate Logical: include individual rate vector in output when theta provided.
//' @param return_df Logical: return a data.frame with event-level rows (if TRUE).
//' 
//' @return If return_df=TRUE: a data.frame (or NULL if no events). Else: list with tte and optionally ind_rate.
//' 
//' @examples
//' rpoisgamma_rcpp(1, rate = 1, obs_time = 1, theta = 1)
//' 
//' 
//' 
//' @export
// [[Rcpp::export]]
List rpoisgamma_rcpp(int n, double rate, Nullable<double> theta = R_NilValue,
                      double obs_time = 1.0, Nullable<int> t_reps = R_NilValue,
                      Nullable<int> seed = R_NilValue, bool return_ind_rate = false,
                      bool return_df = false) {
   if (n < 0) stop("n must be non-negative");
   if (obs_time <= 0.0) stop("obs_time must be positive");
   
   RNGScope rng;  // attach to R's RNG (does NOT set seed)
   
   int trep = 0;
   NumericVector rate_par;          // only used when theta is provided
   const bool has_theta = !theta.isNull();
   
   // Decide t_reps (deterministic, no RNG)
   if (t_reps.isNull()) {
     if (!has_theta) {
       trep = (int) R::qpois(0.9999, rate * obs_time, /*lower_tail=*/1, /*log_p=*/0);
     } else {
       const double theta_val = as<double>(theta);
       const double mu   = rate * obs_time;
       const double prob = theta_val / (theta_val + mu);
       trep = (int) R::qnbinom(0.9999, theta_val, prob, /*lower_tail=*/1, /*log_p=*/0);
     }
   } else {
     trep = as<int>(t_reps);
   }
   if (trep < 0) stop("t_reps must be non-negative");
   
   // NB/PG: draw individual rates *before* any seeding (this matches your R code)
   if (has_theta) {
     const double theta_val = as<double>(theta);
     const double scale = rate / theta_val;  // Gamma scale
     rate_par = NumericVector(n);
     for (int i = 0; i < n; ++i) {
       rate_par[i] = R::rgamma(theta_val, scale);
     }
   }
   
   // Seed ONCE immediately before drawing exponentials (matches your R code)
   if (!seed.isNull()) {
     int sd = as<int>(seed);
     Environment base_env("package:base");
     Function set_seed = base_env["set.seed"];
     set_seed(sd);
   }
   
   // Generate event times
   std::vector<std::vector<double>> ds(n);
   for (int i = 0; i < n; ++i) {
     const double individual_rate = has_theta ? rate_par[i] : rate;
     if (trep == 0) {
       ds[i] = std::vector<double>();
       continue;
     }
     
     std::vector<double> cum;
     cum.reserve(std::min(100, trep));
     double s = 0.0;
     
     // Draw exactly 'trep' exponentials, no early break (R draws rexp(t_reps, ...))
     for (int j = 0; j < trep; ++j) {
       double x = R::rexp(1.0 / individual_rate);   // mean = 1/rate
       s += x;
       if (s <= obs_time) {
         cum.push_back(s);
       }
     }
     ds[i] = std::move(cum);
   }
   
   if (return_df) {
     std::vector<int> ids, evt_num, evt_count;
     std::vector<double> ttes, t_btw_evt, ind_rate_col;
     
     for (int i = 0; i < n; ++i) {
       const auto &v = ds[i];
       const int m = (int)v.size();
       if (m == 0) continue;
       
       for (int j = 0; j < m; ++j) {
         ids.push_back(i + 1);
         ttes.push_back(v[j]);
         const double prev = (j == 0 ? 0.0 : v[j - 1]);
         t_btw_evt.push_back(v[j] - prev);
         evt_num.push_back(j + 1);
         evt_count.push_back(m);
         if (has_theta && return_ind_rate) {
           ind_rate_col.push_back(rate_par[i]);
         }
       }
     }
     
     if (ttes.empty()) {
       // return list(df = NULL) so $df is safe to access
       return List::create(Named("df") = R_NilValue);
     }
     
     DataFrame df = DataFrame::create(
       Named(".id")       = wrap(ids),
       Named("tte")       = wrap(ttes),
       Named("t_btw_evt") = wrap(t_btw_evt),
       Named("evt_num")   = wrap(evt_num),
       Named("evt_count") = wrap(evt_count)
     );
     if (has_theta && return_ind_rate) {
       df["ind_rate"] = wrap(ind_rate_col);
     }
     return List::create(Named("df") = df);
   }
   
   // list output (tte per id)
   List tte_list(n);
   for (int i = 0; i < n; ++i) tte_list[i] = wrap(ds[i]);
   
   if (return_ind_rate && has_theta) {
     return List::create(Named("tte") = tte_list, Named("ind_rate") = rate_par);
   } else if (return_ind_rate && !has_theta) {
     // keep interface identical: return NA vector for ind_rate
     return List::create(Named("tte") = tte_list, Named("ind_rate") = NumericVector(n, NA_REAL));
   } else {
     return List::create(Named("tte") = tte_list);
   }
}
