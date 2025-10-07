#include <Rcpp.h>
using namespace Rcpp;

//' Perform luck adjustment
//'
//' @param prevsurv Value of the previous survival
//' @param cursurv Value of the current survival
//' @param luck Luck used to be adjusted (number between 0 and 1)
//' @param condq Conditional quantile approach or standard approach
//'
//' @export
//' 
//' @return Adjusted luck number between 0 and 1
//'
//' @details
//' This function performs the luck adjustment automatically for the user, returning the adjusted luck number.
//' Luck is interpreted in the same fashion as is standard in R (higher luck, higher time to event).
//' 
//' Note that if TTE is predicted using a conditional quantile function (e.g., conditional gompertz, conditional quantile weibull...) `prevsurv` and `cursurv`
//' are the unconditional survival using the "previous" parametrization but at the previous time for `presurv` and at the current time for `cursurv`.
//' For other distributions, `presurv` is the survival up to current time using the previous parametrization, and `cursurv` 
//' is the survival up to current time using the current parametrization.
//' 
//' Note that the advantage of the conditional quantile function is that it does not need the new parametrization to update the luck, 
//' which makes this approach computationally more efficient.
//' This function can also work with vectors, which could allow to update multiple lucks in a single approach, and it can preserve names
//'
//' @examples
//' luck_adj(prevsurv = 0.8,
//'  cursurv = 0.7,
//'  luck = 0.5,
//'  condq = TRUE)
//'  
//' luck_adj(prevsurv = c(1,0.8,0.7),
//'  cursurv = c(0.7,0.6,0.5),
//'  luck = setNames(c(0.5,0.6,0.7),c("A","B","C")),
//'  condq = TRUE)
//'  
//' luck_adj(prevsurv = 0.8,
//'  cursurv = 0.7,
//'  luck = 0.5,
//'  condq = FALSE) #different results
//' 
//' #Unconditional approach, timepoint of change is 25,
//' # parameter goes from 0.02 at time 10 to 0.025 to 0.015 at time 25,
//' #  starting luck is 0.37
//' new_luck <- luck_adj(prevsurv = 1 - pweibull(q=10,3,1/0.02),
//'  cursurv = 1 - pweibull(q=10,3,1/0.025),
//'  luck = 0.37,
//'  condq = FALSE) #time 10 change
//'  
//' new_luck <- luck_adj(prevsurv = 1 - pweibull(q=25,3,1/0.025),
//'  cursurv = 1 - pweibull(q=25,3,1/0.015),
//'  luck = new_luck,
//'  condq = FALSE) #time 25 change
//'  
//' qweibull(new_luck, 3, 1/0.015) #final TTE 
//' 
//' #Conditional quantile approach 
//' new_luck <- luck_adj(prevsurv = 1-pweibull(q=0,3,1/0.02),
//'                       cursurv = 1- pweibull(q=10,3,1/0.02),
//'                       luck = 0.37,
//'                       condq = TRUE) #time 10 change, previous time is 0 so prevsurv will be 1
//' 
//' new_luck <- luck_adj(prevsurv = 1-pweibull(q=10,3,1/0.025),
//'                       cursurv = 1- pweibull(q=25,3,1/0.025),
//'                       luck = new_luck,
//'                       condq = TRUE) #time 25 change
//' 
//' qcond_weibull(rnd = new_luck,
//'                      shape = 3,
//'                      scale = 1/0.015,
//'                      lower_bound = 25) + 25 #final TTE
// [[Rcpp::export]]
NumericVector luck_adj(NumericVector prevsurv,
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
