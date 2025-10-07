// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// Helper to recycle scalar values
inline NumericVector recycle_or_check(NumericVector x, int n) {
  if (x.size() == 1) return rep(x[0], n);
  if (x.size() != n) stop("Input vectors must have length 1 or equal to 'rnd'");
  return x;
}

//' Conditional quantile function for exponential distribution 
//'
//' @param rnd Vector of quantiles
//' @param rate The rate parameter
//' 
//' Note taht the conditional quantile for an exponential is independent of time due to constant hazard
//'
//' @return Estimate(s) from the conditional exponential distribution based on given parameters
//'
//' @export
//'
//' @examples
//' qcond_exp(rnd = 0.5,rate = 3)
// [[Rcpp::export]]
NumericVector qcond_exp(NumericVector rnd, NumericVector rate) {
  int n = rnd.size();
  NumericVector rate_ = recycle_or_check(rate, n);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    if (rnd[i] < 0.0 || rnd[i] > 1.0) {
      out[i] = NA_REAL;
    } else {
      out[i] = -std::log1p(-rnd[i]) / rate_[i];
    }
  }
  return out;
}

//' Conditional quantile function for weibull distribution 
//'
//' @param rnd Vector of quantiles
//' @param shape The shape parameter as in R stats package weibull
//' @param scale The scale parameter as in R stats package weibull
//' @param lower_bound The lower bound to be used (current time)
//'
//' @return Estimate(s) from the conditional weibull distribution based on given parameters
//'
//' @export
//'
//' @examples
//' qcond_weibull(rnd = 0.5,shape = 3,scale = 66.66,lower_bound = 50)
// [[Rcpp::export]]
NumericVector qcond_weibull(NumericVector rnd, NumericVector shape, NumericVector scale, NumericVector lower_bound = NumericVector::create(0.0)) {
  int n = rnd.size();
  NumericVector shape_ = recycle_or_check(shape, n);
  NumericVector scale_ = recycle_or_check(scale, n);
  NumericVector lb_ = recycle_or_check(lower_bound, n);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    if (rnd[i] < 0.0 || rnd[i] > 1.0) {
      out[i] = NA_REAL;
    } else {
      out[i] = std::pow(std::pow(lb_[i] / scale_[i], shape_[i]) - std::log1p(-rnd[i]), 1.0 / shape_[i]) * scale_[i] - lb_[i];
    }
  }
  return out;
}

//' Conditional quantile function for WeibullPH (flexsurv)
//'
//' @param rnd Vector of quantiles (between 0 and 1)
//' @param shape Shape parameter of WeibullPH
//' @param scale Scale (rate) parameter of WeibullPH (i.e., as in hazard = scale * t^(shape - 1))
//' @param lower_bound Lower bound (current time)
//'
//' @return Estimate(s) from the conditional weibullPH distribution based on given parameters
//'
//' @export
//'
//' @examples
//' qcond_weibullPH(rnd = 0.5, shape = 2, scale = 0.01, lower_bound = 5)
// [[Rcpp::export]]
NumericVector qcond_weibullPH(NumericVector rnd, NumericVector shape, NumericVector scale, NumericVector lower_bound = NumericVector::create(0.0)) {
  int n = rnd.size();
  NumericVector shape_ = recycle_or_check(shape, n);
  NumericVector scale_ = recycle_or_check(scale, n);
  NumericVector lb_ = recycle_or_check(lower_bound, n);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    if (rnd[i] < 0.0 || rnd[i] > 1.0) {
      out[i] = NA_REAL;
    } else {
      double base = std::pow(lb_[i], shape_[i]) - std::log1p(-rnd[i]) / scale_[i];
      out[i] = base < 0 ? NA_REAL : std::pow(base, 1.0 / shape_[i]) - lb_[i];
    }
  }
  return out;
}

//' Conditional quantile function for loglogistic distribution 
//'
//' @param rnd Vector of quantiles
//' @param shape The shape parameter
//' @param scale The scale parameter
//' @param lower_bound The lower bound to be used (current time)
//'
//' @return Estimate(s) from the conditional loglogistic distribution based on given parameters
//'
//' @export
//'
//' @examples
//' qcond_llogis(rnd = 0.5,shape = 1,scale = 1,lower_bound = 1)
// [[Rcpp::export]]
NumericVector qcond_llogis(NumericVector rnd, NumericVector shape, NumericVector scale, NumericVector lower_bound = NumericVector::create(0.0)) {
  int n = rnd.size();
  NumericVector shape_ = recycle_or_check(shape, n);
  NumericVector scale_ = recycle_or_check(scale, n);
  NumericVector lb_ = recycle_or_check(lower_bound, n);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    if (rnd[i] < 0.0 || rnd[i] > 1.0) {
      out[i] = NA_REAL;
    } else {
      out[i] = scale_[i] * std::pow((1 / (1 - rnd[i])) * (rnd[i] + std::pow(lb_[i] / scale_[i], shape_[i])), 1.0 / shape_[i]) - lb_[i];
    }
  }
  return out;
}

//' Quantile function for conditional Gompertz distribution (lower bound only)
//'
//' @param rnd Vector of quantiles
//' @param shape The shape parameter of the Gompertz distribution, defined as in the coef() output on a flexsurvreg object
//' @param rate The rate parameter of the Gompertz distribution, defined as in the coef() output on a flexsurvreg object
//' @param lower_bound The lower bound of the conditional distribution
//'
//' @return Estimate(s) from the conditional Gompertz distribution based on given parameters
//'
//' @export
//'
//' @examples
//' qcond_gompertz(rnd=0.5,shape=0.05,rate=0.01,lower_bound = 50)
// [[Rcpp::export]]
NumericVector qcond_gompertz(NumericVector rnd, NumericVector shape, NumericVector rate, NumericVector lower_bound = NumericVector::create(0.0)) {
  int n = rnd.size();
  NumericVector shape_ = recycle_or_check(shape, n);
  NumericVector rate_ = recycle_or_check(rate, n);
  NumericVector lb_ = recycle_or_check(lower_bound, n);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    if (rnd[i] < 0.0 || rnd[i] > 1.0) {
      out[i] = NA_REAL;
    } else {
      out[i] = (1.0 / shape_[i]) * std::log1p(-((shape_[i] / rate_[i]) * std::log1p(-rnd[i]) / (std::expm1(shape_[i] * lb_[i]) + 1)));
    }
  }
  return out;
}

//' Conditional quantile function for lognormal distribution 
//'
//' @param rnd Vector of quantiles
//' @param meanlog The meanlog parameter
//' @param sdlog The sdlog parameter
//' @param lower_bound The lower bound to be used (current time)
//' @param s_obs is the survival observed up to lower_bound time,
//'  normally defined from time 0 as 1 - plnorm(q = lower_bound, meanlog, sdlog) but may be different if parametrization has changed previously
//'
//' @importFrom stats qnorm
//'
//' @return Estimate(s) from the conditional lognormal distribution based on given parameters
//'
//' @export
//'
//' @examples
//' qcond_lnorm(rnd = 0.5, meanlog = 1,sdlog = 1,lower_bound = 1, s_obs=0.8)
// [[Rcpp::export]]
NumericVector qcond_lnorm(NumericVector rnd, NumericVector meanlog, NumericVector sdlog, NumericVector lower_bound, NumericVector s_obs) {
  int n = rnd.size();
  NumericVector meanlog_ = recycle_or_check(meanlog, n);
  NumericVector sdlog_ = recycle_or_check(sdlog, n);
  NumericVector lb_ = recycle_or_check(lower_bound, n);
  NumericVector s_obs_ = recycle_or_check(s_obs, n);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    if (rnd[i] < 0.0 || rnd[i] > 1.0) {
      out[i] = NA_REAL;
    } else {
      out[i] = std::expm1(meanlog_[i] + sdlog_[i] * R::qnorm(1 - s_obs_[i] * (1 - rnd[i]), 0.0, 1.0, 1, 0)) + 1 - lb_[i];
    }
  }
  return out;
}

//' Conditional quantile function for normal distribution 
//'
//' @param rnd Vector of quantiles
//' @param mean The mean parameter
//' @param sd The sd parameter
//' @param lower_bound The lower bound to be used (current time)
//' @param s_obs is the survival observed up to lower_bound time,
//'  normally defined from time 0 as 1 - pnorm(q = lower_bound, mean, sd) but may be different if parametrization has changed previously
//'
//' @importFrom stats qnorm
//'
//' @return Estimate(s) from the conditional normal distribution based on given parameters
//'
//' @export
//'
//' @examples
//' qcond_norm(rnd = 0.5, mean = 1,sd = 1,lower_bound = 1, s_obs=0.8)
// [[Rcpp::export]]
NumericVector qcond_norm(NumericVector rnd, NumericVector mean, NumericVector sd, NumericVector lower_bound, NumericVector s_obs) {
  int n = rnd.size();
  NumericVector mean_ = recycle_or_check(mean, n);
  NumericVector sd_ = recycle_or_check(sd, n);
  NumericVector lb_ = recycle_or_check(lower_bound, n);
  NumericVector s_obs_ = recycle_or_check(s_obs, n);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    if (rnd[i] < 0.0 || rnd[i] > 1.0) {
      out[i] = NA_REAL;
    } else {
      out[i] = R::qnorm(1 - s_obs_[i] * (1 - rnd[i]), mean_[i], sd_[i], 1, 0) - lb_[i];
    }
  }
  return out;
}

//' Conditional quantile function for gamma distribution 
//'
//' @param rnd Vector of quantiles
//' @param rate The rate parameter
//' @param shape The shape parameter
//' @param lower_bound The lower bound to be used (current time)
//' @param s_obs is the survival observed up to lower_bound time,
//'  normally defined from time 0 as 1 - pgamma(q = lower_bound, rate, shape) but may be different if parametrization has changed previously
//'
//' @importFrom stats qgamma
//'
//' @return Estimate(s) from the conditional gamma distribution based on given parameters
//'
//' @export
//'
//' @examples
//' qcond_gamma(rnd = 0.5, shape = 1.06178, rate = 0.01108,lower_bound = 1, s_obs=0.8)
// [[Rcpp::export]]
NumericVector qcond_gamma(NumericVector rnd, NumericVector shape, NumericVector rate, NumericVector lower_bound, NumericVector s_obs) {
  int n = rnd.size();
  NumericVector shape_ = recycle_or_check(shape, n);
  NumericVector rate_ = recycle_or_check(rate, n);
  NumericVector lb_ = recycle_or_check(lower_bound, n);
  NumericVector s_obs_ = recycle_or_check(s_obs, n);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    if (rnd[i] < 0.0 || rnd[i] > 1.0) {
      out[i] = NA_REAL;
    } else {
      out[i] = R::qgamma(1 - s_obs_[i] * (1 - rnd[i]), shape_[i], 1.0 / rate_[i], 1, 0) - lb_[i];
    }
  }
  return out;
}


// --- Survival functions wrappers --- //
// We'll call R's survival functions via Rcpp::Function


// Exponential survival function: S(t) = 1 - F(t)
// R::pexp parameters: p, rate, lower.tail, log.p
double sf_exp(double time, double rate) {
  return R::pexp(time, 1/rate, /*lower.tail=*/false, /*log.p=*/false);
}

// Gamma survival function: S(t) = 1 - F(t)
// R::pgamma parameters: p, shape, scale, lower.tail, log.p
// Note: R uses scale = 1/rate, so convert rate -> scale
double sf_gamma(double time, double shape, double rate) {
  double scale = 1.0 / rate;
  return R::pgamma(time, shape, scale, /*lower.tail=*/false, /*log.p=*/false);
}

// Lognormal survival function: S(t) = 1 - F(t)
// R::plnorm parameters: q, meanlog, sdlog, lower.tail, log.p
double sf_lnorm(double time, double meanlog, double sdlog) {
  return R::plnorm(time, meanlog, sdlog, /*lower.tail=*/false, /*log.p=*/false);
}

// Normal survival function: S(t) = 1 - F(t)
// R::pnorm parameters: q, mean, sd, lower.tail, log.p
double sf_norm(double time, double mean, double sd) {
  return R::pnorm(time, mean, sd, /*lower.tail=*/false, /*log.p=*/false);
}

// Weibull survival function: S(t) = 1 - F(t)
// R::pweibull parameters: q, shape, scale, lower.tail, log.p
double sf_weibull(double time, double shape, double scale) {
  return R::pweibull(time, shape, scale, /*lower.tail=*/false, /*log.p=*/false);
}

// For WeibullPH, Loglogistic and Gompertz: call flexsurv via R (no direct C API available)

// Weibull Proportional Hazards survival (flexsurv)
// S(t) = 1 - F(t), so subtract from 1
double sf_weibullPH(double time, double shape, double scale) {
  if (time < 0) return 1.0;
  return (std::expm1(- scale * std::pow(time, shape)) +1);
}

// Loglogistic survival (flexsurv), already returns survival with lower.tail=FALSE
double sf_llogis(double time, double shape, double scale) {
  if (time < 0) return 1.0;
  double ratio = time / scale;
  return 1.0 / (1.0 + std::pow(ratio, shape));
}

// Helper functions from flexsurv
inline double exprel(const double x) {
  if (x != 0) {
    return expm1(x) / x;
  } else {
    return 1;
  }
}

inline double safe_coeff(const double q, const double shape, const double rate) {
  if (!std::isinf(q)) {
    const double scale_q = shape * q;
    return -rate * q * exprel(scale_q);
  } else {
    if (shape < 0) {
      return rate / shape;
    } else {
      return R_NegInf;
    }
  }
}
// Gompertz survival (flexsurv), also survival with lower.tail=FALSE
double sf_gompertz(double time, double shape, double rate) {
  // Check for invalid parameters
  if (rate < 0) {
    Rcpp::warning("Negative rate parameter");
    return NA_REAL;
  }
  
  // If time < 0, survival = 1
  if (time < 0) {
    return 1.0;
  }
  
  if (shape != 0) {
    const double coeff = safe_coeff(time, shape, rate);
    return std::exp(coeff);  // This is the (!give_log) & (!lower_tail) case
  } else {
    // When shape = 0, becomes exponential: survival = exp(-rate * time)
    return std::exp(-rate * time);
  }
}

NumericVector luck_adj(NumericVector prevsurv,
                            NumericVector cursurv,
                            NumericVector luck,
                            bool condq);

//' Draw Time-to-Event with Time-Dependent Covariates and Luck Adjustment
//'
//' Simulate a time-to-event (TTE) from a parametric distribution with parameters varying over time.
//' User provides parameter functions and distribution name. The function uses internal survival and
//' conditional quantile functions, plus luck adjustment to simulate the event time. See
//' the vignette on avoiding cycles for an example in a model.
//' 
//'
//' @param luck Numeric between 0 and 1. Initial random quantile (luck).
//' @param a_fun Function of time .time returning the first distribution parameter (e.g., rate, shape, meanlog).
//' @param b_fun Function of time .time returning the second distribution parameter (e.g., scale, sdlog). Defaults to a function returning NA.
//' @param dist Character string specifying the distribution. Supported: "exp", "gamma", "lnorm", "norm", "weibull", "llogis", "gompertz".
//' @param dt Numeric. Time step increment to update parameters and survival. Default 0.1.
//' @param max_time Numeric. Max allowed event time to prevent infinite loops. Default 100.
//' @param start_time Numeric. Time to use as a starting point of reference (e.g., curtime).
//'
//' @return List with simulated time-to-event and final luck value.
//' 
//' 
//' @importFrom flexsurv pllogis pgompertz pweibullPH
//' @importFrom stats pexp pgamma plnorm pnorm pweibull
//' 
//' @export
//' 
//' @details The objective of this function is to avoid the user to have cycle events
//' with the only scope of updating some variables that depend on time and re-evaluate
//' a TTE. The idea is that this function should only be called at start and when 
//' an event impacts a variable (e.g., stroke event impacting death TTE), in which case
//' it would need to be called again at that point. In that case, the user would need to 
//' call e.g., `a <- qtimecov` with `max_time = curtime` arguments,
//' and then call it again with no max_time, and
//' `luck = a$luck, start_time=a$tte` (so there is no need to add curtime to the resulting time).
//' 
//' It's recommended to play with `dt` argument to balance running time and precision of the estimates.
//' For example, if we know we only update the equation annually (not continuously),
//' then we could just set `dt = 1`, which would make computations faster.
//'
//' @examples
//'
//' param_fun_factory <- function(p0, p1, p2, p3) {
//'   function(.time) p0 + p1*.time + p2*.time^2 + p3*(floor(.time) + 1)
//' }
//' 
//' set.seed(42)
//' 
//' # 1. Exponential Example
//' rate_exp <- param_fun_factory(0.1, 0, 0, 0)
//' qtimecov(
//'   luck = runif(1),
//'   a_fun = rate_exp,
//'   dist = "exp"
//' )
//' 
//' 
//' # 2. Gamma Example
//' shape_gamma <- param_fun_factory(2, 0, 0, 0)
//' rate_gamma <- param_fun_factory(0.2, 0, 0, 0)
//' qtimecov(
//'   luck = runif(1),
//'   a_fun = shape_gamma,
//'   b_fun = rate_gamma,
//'   dist = "gamma"
//' )
//' 
//' 
//' # 3. Lognormal Example
//' meanlog_lnorm <- param_fun_factory(log(10) - 0.5*0.5^2, 0, 0, 0)
//' sdlog_lnorm <- param_fun_factory(0.5, 0, 0, 0)
//' qtimecov(
//'   luck = runif(1),
//'   a_fun = meanlog_lnorm,
//'   b_fun = sdlog_lnorm,
//'   dist = "lnorm"
//' )
//' 
//' 
//' # 4. Normal Example
//' mean_norm <- param_fun_factory(10, 0, 0, 0)
//' sd_norm <- param_fun_factory(2, 0, 0, 0)
//' qtimecov(
//'   luck = runif(1),
//'   a_fun = mean_norm,
//'   b_fun = sd_norm,
//'   dist = "norm"
//' )
//' 
//' 
//' # 5. Weibull Example
//' shape_weibull <- param_fun_factory(2, 0, 0, 0)
//' scale_weibull <- param_fun_factory(10, 0, 0, 0)
//' qtimecov(
//'   luck = runif(1),
//'   a_fun = shape_weibull,
//'   b_fun = scale_weibull,
//'   dist = "weibull"
//' )
//' 
//' 
//' # 6. Loglogistic Example
//' shape_llogis <- param_fun_factory(2.5, 0, 0, 0)
//' scale_llogis <- param_fun_factory(7.6, 0, 0, 0)
//' qtimecov(
//'   luck = runif(1),
//'   a_fun = shape_llogis,
//'   b_fun = scale_llogis,
//'   dist = "llogis"
//' )
//' 
//' 
//' # 7. Gompertz Example
//' shape_gomp <- param_fun_factory(0.01, 0, 0, 0)
//' rate_gomp <- param_fun_factory(0.091, 0, 0, 0)
//' qtimecov(
//'   luck = runif(1),
//'   a_fun = shape_gomp,
//'   b_fun = rate_gomp,
//'   dist = "gompertz"
//' )
//'
//' #Time varying example, with change at time 8
//' rate_exp <- function(.time) 0.1 + 0.01*.time * 0.00001*.time^2
//' rate_exp2 <- function(.time) 0.2 + 0.02*.time
//' time_change <- 8
//' init_luck <- 0.95
//'
//' a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005,
//'                       max_time = time_change)
//' qtimecov(luck = a$luck,a_fun = rate_exp2,dist = "exp", dt = 0.005, start_time=a$tte)
//' 
//' 
//' #An example of how it would work in the model, this would also work with time varying covariates!
//' rate_exp <- function(.time) 0.1
//' rate_exp2 <- function(.time) 0.2
//' rate_exp3 <- function(.time) 0.3
//' time_change <- 10 #evt 1
//' time_change2 <- 15 #evt2
//' init_luck <- 0.95
//' #at start, we would just draw TTE
//' qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005)
//' 
//' #at event in which rate changes (at time 10) we need to do this:
//' a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005,
//'                       max_time = time_change)
//' new_luck <- a$luck
//' qtimecov(luck = new_luck,a_fun = rate_exp2,dist = "exp", dt = 0.005, start_time=a$tte)
//' 
//' #at second  event in which rate changes again (at time 15) we need to do this:
//' a <- qtimecov(luck = new_luck,a_fun = rate_exp2,dist = "exp", dt = 0.005,
//'                       max_time = time_change2, start_time=a$tte)
//' new_luck <- a$luck
//' #final TTE is
//' qtimecov(luck = new_luck,a_fun = rate_exp3,dist = "exp", dt = 0.005, start_time=a$tte)
// [[Rcpp::export]]
List qtimecov(double luck,
                  Function a_fun,
                  RObject b_fun = R_NilValue,
                  std::string dist = "exp",
                  double dt = 0.1,
                  double max_time = 100,
                  double start_time = 0) {
  
  if (a_fun.isNULL()) stop("a_fun must be a valid function");
  
  // Define function pointers for survival and conditional quantile, avoiding repeated if-else
  
  typedef double (*SfFun)(double, double, double);
  typedef NumericVector (*QCondFun)(NumericVector, NumericVector, NumericVector, NumericVector, NumericVector);
  
  SfFun sf_fun = nullptr;
  QCondFun qcond_fun = nullptr;
  
  if (dist == "exp") {
    sf_fun = [](double t, double rate, double unused) {
      return sf_exp(t, rate);
    };
    qcond_fun = [](NumericVector luck, NumericVector a, NumericVector b, NumericVector t, NumericVector surv_prev) {
      return qcond_exp(luck, a);
    };
  } else if (dist == "gamma") {
    sf_fun = sf_gamma;
    qcond_fun = qcond_gamma;
  } else if (dist == "lnorm") {
    sf_fun = sf_lnorm;
    qcond_fun = qcond_lnorm;
  } else if (dist == "norm") {
    sf_fun = sf_norm;
    qcond_fun = qcond_norm;
  } else if (dist == "weibull") {
    sf_fun = sf_weibull;
    qcond_fun = [](NumericVector luck, NumericVector a, NumericVector b, NumericVector t, NumericVector surv_prev) {
      return qcond_weibull(luck, a, b, t);
    };
  } else if (dist == "weibullPH") {
    sf_fun = sf_weibullPH;
    qcond_fun = [](NumericVector luck, NumericVector a, NumericVector b, NumericVector t, NumericVector surv_prev) {
      return qcond_weibullPH(luck, a, b, t);
    };
  } else if (dist == "llogis") {
    sf_fun = sf_llogis;
    qcond_fun = [](NumericVector luck, NumericVector a, NumericVector b, NumericVector t, NumericVector surv_prev) {
      return qcond_llogis(luck, a, b, t);
    };
  } else if (dist == "gompertz") {
    sf_fun = sf_gompertz;
    qcond_fun = [](NumericVector luck, NumericVector a, NumericVector b, NumericVector t, NumericVector surv_prev) {
      return qcond_gompertz(luck, a, b, t);
    };
  } else {
    stop("Unsupported distribution");
  }
  
  double time = start_time;
  
  // Determine if b_fun is a callable function
   
  // Evaluate a and b at start time
  double a_curr = as<double>(a_fun(time));
  bool type_b_fun_f;
  type_b_fun_f = TYPEOF(b_fun) == CLOSXP;
  double b_curr;
  if (b_fun.isNULL()) {
    b_curr = NA_REAL;
  } else if (type_b_fun_f) {
    Function b_fun_func(b_fun); // safe conversion now
    b_curr = as<double>(b_fun_func(time));
  } else {
    stop("b_fun is not a function or NULL");
  }
  // Compute initial survival at current time
  double surv_prev = sf_fun(time, a_curr, b_curr);
  
  // Compute residual time-to-event based on distribution
  double residual_tte = qcond_fun(NumericVector::create(luck),
                                  NumericVector::create(a_curr),
                                  NumericVector::create(b_curr),
                                  NumericVector::create(time),
                                  NumericVector::create(surv_prev))[0];
  
  
  if (residual_tte <= dt) {
    return List::create(Named("tte") = residual_tte + time, Named("luck") = luck);
  }
  
  while (true) {
    time += dt;
    
    double surv_curr = sf_fun(time, a_curr, b_curr);
    double surv_prev_dt = sf_fun(time - dt, a_curr, b_curr);
    
    a_curr = as<double>(a_fun(time));
    if (type_b_fun_f) {
      Function b_fun_func(b_fun); // safe conversion now
      b_curr = as<double>(b_fun_func(time));
    }

    // Update luck with survival adjustment
    luck = luck_adj(NumericVector::create(surv_prev_dt),
                         NumericVector::create(surv_curr),
                         NumericVector::create(luck),
                         true)[0];
    
    // Calculate residual time-to-event
    residual_tte = qcond_fun(NumericVector::create(luck),
                             NumericVector::create(a_curr),
                             NumericVector::create(b_curr),
                             NumericVector::create(time - dt),
                             NumericVector::create(surv_prev_dt))[0];
    // Rcout << "time: " << time << ", luck: " << luck << ", residual_tte: " << residual_tte
    //       << ", surv_prev_dt: " << surv_prev_dt << ", surv_curr: " << surv_curr << std::endl;
    
    double total_tte = time + residual_tte;
    
    if (residual_tte <= dt || total_tte <= time || time >= max_time) {
      return List::create(Named("tte") = std::min(total_tte, max_time), Named("luck") = luck);
    }
  }
  
  return NA_REAL;
}