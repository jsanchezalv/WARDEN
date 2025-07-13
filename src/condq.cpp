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

// [[Rcpp::export]]
NumericVector qcond_exp_cpp(NumericVector rnd, NumericVector rate) {
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

// [[Rcpp::export]]
NumericVector qcond_weibull_cpp(NumericVector rnd, NumericVector shape, NumericVector scale, NumericVector lower_bound = NumericVector::create(0.0)) {
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

// [[Rcpp::export]]
NumericVector qcond_weibullPH_cpp(NumericVector rnd, NumericVector shape, NumericVector scale, NumericVector lower_bound = NumericVector::create(0.0)) {
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

// [[Rcpp::export]]
NumericVector qcond_llogis_cpp(NumericVector rnd, NumericVector shape, NumericVector scale, NumericVector lower_bound = NumericVector::create(0.0)) {
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

// [[Rcpp::export]]
NumericVector qcond_gompertz_cpp(NumericVector rnd, NumericVector shape, NumericVector rate, NumericVector lower_bound = NumericVector::create(0.0)) {
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

// [[Rcpp::export]]
NumericVector qcond_lnorm_cpp(NumericVector rnd, NumericVector meanlog, NumericVector sdlog, NumericVector lower_bound, NumericVector s_obs) {
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

// [[Rcpp::export]]
NumericVector qcond_norm_cpp(NumericVector rnd, NumericVector mean, NumericVector sd, NumericVector lower_bound, NumericVector s_obs) {
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

// [[Rcpp::export]]
NumericVector qcond_gamma_cpp(NumericVector rnd, NumericVector shape, NumericVector rate, NumericVector lower_bound, NumericVector s_obs) {
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

// Gompertz survival (flexsurv), also survival with lower.tail=FALSE
double sf_gompertz(double time, double shape, double rate) {
  Environment flexsurv_env = Environment::namespace_env("flexsurv");
  
  // Access the internal pgompertz_work symbol
  Function pgompertz_work = flexsurv_env["pgompertz_work"];
  
  // Arguments: q, shape, rate, lower.tail, log.p
  NumericVector res = pgompertz_work(
    NumericVector::create(time),
    NumericVector::create(shape),
    NumericVector::create(rate),
    LogicalVector::create(false),  // lower.tail = FALSE
    LogicalVector::create(false)   // log.p = FALSE
  );
  
  return res[0];
}

NumericVector luck_adj_rcpp(NumericVector prevsurv,
                            NumericVector cursurv,
                            NumericVector luck,
                            bool condq);

// [[Rcpp::export]]
List qtimecov_cpp(double luck,
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
      return qcond_exp_cpp(luck, a);
    };
  } else if (dist == "gamma") {
    sf_fun = sf_gamma;
    qcond_fun = qcond_gamma_cpp;
  } else if (dist == "lnorm") {
    sf_fun = sf_lnorm;
    qcond_fun = qcond_lnorm_cpp;
  } else if (dist == "norm") {
    sf_fun = sf_norm;
    qcond_fun = qcond_norm_cpp;
  } else if (dist == "weibull") {
    sf_fun = sf_weibull;
    qcond_fun = [](NumericVector luck, NumericVector a, NumericVector b, NumericVector t, NumericVector surv_prev) {
      return qcond_weibull_cpp(luck, a, b, t);
    };
  } else if (dist == "weibullPH") {
    sf_fun = sf_weibullPH;
    qcond_fun = [](NumericVector luck, NumericVector a, NumericVector b, NumericVector t, NumericVector surv_prev) {
      return qcond_weibullPH_cpp(luck, a, b, t);
    };
  } else if (dist == "llogis") {
    sf_fun = sf_llogis;
    qcond_fun = [](NumericVector luck, NumericVector a, NumericVector b, NumericVector t, NumericVector surv_prev) {
      return qcond_llogis_cpp(luck, a, b, t);
    };
  } else if (dist == "gompertz") {
    sf_fun = sf_gompertz;
    qcond_fun = [](NumericVector luck, NumericVector a, NumericVector b, NumericVector t, NumericVector surv_prev) {
      return qcond_gompertz_cpp(luck, a, b, t);
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
    luck = luck_adj_rcpp(NumericVector::create(surv_prev_dt),
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