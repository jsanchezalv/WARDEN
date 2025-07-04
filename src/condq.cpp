// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
double qcond_exp_cpp(double rnd, double rate) {
  if (rnd < 0.0 || rnd > 1.0) stop("rnd is <0 or >1");
  return -std::log(1 - rnd) / rate;
}

// [[Rcpp::export]]
double qcond_weibull_cpp(double rnd, double shape, double scale, double lower_bound = 0.0) {
  if (rnd < 0.0 || rnd > 1.0) stop("rnd is <0 or >1");
  return std::pow(std::pow(lower_bound / scale, shape) - std::log(1 - rnd), 1.0 / shape) * scale - lower_bound;
}

// [[Rcpp::export]]
double qcond_weibullPH_cpp(double rnd, double shape, double scale, double lower_bound = 0.0) {
  if (rnd < 0.0 || rnd > 1.0) stop("rnd is <0 or >1");
  double base = std::pow(lower_bound, shape) - std::log(1 - rnd) / scale;
  if (base < 0) return NA_REAL;
  return std::pow(base, 1.0 / shape) - lower_bound;
}

// [[Rcpp::export]]
double qcond_llogis_cpp(double rnd, double shape, double scale, double lower_bound = 0.0) {
  if (rnd < 0.0 || rnd > 1.0) stop("rnd is <0 or >1");
  return scale * std::pow((1 / (1 - rnd)) * (rnd + std::pow(lower_bound / scale, shape)), 1.0 / shape) - lower_bound;
}

// [[Rcpp::export]]
double qcond_gompertz_cpp(double rnd, double shape, double rate, double lower_bound = 0.0) {
  if (rnd < 0.0 || rnd > 1.0) stop("rnd is <0 or >1");
  return (1.0 / shape) * std::log(1 - ((shape / rate) * std::log(1 - rnd) / std::exp(shape * lower_bound)));
}

// [[Rcpp::export]]
double qcond_lnorm_cpp(double rnd, double meanlog, double sdlog, double lower_bound, double s_obs) {
  if (rnd < 0.0 || rnd > 1.0) stop("rnd is <0 or >1");
  return std::exp(meanlog + sdlog * R::qnorm(1 - s_obs * (1 - rnd), 0.0, 1.0, 1, 0)) - lower_bound;
}

// [[Rcpp::export]]
double qcond_norm_cpp(double rnd, double mean, double sd, double lower_bound, double s_obs) {
  if (rnd < 0.0 || rnd > 1.0) stop("rnd is <0 or >1");
  return R::qnorm(1 - s_obs * (1 - rnd), mean, sd, 1, 0) - lower_bound;
}

// [[Rcpp::export]]
double qcond_gamma_cpp(double rnd, double shape, double rate, double lower_bound, double s_obs) {
  if (rnd < 0.0 || rnd > 1.0) stop("rnd is <0 or >1");
  return R::qgamma(1 - s_obs * (1 - rnd), shape, 1.0 / rate, 1, 0) - lower_bound;
}


// --- Survival functions wrappers --- //
// We'll call R's survival functions via Rcpp::Function

double sf_exp(double time, double rate) {
  Function pexp("pexp");
  return 1.0 - as<double>(pexp(time, Named("rate") = rate));
}

double sf_gamma(double time, double shape, double rate) {
  Function pgamma("pgamma");
  return 1.0 - as<double>(pgamma(time, Named("shape") = shape, Named("rate") = rate));
}

double sf_lnorm(double time, double meanlog, double sdlog) {
  Function plnorm("plnorm");
  return 1.0 - as<double>(plnorm(time, meanlog, sdlog));
}

double sf_norm(double time, double mean, double sd) {
  Function pnorm("pnorm");
  return 1.0 - as<double>(pnorm(time, mean, sd));
}

double sf_weibull(double time, double shape, double scale) {
  Function pweibull("pweibull");
  return 1.0 - as<double>(pweibull(time, Named("shape") = shape, Named("scale") = scale));
}

// For WeibullPH, Loglogistic and Gompertz, call flexsurv from R via Rcpp

double sf_weibullPH(double time, double shape, double scale) {
  Function pweibullPH("flexsurv::pweibullPH");
  return 1.0 - as<double>(pweibullPH(time, shape, scale));
}

double sf_llogis(double time, double shape, double scale) {
  Function pllogis("flexsurv::pllogis");
  return as<double>(pllogis(time, Named("shape") = shape, Named("scale") = scale, Named("lower.tail") = false));
}

double sf_gompertz(double time, double shape, double rate) {
  Function pgompertz("flexsurv::pgompertz");
  return as<double>(pgompertz(time, Named("shape") = shape, Named("rate") = rate, Named("lower.tail") = false));
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
                  bool return_luck = false,
                  double start_time = 0) {
  
  if (a_fun.isNULL()) stop("a_fun must be a valid function");
  
  // Define function pointers for survival and conditional quantile, avoiding repeated if-else
  
  typedef double (*SfFun)(double, double, double);
  typedef double (*QCondFun)(double, double, double, double, double);
  typedef double (*QCondFunSimple)(double, double); // for simpler qcond functions (like exp)
  
  SfFun sf_fun = nullptr;
  QCondFun qcond_fun = nullptr;
  QCondFunSimple qcond_fun_simple = nullptr;
  
  if (dist == "exp") {
    sf_fun = [](double t, double a, double b){ return sf_exp(t, a); };
    qcond_fun_simple = qcond_exp_cpp;
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
    qcond_fun = [](double luck, double a, double b, double t, double surv_prev) {
      return qcond_weibull_cpp(luck, a, b, t);
    };
  } else if (dist == "weibullPH") {
    sf_fun = sf_weibullPH;
    qcond_fun = [](double luck, double a, double b, double t, double surv_prev) {
      return qcond_weibullPH_cpp(luck, a, b, t);
    };
  } else if (dist == "llogis") {
    sf_fun = sf_llogis;
    qcond_fun = [](double luck, double a, double b, double t, double surv_prev) {
      return qcond_llogis_cpp(luck, a, b, t);
    };
  } else if (dist == "gompertz") {
    sf_fun = sf_gompertz;
    qcond_fun = [](double luck, double a, double b, double t, double surv_prev) {
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
  double residual_tte = 0;
  if (dist == "exp") {
    residual_tte = qcond_fun_simple(luck, a_curr);
  } else {
    residual_tte = qcond_fun(luck, a_curr, b_curr, time, surv_prev);
  }
  
  if (residual_tte <= dt) {
    if (return_luck) return List::create(Named("tte") = residual_tte + time, Named("luck") = luck);
    else return List::create(Named("tte") = residual_tte);
  }
  
  while (true) {
    time += dt;
    
    a_curr = as<double>(a_fun(time));
    if (type_b_fun_f) {
      Function b_fun_func(b_fun); // safe conversion now
      b_curr = as<double>(b_fun_func(time));
    }
    
    double surv_curr = sf_fun(time, a_curr, b_curr);
    double surv_prev_dt = sf_fun(time - dt, a_curr, b_curr);
    
    // Update luck with survival adjustment
    luck = luck_adj_rcpp(NumericVector::create(surv_prev_dt),
                         NumericVector::create(surv_curr),
                         NumericVector::create(luck),
                         true)[0];
    
    // Calculate residual time-to-event
    if (dist == "exp") {
      residual_tte = qcond_fun_simple(luck, a_curr);
    } else {
      residual_tte = qcond_fun(luck, a_curr, b_curr, time - dt, surv_prev_dt);
    }
    
    double total_tte = time + residual_tte;
    
    if (residual_tte <= dt || total_tte <= time || time >= max_time) {
      if (return_luck) return List::create(Named("tte") = std::min(total_tte, max_time), Named("luck") = luck);
      else return List::create(Named("tte") = std::min(total_tte, max_time));
    }
  }
  
  return List::create(Named("tte") = NA_REAL);
}