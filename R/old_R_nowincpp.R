# Perform luck adjustment when using conditional quantile -------------------------------------------------------------------------------------------------------------------------------

#' Perform luck adjustment
#'
#' @param prevsurv Value of the previous survival
#' @param cursurv Value of the current survival
#' @param luck Luck used to be adjusted (number between 0 and 1)
#' @param condq Conditional quantile approach or standard approach
#'
#'
#' @return Adjusted luck number between 0 and 1
#'
#' @details
#' This function performs the luck adjustment automatically for the user, returning the adjusted luck number.
#' Luck is interpreted in the same fashion as is standard in R (higher luck, higher time to event).
#'
#' Note that if TTE is predicted using a conditional quantile function (e.g., conditional gompertz, conditional quantile weibull...) `prevsurv` and `cursurv`
#' are the unconditional survival using the "previous" parametrization but at the previous time for `presurv` and at the current time for `cursurv`.
#' For other distributions, `presurv` is the survival up to current time using the previous parametrization, and `cursurv`
#' is the survival up to current time using the current parametrization.
#'
#' Note that the advantage of the conditional quantile function is that it does not need the new parametrization to update the luck,
#' which makes this approach computationally more efficient.
#' This function can also work with vectors, which could allow to update multiple lucks in a single approach, and it can preserve names
#'
#' @examples
#' luck_adj_old(prevsurv = 0.8,
#'  cursurv = 0.7,
#'  luck = 0.5,
#'  condq = TRUE)
#'
#' luck_adj_old(prevsurv = c(1,0.8,0.7),
#'  cursurv = c(0.7,0.6,0.5),
#'  luck = setNames(c(0.5,0.6,0.7),c("A","B","C")),
#'  condq = TRUE)
#'
#' luck_adj_old(prevsurv = 0.8,
#'  cursurv = 0.7,
#'  luck = 0.5,
#'  condq = FALSE) #different results
#'
#' #Unconditional approach, timepoint of change is 25,
#' # parameter goes from 0.02 at time 10 to 0.025 to 0.015 at time 25,
#' #  starting luck is 0.37
#' new_luck <- luck_adj_old(prevsurv = 1 - pweibull(q=10,3,1/0.02),
#'  cursurv = 1 - pweibull(q=10,3,1/0.025),
#'  luck = 0.37,
#'  condq = FALSE) #time 10 change
#'
#' new_luck <- luck_adj_old(prevsurv = 1 - pweibull(q=25,3,1/0.025),
#'  cursurv = 1 - pweibull(q=25,3,1/0.015),
#'  luck = new_luck,
#'  condq = FALSE) #time 25 change
#'
#' qweibull(new_luck, 3, 1/0.015) #final TTE
#'
#' #Conditional quantile approach
#' new_luck <- luck_adj_old(prevsurv = 1-pweibull(q=0,3,1/0.02),
#'                       cursurv = 1- pweibull(q=10,3,1/0.02),
#'                       luck = 0.37,
#'                       condq = TRUE) #time 10 change, previous time is 0 so prevsurv will be 1
#'
#' new_luck <- luck_adj_old(prevsurv = 1-pweibull(q=10,3,1/0.025),
#'                       cursurv = 1- pweibull(q=25,3,1/0.025),
#'                       luck = new_luck,
#'                       condq = TRUE) #time 25 change
#'
#' qcond_weibull(rnd = new_luck,
#'                      shape = 3,
#'                      scale = 1/0.015,
#'                      lower_bound = 25) + 25 #final TTE
#' @keywords internal
#' @noRd
luck_adj_old <- function(prevsurv,cursurv,luck,condq=TRUE){
  
  #If length is 1, go for faster approach, otherwise use vectorial approach
  if(length(prevsurv)==1){
    if(condq==TRUE){
      adj_luck <-  max(1e-7, min(0.9999999, if(prevsurv != 0){ 1 - ((1 - luck) * (prevsurv/cursurv)) } else{ luck }))
    }else{
      adj_luck <-  max(1e-7, min(0.9999999, if(prevsurv != 0){ 1 - ((1 - luck) * (cursurv/prevsurv)) } else{ luck }))
    }
  } else{
    if(condq==TRUE){
      adj_luck <-  pmax(1e-7, pmin(0.9999999, ifelse(prevsurv != 0, 1 - ((1 - luck) * (prevsurv/cursurv)), luck )))
    }else{
      adj_luck <-   pmax(1e-7, pmin(0.9999999, ifelse(prevsurv != 0, 1 - ((1 - luck) * (cursurv/prevsurv)), luck )))
    }
  }
  
  #preserve names if existing
  names_luck <- names(luck)
  if(!is.null(names_luck)){
    names(adj_luck) <-names_luck
  }
  
  return(adj_luck)
}



#' Calculate discounted costs and qalys between events for vectors
#'
#' @param lcldr The discount rate 
#' @param lclprvtime The time of the previous event in the simulation
#' @param lclcurtime The time of the current event in the simulation
#' @param lclval The value to be discounted
#'
#' @return Double based on continuous time discounting
#'
#' @examples 
#' disc_ongoing_v_old(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)
#' 
#' @keywords internal
#' @noRd

disc_ongoing_v_old <- function(lcldr=0.035, lclprvtime, lclcurtime, lclval){
  
  Instantdr <- log(1+lcldr)
  
  # calculate additional qalys
  if(lcldr==0) {
    add <- lclval*(lclcurtime - lclprvtime)
  } else{
    add <- ((lclval)/(0 - Instantdr)) * (exp(lclcurtime * ( 0 - Instantdr)) - exp(lclprvtime * (0 - Instantdr)))
    
  }
  
  return(add)
}

#' Calculate instantaneous discounted costs or qalys for vectors
#'
#' @param lcldr The discount rate
#' @param lclcurtime The time of the current event in the simulation
#' @param lclval The value to be discounted
#'
#' @return Double based on discrete time discounting
#'
#' @examples
#' disc_instant_v_old(lcldr=0.035, lclcurtime=3, lclval=2500)
#' 
#'
#' @keywords internal
#' @noRd
disc_instant_v_old <- function(lcldr=0.035, lclcurtime, lclval){
  
  
  addinst <- lclval * ((1+lcldr)^(-lclcurtime))    
  
  return(addinst)
}


#' Cycle discounting for vectors
#'
#' @param lcldr The discount rate
#' @param lclprvtime The time of the previous event in the simulation
#' @param cyclelength The cycle length
#' @param lclcurtime The time of the current event in the simulation
#' @param lclval The  value to be discounted
#' @param starttime The start time for accrual of cycle costs (if not 0)
#' @param max_cycles The maximum number of cycles
#'
#' @return Double vector based on cycle discounting
#' 
#' @details
#' This function per cycle discounting, i.e., considers that the cost/qaly is accrued
#' per cycles, and performs it automatically without needing to create new events.
#' It can accommodate changes in cycle length/value/starttime (e.g., in the case of 
#' induction and maintenance doses) within the same item.
#'
#' @examples 
#' disc_cycle_v_old(lcldr=0.03, lclprvtime=0, cyclelength=1/12, lclcurtime=2, lclval=500,starttime=0)
#' disc_cycle_v_old(
#'  lcldr=0.000001,
#'  lclprvtime=0,
#'  cyclelength=1/12,
#'  lclcurtime=2,
#'  lclval=500,
#'  starttime=0,
#'  max_cycles = 4)
#' 
#' #Here we have a change in cycle length, max number of cylces and starttime at time 2
#'  #(e.g., induction to maintenance)
#' #In the model, one would do this by redifining cycle_l, max_cycles and starttime
#'  #of the corresponding item at a given event time. 
#' disc_cycle_v_old(lcldr=0,
#'  lclprvtime=c(0,1,2,2.5),
#'  cyclelength=c(1/12, 1/12,1/2,1/2),
#'  lclcurtime=c(1,2,2.5,4), lclval=c(500,500,500,500),
#'  starttime=c(0,0,2,2), max_cycles = c(24,24,2,2)
#'   )
#' @keywords internal
#' @noRd
disc_cycle_v_old <- function(
    lcldr = 0.035,
    lclprvtime = 0,
    cyclelength,
    lclcurtime,
    lclval,
    starttime = 0,
    max_cycles = NULL
) {
  n <- length(lclval)
  
  # Recycle scalar inputs if necessary
  if (length(lclprvtime) == 1) lclprvtime <- rep(lclprvtime, n)
  if (length(cyclelength) == 1) cyclelength <- rep(cyclelength, n)
  if (length(lclcurtime) == 1) lclcurtime <- rep(lclcurtime, n)
  if (length(starttime) == 1) starttime <- rep(starttime, n)
  if (is.null(max_cycles)) {
    max_cycles <- rep(NA, n)
  } else if (length(max_cycles) == 1) {
    max_cycles <- rep(max_cycles, n)
  }
  
  # Compute total possible cycles in this step
  total_cycles <- pmax(0, ceiling((lclcurtime - starttime) / cyclelength),na.rm=TRUE)
  prev_cycles <-  pmax(0, ceiling((lclprvtime - starttime) / cyclelength),na.rm=TRUE)
  n_cycles <- total_cycles - prev_cycles
  
  # Adjust n_cycles to account for the starting condition
  n_cycles <- ifelse(lclprvtime == lclcurtime & (lclcurtime - starttime) %% cyclelength == 0, n_cycles+1, n_cycles)
  
  # Cap using max_cycles
  remaining_cycles <- ifelse(is.na(max_cycles), n_cycles, pmax(0, max_cycles - prev_cycles))
  effective_cycles <- pmin(n_cycles, remaining_cycles)
  
  # Compute discount factor per row
  s <- (1 + lcldr)^cyclelength - 1
  d <- ifelse(lclprvtime==0,0,lclprvtime / cyclelength)-1
  
  discounted <- numeric(n)
  
  if (lcldr == 0) {
    discounted <- lclval * effective_cycles
  } else {
    discounted <- lclval * (1 - (1 + s)^(-effective_cycles)) / (s * (1 + s)^d)
  }
  
  return(discounted)
}


#' Quantile function for conditional Gompertz distribution (lower bound only)
#'
#' @param rnd Vector of quantiles
#' @param shape The shape parameter of the Gompertz distribution, defined as in the coef() output on a flexsurvreg object
#' @param rate The rate parameter of the Gompertz distribution, defined as in the coef() output on a flexsurvreg object
#' @param lower_bound The lower bound of the conditional distribution
#'
#' @return Estimate(s) from the conditional Gompertz distribution based on given parameters
#'
#' @export
#'
#' @examples
#' qcond_gompertz_old(rnd=0.5,shape=0.05,rate=0.01,lower_bound = 50)
#' @keywords internal
#' @noRd
qcond_gompertz_old <- function(rnd=0.5,shape,rate,lower_bound=0){
  if(rnd <0 | rnd > 1){
    stop("rnd is <0 or >1")
  }
  
  out <- (1/shape)*log(1- ((shape/rate)*log(1-rnd)/exp(shape*lower_bound)))
  
  return(out)
}


#' Conditional quantile function for exponential distribution 
#'
#' @param rnd Vector of quantiles
#' @param rate The rate parameter
#' 
#' Note taht the conditional quantile for an exponential is independent of time due to constant hazard
#'
#' @return Estimate(s) from the conditional exponential distribution based on given parameters
#'
#' @export
#'
#' @examples
#' qcond_exp_old(rnd = 0.5,rate = 3)
#' @keywords internal
#' @noRd
qcond_exp_old <- function(rnd = 0.5, rate) {
  if(rnd <0 | rnd > 1){
    stop("rnd is <0 or >1")
  }
  
  -log(1-rnd)/rate
}


#' Conditional quantile function for weibull distribution 
#'
#' @param rnd Vector of quantiles
#' @param shape The shape parameter as in R stats package weibull
#' @param scale The scale parameter as in R stats package weibull
#' @param lower_bound The lower bound to be used (current time)
#'
#' @return Estimate(s) from the conditional weibull distribution based on given parameters
#'
#' @export
#'
#' @examples
#' qcond_weibull_old(rnd = 0.5,shape = 3,scale = 66.66,lower_bound = 50)
#' @keywords internal
#' @noRd
qcond_weibull_old <- function(rnd = 0.5, shape, scale, lower_bound=0) {
  if(rnd <0 | rnd > 1){
    stop("rnd is <0 or >1")
  }
  
  ((lower_bound/scale)^shape - log(1-rnd))^(1/shape)*scale - lower_bound
}

#' Conditional quantile function for WeibullPH (flexsurv)
#'
#' @param rnd Vector of quantiles (between 0 and 1)
#' @param shape Shape parameter of WeibullPH
#' @param scale Scale (rate) parameter of WeibullPH (i.e., as in hazard = scale * t^(shape - 1))
#' @param lower_bound Lower bound (current time)
#'
#' @return Estimate(s) from the conditional weibullPH distribution based on given parameters
#'
#' @export
#'
#' @examples
#' qcond_weibullPH_old(rnd = 0.5, shape = 2, scale = 0.01, lower_bound = 5)
#' @keywords internal
#' @noRd
qcond_weibullPH_old <- function(rnd = 0.5, shape, scale, lower_bound = 0) {
  if (any(rnd < 0 | rnd > 1)) {
    stop("rnd must be between 0 and 1")
  }
  
  base <- lower_bound^shape - log(1 - rnd) / scale
  base[base < 0] <- NA  # to avoid complex numbers due to numerical issues
  
  residual <- base^(1 / shape) - lower_bound
  return(residual)
}


#' Conditional quantile function for loglogistic distribution 
#'
#' @param rnd Vector of quantiles
#' @param shape The shape parameter
#' @param scale The scale parameter
#' @param lower_bound The lower bound to be used (current time)
#'
#' @return Estimate(s) from the conditional loglogistic distribution based on given parameters
#'
#' @export
#'
#' @examples
#' qcond_llogis_old(rnd = 0.5,shape = 1,scale = 1,lower_bound = 1)
#' @keywords internal
#' @noRd
qcond_llogis_old <- function(rnd = 0.5, shape, scale, lower_bound=0) {
  if(rnd <0 | rnd > 1){
    stop("rnd is <0 or >1")
  }
  
  scale * ((1/(1-rnd))*(rnd + (lower_bound/scale)^shape))^(1/shape) - lower_bound
}

#' Conditional quantile function for lognormal distribution 
#'
#' @param rnd Vector of quantiles
#' @param meanlog The meanlog parameter
#' @param sdlog The sdlog parameter
#' @param lower_bound The lower bound to be used (current time)
#' @param s_obs is the survival observed up to lower_bound time,
#'  normally defined from time 0 as 1 - plnorm(q = lower_bound, meanlog, sdlog) but may be different if parametrization has changed previously
#'
#' @importFrom stats qnorm
#'
#' @return Estimate(s) from the conditional lognormal distribution based on given parameters
#'
#' @export
#'
#' @examples
#' qcond_lnorm_old(rnd = 0.5, meanlog = 1,sdlog = 1,lower_bound = 1, s_obs=0.8)
#' @keywords internal
#' @noRd
qcond_lnorm_old <- function(rnd = 0.5, meanlog, sdlog, lower_bound=0, s_obs) {
  if(rnd <0 | rnd > 1){
    stop("rnd is <0 or >1")
  }
  
  exp(meanlog + sdlog * qnorm(1 - s_obs * (1-rnd))) - lower_bound
}


#' Conditional quantile function for normal distribution 
#'
#' @param rnd Vector of quantiles
#' @param mean The mean parameter
#' @param sd The sd parameter
#' @param lower_bound The lower bound to be used (current time)
#' @param s_obs is the survival observed up to lower_bound time,
#'  normally defined from time 0 as 1 - pnorm(q = lower_bound, mean, sd) but may be different if parametrization has changed previously
#'
#' @importFrom stats qnorm
#'
#' @return Estimate(s) from the conditional normal distribution based on given parameters
#'
#' @export
#'
#' @examples
#' qcond_norm_old(rnd = 0.5, mean = 1,sd = 1,lower_bound = 1, s_obs=0.8)
#' @keywords internal
#' @noRd
qcond_norm_old <- function(rnd = 0.5, mean,sd, lower_bound=0, s_obs) {
  if(rnd <0 | rnd > 1){
    stop("rnd is <0 or >1")
  }
  
  qnorm(1 - s_obs*(1-rnd),mean,sd) - lower_bound
}


#' Conditional quantile function for gamma distribution 
#'
#' @param rnd Vector of quantiles
#' @param rate The rate parameter
#' @param shape The shape parameter
#' @param lower_bound The lower bound to be used (current time)
#' @param s_obs is the survival observed up to lower_bound time,
#'  normally defined from time 0 as 1 - pgamma(q = lower_bound, rate, shape) but may be different if parametrization has changed previously
#'
#' @importFrom stats qgamma
#'
#' @return Estimate(s) from the conditional gamma distribution based on given parameters
#'
#' @export
#'
#' @examples
#' qcond_gamma_old(rnd = 0.5, shape = 1.06178, rate = 0.01108,lower_bound = 1, s_obs=0.8)
#' @keywords internal
#' @noRd
qcond_gamma_old <- function(rnd = 0.5, shape,rate, lower_bound=0, s_obs) {
  if(rnd <0 | rnd > 1){
    stop("rnd is <0 or >1")
  }
  
  qgamma(1 - s_obs*(1-rnd),rate = rate, shape = shape) - lower_bound
}



#' Draw Time-to-Event with Time-Dependent Covariates and Luck Adjustment
#'
#' Simulate a time-to-event (TTE) from a parametric distribution with parameters varying over time.
#' User provides parameter functions and distribution name. The function uses internal survival and
#' conditional quantile functions, plus luck adjustment to simulate the event time. See
#' the vignette on avoiding cycles for an example in a model.
#' 
#'
#' @param luck Numeric between 0 and 1. Initial random quantile (luck).
#' @param a_fun Function of time .time returning the first distribution parameter (e.g., rate, shape, meanlog).
#' @param b_fun Function of time .time returning the second distribution parameter (e.g., scale, sdlog). Defaults to a function returning NA.
#' @param dist Character string specifying the distribution. Supported: "exp", "gamma", "lnorm", "norm", "weibull", "llogis", "gompertz".
#' @param dt Numeric. Time step increment to update parameters and survival. Default 0.1.
#' @param max_time Numeric. Max allowed event time to prevent infinite loops. Default 100.
#' @param return_luck Boolean. If TRUE, returns a list with tte and luck (useful if max_time caps TTE)
#' @param start_time Numeric. Time to use as a starting point of reference (e.g., curtime).
#'
#' @return Numeric. Simulated time-to-event.
#' 
#' 
#' @importFrom flexsurv pllogis pgompertz pweibullPH
#' @importFrom stats pexp pgamma plnorm pnorm pweibull
#' 
#' @export
#' 
#' @details The objective of this function is to avoid the user to have cycle events
#' with the only scope of updating some variables that depend on time and re-evaluate
#' a TTE. The idea is that this function should only be called at start and when 
#' an event impacts a variable (e.g., stroke event impacting death TTE), in which case
#' it would need to be called again at that point. In that case, the user would need to 
#' call e.g., `a <- qtimecov` with `max_time = curtime, return_luck = TRUE` arguments,
#' and then call it again with no max_time, return_luck = FALSE, and
#' `luck = a$luck, start_time=a$tte` (so there is no need to add curtime to the resulting time).
#' 
#' It's recommended to play with `dt` argument to balance running time and precision of the estimates.
#' For example, if we know we only update the equation annually (not continuously),
#' then we could just set `dt = 1`, which would make computations faster.
#'
#' @examples
#'
#' param_fun_factory <- function(p0, p1, p2, p3) {
#'   function(.time) p0 + p1*.time + p2*.time^2 + p3*(floor(.time) + 1)
#' }
#' 
#' set.seed(42)
#' 
#' # 1. Exponential Example
#' rate_exp <- param_fun_factory(0.1, 0, 0, 0)
#' qtimecov_old(
#'   luck = runif(1),
#'   a_fun = rate_exp,
#'   dist = "exp"
#' )
#' 
#' 
#' # 2. Gamma Example
#' shape_gamma <- param_fun_factory(2, 0, 0, 0)
#' rate_gamma <- param_fun_factory(0.2, 0, 0, 0)
#' qtimecov_old(
#'   luck = runif(1),
#'   a_fun = shape_gamma,
#'   b_fun = rate_gamma,
#'   dist = "gamma"
#' )
#' 
#' 
#' # 3. Lognormal Example
#' meanlog_lnorm <- param_fun_factory(log(10) - 0.5*0.5^2, 0, 0, 0)
#' sdlog_lnorm <- param_fun_factory(0.5, 0, 0, 0)
#' qtimecov_old(
#'   luck = runif(1),
#'   a_fun = meanlog_lnorm,
#'   b_fun = sdlog_lnorm,
#'   dist = "lnorm"
#' )
#' 
#' 
#' # 4. Normal Example
#' mean_norm <- param_fun_factory(10, 0, 0, 0)
#' sd_norm <- param_fun_factory(2, 0, 0, 0)
#' qtimecov_old(
#'   luck = runif(1),
#'   a_fun = mean_norm,
#'   b_fun = sd_norm,
#'   dist = "norm"
#' )
#' 
#' 
#' # 5. Weibull Example
#' shape_weibull <- param_fun_factory(2, 0, 0, 0)
#' scale_weibull <- param_fun_factory(10, 0, 0, 0)
#' qtimecov_old(
#'   luck = runif(1),
#'   a_fun = shape_weibull,
#'   b_fun = scale_weibull,
#'   dist = "weibull"
#' )
#' 
#' 
#' # 6. Loglogistic Example
#' shape_llogis <- param_fun_factory(2.5, 0, 0, 0)
#' scale_llogis <- param_fun_factory(7.6, 0, 0, 0)
#' qtimecov_old(
#'   luck = runif(1),
#'   a_fun = shape_llogis,
#'   b_fun = scale_llogis,
#'   dist = "llogis"
#' )
#' 
#' 
#' # 7. Gompertz Example
#' shape_gomp <- param_fun_factory(0.01, 0, 0, 0)
#' rate_gomp <- param_fun_factory(0.091, 0, 0, 0)
#' qtimecov_old(
#'   luck = runif(1),
#'   a_fun = shape_gomp,
#'   b_fun = rate_gomp,
#'   dist = "gompertz"
#' )
#'
#' #Time varying example, with change at time 8
#' rate_exp <- function(.time) 0.1 + 0.01*.time * 0.00001*.time^2
#' rate_exp2 <- function(.time) 0.2 + 0.02*.time
#' time_change <- 8
#' init_luck <- 0.95
#'
#' a <- qtimecov_old(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005,
#'                       max_time = time_change, return_luck = TRUE)
#' qtimecov_old(luck = a$luck,a_fun = rate_exp2,dist = "exp", dt = 0.005, start_time=a$tte)
#' 
#' 
#' #An example of how it would work in the model, this would also work with time varying covariates!
#' rate_exp <- function(.time) 0.1
#' rate_exp2 <- function(.time) 0.2
#' rate_exp3 <- function(.time) 0.3
#' time_change <- 10 #evt 1
#' time_change2 <- 15 #evt2
#' init_luck <- 0.95
#' #at start, we would just draw TTE
#' qtimecov_old(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005)
#' 
#' #at event in which rate changes (at time 10) we need to do this:
#' a <- qtimecov_old(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005,
#'                       max_time = time_change, return_luck = TRUE)
#' new_luck <- a$luck
#' qtimecov_old(luck = new_luck,a_fun = rate_exp2,dist = "exp", dt = 0.005, start_time=a$tte)
#' 
#' #at second  event in which rate changes again (at time 15) we need to do this:
#' a <- qtimecov_old(luck = new_luck,a_fun = rate_exp2,dist = "exp", dt = 0.005,
#'                       max_time = time_change2, return_luck = TRUE, start_time=a$tte)
#' new_luck <- a$luck
#' #final TTE is
#' qtimecov_old(luck = new_luck,a_fun = rate_exp3,dist = "exp", dt = 0.005, start_time=a$tte)
#'   
#' @keywords internal
#' @noRd
qtimecov <- function(
    luck,
    a_fun,
    b_fun = function(.time) NA,
    dist,
    dt = 0.1,
    max_time = 100,
    return_luck = FALSE,
    start_time = 0
) {
  # Internal survival and conditional quantile functions (require flexsurv for some)
  sf_fun <- switch(dist,
                   exp       = function(.time, rate, ...) 1 - pexp(.time, rate),
                   gamma     = function(.time, shape, rate) 1 - pgamma(.time, shape = shape, rate = rate),
                   lnorm     = function(.time, meanlog, sdlog) 1 - plnorm(.time, meanlog, sdlog),
                   norm      = function(.time, mean, sd) 1 - pnorm(.time, mean, sd),
                   weibull   = function(.time, shape, scale) 1 - pweibull(.time, shape, scale),
                   weibullPH = function(.time, shape, scale) 1 - flexsurv::pweibullPH(.time, shape, scale),
                   llogis    = function(.time, shape, scale) flexsurv::pllogis(q = .time, shape = shape, scale = scale, lower.tail = FALSE),
                   gompertz  = function(.time, shape, rate) flexsurv::pgompertz(q = .time, shape = shape, rate = rate, lower.tail = FALSE),
                   stop("Unsupported distribution")
  )
  
  qcond_fun <- switch(dist,
                      exp       = qcond_exp,
                      gamma     = qcond_gamma,
                      lnorm     = qcond_lnorm,
                      norm      = qcond_norm,
                      weibull   = qcond_weibull,
                      weibullPH = qcond_weibullPH,
                      llogis    = qcond_llogis,
                      gompertz  = qcond_gompertz,
                      stop("Unsupported distribution")
  )
  
  .time <- start_time
  a_curr <- a_fun(.time)
  b_curr <- b_fun(.time)
  surv_prev <- sf_fun(.time, a_curr, b_curr)
  
  
  # Early exit: check whether event occurs in [0, dt]
  residual_tte <- switch(dist,
                         exp       = qcond_exp(luck, rate = a_curr),
                         gamma     = qcond_gamma(luck, rate = b_curr, shape = a_curr, lower_bound = start_time, s_obs = surv_prev),
                         lnorm     = qcond_lnorm(luck, meanlog = a_curr, sdlog = b_curr, lower_bound = start_time, s_obs = surv_prev),
                         norm      = qcond_norm(luck, mean = a_curr, sd = b_curr, lower_bound = start_time, s_obs = surv_prev),
                         weibull   = qcond_weibull(luck, shape = a_curr, scale = b_curr, lower_bound = start_time),
                         weibullPH = qcond_weibullPH(luck, shape = a_curr, scale = b_curr, lower_bound = start_time),
                         llogis    = qcond_llogis(luck, shape = a_curr, scale = b_curr, lower_bound = start_time),
                         gompertz  = qcond_gompertz(luck, shape = a_curr, rate = b_curr, lower_bound = start_time),
                         stop("Unsupported distribution")
  )
  
  if (residual_tte <= dt) {
    if (return_luck == TRUE) {
      return(list(tte = residual_tte + .time, luck = luck))
    } else {
      return(residual_tte)
    }
  }
  
  repeat {
    .time <- .time + dt
    
    # Correctly compute survival at .time - dt for luck adjustment
    surv_prev <- sf_fun(.time - dt, a_curr, b_curr)
    surv_curr <- sf_fun(.time, a_curr, b_curr)
    a_curr <- a_fun(.time)
    b_curr <- b_fun(.time)
    
    # Use condq = TRUE for conditional quantiles
    luck <- luck_adj(prevsurv = surv_prev, cursurv = surv_curr, luck = luck, condq = TRUE)
    
    residual_tte <- switch(dist,
                           exp       = qcond_exp(luck, rate = a_curr),
                           gamma     = qcond_gamma(luck, rate = b_curr, shape = a_curr, lower_bound = .time - dt, s_obs = surv_prev),
                           lnorm     = qcond_lnorm(luck, meanlog = a_curr, sdlog = b_curr, lower_bound = .time - dt, s_obs = surv_prev),
                           norm      = qcond_norm(luck, mean = a_curr, sd = b_curr, lower_bound = .time - dt, s_obs = surv_prev),
                           weibull   = qcond_weibull(luck, shape = a_curr, scale = b_curr, lower_bound = .time - dt),
                           weibullPH = qcond_weibullPH(luck, shape = a_curr, scale = b_curr, lower_bound = .time - dt),
                           llogis    = qcond_llogis(luck, shape = a_curr, scale = b_curr, lower_bound = .time - dt),
                           gompertz  = qcond_gompertz(luck, shape = a_curr, rate = b_curr, lower_bound = .time - dt),
                           stop("Unsupported distribution")
    )
    
    total_tte <- .time + residual_tte
    
    if (residual_tte <= dt || total_tte <= .time || .time >= max_time) {
      if (return_luck == TRUE) {
        return(list(tte = min(total_tte, max_time), luck = luck))
      } else {
        return(min(total_tte, max_time))
      }
    }
  }
}


