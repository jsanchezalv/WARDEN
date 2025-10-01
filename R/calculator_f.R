#The functions here are not worth to convert to C++ based on internal experiments
#except for rpoisgamma, which was converted as rpoisgamma_rcpp, but given its complexity,
#the R version is left as well (but it's 5-6x slower)

# Global variables for CRAN check -----------------------------------------

if(getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      c("t_btw_evt", #rpoisgamma
        "tte",
        "evt_num",
        ".id",
        "evt_count")
    )) 
}

# Draw TTE -------------------------------------------------------------------------

#' Draw a time to event from a list of parametric survival functions
#'
#' @param n_chosen The number of observations to be drawn
#' @param dist The distribution; takes values 'lnorm','norm','mvnorm','weibullPH','weibull','llogis','gompertz','gengamma','gamma','exp','beta','poisgamma'
#' @param coef1 First coefficient of the distribution, defined as in the coef() output on a flexsurvreg object (rate in "rpoisgamma")
#' @param coef2 Second coefficient of the distribution, defined as in the coef() output on a flexsurvreg object (theta in "rpoisgamma")
#' @param coef3 Third coefficient of the distribution, defined as in the coef() output on a flexsurvreg object (not used in "rpoisgamma")
#' @param beta_tx Parameter in natural scale applied in addition to the scale/rate coefficient -e.g., a HR if used in an exponential- (not used in "rpoisgamma" nor "beta")
#' @param seed An integer which will be used to set the seed for this draw.
#' @param ... Additional arguments to be used by the specific distribution (e.g., return_ind_rate if dist = "poisgamma")
#'
#' @return A vector of time to event estimates from the given parameters
#'
#' @importFrom stats rlnorm rnorm rweibull rgamma rexp rbeta
#' @importFrom flexsurv rweibullPH rllogis rgompertz rgengamma
#' @importFrom MASS mvrnorm
#'
#' @export
#' 
#' @details Other arguments relevant to each function can be called directly
#'
#' @examples
#' draw_tte(n_chosen=1,dist='exp',coef1=1,beta_tx=1)
#' draw_tte(n_chosen=10,"poisgamma",coef1=1,coef2=1,obs_time=1,return_ind_rate=FALSE)

draw_tte <- function(n_chosen,dist,coef1=NULL,coef2=NULL,coef3=NULL,...,beta_tx=1,seed=NULL) {

  if(!is.null(seed)){
    set.seed(seed)
  }

  if (any(length(coef1)>1,length(coef2)>1,length(coef3)>1) & dist!="mvnorm") {
    message("Provided a coefficient parameter that is a vector")
  }
  
  if (!(dist %in% c("lnorm", "norm", "mvnorm", "weibullPH", "weibull", "llogis", "gompertz", "gengamma", "gamma", "exp", "beta", "poisgamma"))) {
    stop("Invalid distribution. Distribution must be one of: 'lnorm','norm','mvnorm','weibullPH','weibull','llogis','gompertz','gengamma','gamma','exp','beta','poisgamma'")
  }

  draw.out <- switch(dist, 
                     "lnorm" = rlnorm(n_chosen, meanlog=coef1 - log(beta_tx), sdlog=exp(coef2),...), 
                     "norm" = rnorm(n_chosen, mean=coef1 - beta_tx, sd=coef2,...), 
                     "mvnorm" = mvrnorm(n_chosen, mu=coef1 - beta_tx, Sigma=coef2,...), 
                     "weibullPH" = rweibullPH(n_chosen, shape = exp(coef1), scale = exp(coef2 + log(beta_tx))),
                     "weibull" = rweibull(n_chosen, shape = exp(coef1), scale = exp(coef2 + log(beta_tx))),
                     "llogis" = rllogis(n_chosen, shape = exp(coef1), scale = exp(coef2 + log(beta_tx))),
                     "gompertz" = rgompertz(n_chosen, shape = coef1, rate = exp(coef2 + log(beta_tx))),
                     "gengamma" = rgengamma(n_chosen, mu = coef1 + log(beta_tx), sigma = exp(coef2), Q = coef3),
                     "gamma" = rgamma(n_chosen, shape = exp(coef1), rate = exp(coef2 + log(beta_tx))),
                     "exp" = rexp(n_chosen, rate = exp(coef1 + log(beta_tx))),
                     "beta" = rbeta(n_chosen, shape1 = coef1, shape2 = coef2),
                     "poisgamma" = rpoisgamma(n_chosen, rate = coef1, theta=coef2,...),
                     stop("Invalid distribution. Distribution must be one of: 'lnorm','weibullPH','weibull','llogis','gompertz','gengamma','gamma','exp','beta','poisgamma'")
  )
  
  return(draw.out)
}

# Draws dirichlet distribution --------

#' Draw from a dirichlet distribution based on number of counts in transition. Adapted from brms::rdirichlet
#'
#' @param n Number of draws (must be >= 1). If n>1, it will return a list of matrices.
#' @param alpha A matrix of alphas (transition counts)
#' @param seed An integer which will be used to set the seed for this draw.
#'
#' @return A transition matrix. If n>1, it will return a list of matrices.
#'
#' @importFrom stats rgamma
#'
#' @export
#'
#' @examples
#' rdirichlet(n=1,alpha= matrix(c(1251, 0, 350, 731),2,2))
#' rdirichlet(n=2,alpha= matrix(c(1251, 0, 350, 731),2,2))

rdirichlet <- function(n=1,alpha,seed=NULL) {
  out <- NULL
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
    if (!is.matrix(alpha)) {
      alpha <- matrix(alpha, nrow = 1)
    }
    if (prod(dim(alpha)) == 0) {
      stop("alpha should be non-empty.")
    }
    if (isTRUE(any(alpha < 0))) {
      stop("alpha must be positive.")
    }
    if (n == 1) {
      out <- matrix(rgamma(ncol(alpha) * nrow(alpha), alpha), ncol = ncol(alpha))
      out <- out/rowSums(out)
      
    } else{
     out <-  lapply(1:n,function(x) {
      out <- matrix(rgamma(ncol(alpha) * nrow(alpha), alpha), ncol = ncol(alpha)) 
      out/rowSums(out)} 
      )
    }
  
  return(out)
}

# Draws dirichlet distribution --------

#' Draw from a dirichlet distribution based on mean transition probabilities and standard errors
#'
#' @param n Number of draws (must be >= 1). If n>1, it will return a list of matrices.
#' @param alpha A matrix of transition probabilities
#' @param se A matrix of standard errors 
#' @param seed An integer which will be used to set the seed for this draw.
#'
#' @return A transition matrix. If n>1, it will return a list of matrices.
#'
#' @importFrom stats rgamma
#'
#' @export
#'
#' @examples
#' rdirichlet_prob(n=1,alpha= matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7),3,3),
#' se=matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7)/10,3,3))
#' 
#' rdirichlet_prob(n=2,alpha= matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7),3,3),
#' se=matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7)/10,3,3))

rdirichlet_prob <- function(n=1,alpha,se,seed=NULL) {
  out <- NULL
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if (!is.matrix(alpha)) {
    alpha <- matrix(alpha, nrow = 1)
  }
  if (prod(dim(alpha)) == 0) {
    stop("alpha should be non-empty.")
  }
  if (isTRUE(any(alpha < 0))) {
    stop("alpha must be positive.")
  }
  if (n == 1) {
    out <- matrix(suppressWarnings(rgamma_mse(1, c(alpha),c(se))), ncol = ncol(alpha))
    out <- out/rowSums(out)
    
  } else{
    out <-  lapply(1:n,function(x) {
      out <- matrix(suppressWarnings(rgamma_mse(1, c(alpha),c(se))), ncol = ncol(alpha)) 
      out/rowSums(out)} 
    )
  }
  
  return(out)
}


# Draws beta distribution - accepts vectors--------

#' Draw from a beta distribution based on mean and se
#'
#' @param n Number of draws (must be >= 1)
#' @param mean_v A vector of the mean values
#' @param se A vector of the standard errors of the means
#' @param seed An integer which will be used to set the seed for this draw.
#'
#' @return A single estimate from the beta distribution based on given parameters
#'
#' @importFrom stats rbeta
#'
#' @export
#'
#' @examples
#' rbeta_mse(n=1,mean_v=0.8,se=0.2)

rbeta_mse <- function(n=1,mean_v,se,seed=NULL) {

  if(!is.null(seed)){
    set.seed(seed)
  }

    alpha <- ((1 - mean_v) / (se^2) - (1 / mean_v)) * mean_v ^ 2
    beta <- alpha * ((1 / mean_v) - 1)
    out <- rbeta(length(mean_v)*n,alpha,beta)

  return(out)
}


#' Draw from a beta distribution based on mean and se (quantile)
#'
#' @param q Quantiles to be used
#' @param mean_v A vector of the mean values
#' @param se A vector of the standard errors of the means
#'
#' @return A single estimate from the beta distribution based on given parameters
#'
#' @importFrom stats qbeta
#'
#' @export
#'
#' @examples
#' qbeta_mse(q=0.5,mean_v=0.8,se=0.2)

qbeta_mse <- function(q,mean_v,se) {
  alpha <- ((1 - mean_v) / (se^2) - (1 / mean_v)) * mean_v ^ 2
  beta <- alpha * ((1 / mean_v) - 1)
  out <- qbeta(q,alpha,beta)
  
  return(out)
}

# Draws gamma distribution - accepts vectors--------

#' Draw from a gamma distribution based on mean and se
#'
#' @param n Number of draws (must be >= 1)
#' @param mean_v A vector of the mean values
#' @param se A vector of the standard errors of the means
#' @param seed An integer which will be used to set the seed for this draw.
#'
#' @return A single estimate from the gamma distribution based on given parameters
#'
#' @importFrom stats rgamma
#'
#' @export
#'
#' @examples
#' rgamma_mse(n=1,mean_v=0.8,se=0.2)
#'

rgamma_mse <- function(n=1, mean_v, se, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  len <- length(mean_v)
  scale <- se^2 / mean_v
  shape <- mean_v / scale
  out <- stats::rgamma(len * n, shape, scale = scale)
  det <- (se == 0 | mean_v == 0)
  det[is.na(det)] <- FALSE
  if (any(det)) {
    for (i in which(det)) {
      out[seq.int(i, len * n, by = len)] <- mean_v[i]
    }
  }
  out
}

# Quantile gamma distribution - accepts vectors--------

#' Use quantiles from a gamma distribution based on mean and se
#'
#' @param q Quantile to draw
#' @param mean_v A vector of the mean values
#' @param se A vector of the standard errors of the means
#' @param seed An integer which will be used to set the seed for this draw.
#'
#' @return A single estimate from the gamma distribution based on given parameters
#'
#' @importFrom stats rgamma
#'
#' @export
#'
#' @examples
#' qgamma_mse(q=0.5,mean_v=0.8,se=0.2)
#'
qgamma_mse <- function(q=1, mean_v, se, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  len <- length(mean_v)
  scale <- se^2 / mean_v
  shape <- mean_v / scale
  out <- stats::qgamma(q, shape, scale = scale)
  det <- (se == 0 | mean_v == 0)
  det[is.na(det)] <- FALSE
  out[det] <- mean_v[det]
  out
}

#' Draw from a Conditional Gompertz distribution (lower and upper bound)
#'
#' @param n The number of observations to be drawn
#' @param shape The shape parameter of the Gompertz distribution, defined as in the coef() output on a flexsurvreg object
#' @param rate The rate parameter of the Gompertz distribution, defined as in the coef() output on a flexsurvreg object
#' @param lower_bound The lower bound of the conditional distribution
#' @param upper_bound The upper bound of the conditional distribution
#' @param seed An integer which will be used to set the seed for this draw.
#'
#' @return Estimate(s) from the Conditional Gompertz distribution based on given parameters
#'
#' @importFrom stats runif
#'
#' @export
#'
#' @examples
#' rcond_gompertz_lu(1,shape=0.05,rate=0.01,lower_bound = 50)

rcond_gompertz_lu <- function(n, shape, rate , lower_bound = 0, upper_bound = Inf, seed=NULL){

  if(!is.null(seed)){
    set.seed(seed)
  }

  surv_prob <- flexsurv::pgompertz(c(lower_bound, upper_bound),shape, rate)
  uniform_random_numbers <- stats::runif(n, surv_prob[1], surv_prob[2])
  flexsurv::qgompertz(uniform_random_numbers, shape, rate ) - lower_bound
}

#' Draw from a conditional Gompertz distribution (lower bound only)
#'
#' @param n The number of observations to be drawn
#' @param shape The shape parameter of the Gompertz distribution, defined as in the coef() output on a flexsurvreg object
#' @param rate The rate parameter of the Gompertz distribution, defined as in the coef() output on a flexsurvreg object
#' @param lower_bound The lower bound of the conditional distribution
#' @param seed An integer which will be used to set the seed for this draw.
#'
#' @return Estimate(s) from the conditional Gompertz distribution based on given parameters
#'
#' @importFrom stats runif
#'
#' @export
#'
#' @examples
#' rcond_gompertz(1,shape=0.05,rate=0.01,lower_bound = 50)

rcond_gompertz <- function(n=1,shape,rate,lower_bound=0,seed=NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  out <- (1/shape)*log(1- ((shape/rate)*log(runif(n))/exp(shape*lower_bound)))
  
  return(out)
}



#' Survival Probaility function for conditional Gompertz distribution (lower bound only)
#'
#' @param time Vector of times
#' @param shape The shape parameter of the Gompertz distribution, defined as in the coef() output on a flexsurvreg object
#' @param rate The rate parameter of the Gompertz distribution, defined as in the coef() output on a flexsurvreg object
#' @param lower_bound The lower bound of the conditional distribution
#'
#' @return Estimate(s) from the conditional Gompertz distribution based on given parameters
#'
#'
#' @export
#'
#' @examples
#' pcond_gompertz(time=1,shape=0.05,rate=0.01,lower_bound = 50)

pcond_gompertz <- function(time=1,shape,rate,lower_bound=0){
  
  out <- 1- exp(-(rate/shape) * exp(shape*lower_bound)*(exp(shape * time)-1))
    
  return(out)
}

#' Draw time to event (tte) from a Poisson or Poisson-Gamma (PG) Mixture/Negative Binomial (NB) Process
#'
#' @param n The number of observations to be drawn
#' @param rate rate of the event (in terms of events per observation-time)
#' @param obs_time period over which events are observable
#' @param theta Optional.  When omitted, the function simulates times for a Poisson process.
#'               Represents the shape of the gamma mixture distribution. 
#'               Estimated and reported as theta in negative binomial regression analyses in r. 
#' @param t_reps Optional. Number of TBEs to be generated to capture events within the observation window.  
#'               When omitted, the function sets t_reps to the 99.99th quantile of the Poisson (if no theta is provided)
#'               or negative binomial (if theta is provided). Thus, the risk of missing possible events in the observation window
#'               is 0.01%.
#' @param seed An integer which will be used to set the seed for this draw.
#' @param return_ind_rate A boolean that indicates whether an additional vector with the rate parameters used per observation is used.
#'          It will alter the structure of the results to two lists, one storing tte with name tte, and the other with name ind_rate
#' @param return_df A boolean that indicates whether a data.table object should be returned
#'
#' @return Estimate(s) from the time to event based on poisson/Poisson-Gamma (PG) Mixture/Negative Binomial (NB) distribution based on given parameters
#'
#' @importFrom stats qpois
#' @importFrom stats rexp
#' @importFrom stats rgamma
#' @importFrom stats qnbinom
#' @importFrom data.table rbindlist
#'
#' @export
#' 
#' @details
#' Function to simulate event times from a Poisson or Poisson-Gamma (PG) Mixture/Negative Binomial (NB) Process
#' Event times are determined by sampling times between events (TBEs) from an exponential distribution, and cumulating 
#' these to derive the event times. Events occurring within the set observation time window are retained and returned.
#' For times for a Poisson process, the provided rate is assumed constant.
#' For a PG or NB, the individual rates are sampled from a Gamma distribution with shape = theta and scale = rate/theta.
#'
#' @examples
#' rpoisgamma(1,rate=1,obs_time=1,theta=1)

rpoisgamma <- function(n, rate, theta=NULL, obs_time=1, t_reps, seed=NULL,return_ind_rate=FALSE, return_df=FALSE){
  # Create data with sampled event times for n observations and t_reps replications                      
  # Approach is different for Poisson and PG to optimize run time
  # If t_reps not provided, derive based on 99.9th quantile of Poisson or negative binomial distribution

  if(is.null(theta)){
    # Missing theta produces event time draws for a Poisson process
    if (missing(t_reps)) t_reps <- qpois(0.9999, lambda=rate*obs_time)
    # Determine the number of time replications needed                       
    # Based on 99.99th quantile of Poisson  distribution 
    
    #Draw the relevant tte that we will after filter for those which actually occured
    
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    ds <- lapply(1:n, function(x) cumsum(rexp(t_reps, rate=rate))) 
    
  } else{
    # Theta specified; produce event times for Poisson-gamma/NB
    
    if (missing(t_reps)) t_reps <- qnbinom(0.9999, size=theta, mu=rate*obs_time)
    # Determine the number of time replications needed                       
    # Based on 99.99th quantile of NB  distribution 
    
    # For Poisson-Gamma, individual rates are first sampled from a Gamma distribution
      rate_par <- rgamma(n = n, shape = theta, scale = rate/theta)
      
      if(!is.null(seed)){
        set.seed(seed)
      }
      
      ds <- lapply(1:n, function(x) 
        cumsum(
            rexp(t_reps,
               rate=rate_par[x]
          )
        )
      )
  }

  # Process sampled times between events to derive actual event times (as cumulative of time between events)
  # Determine which are observable - i.e., within the obs_time set in argument
  # Flag observable events and derive count of observable events per individual
  # Retain only observable events 
  
  #Get which observations occur within time period and remove unobserved events
  ds <- lapply(ds, function(x) x[x<=obs_time])

  
  #If return as data frame is activated
  if (return_df==TRUE) {
    ds_df <- lapply(ds, function(x) data.frame(tte=x))
    ds_df <- rbindlist(ds_df,idcol=TRUE)
    if(nrow(ds_df)==0){
      return(NULL)
    }
    ds_df[,t_btw_evt:= tte-shift(tte,fill = 0),by=.id]
    ds_df[,evt_num:= 1:.N,by=.id]
    ds_df[,evt_count:= .N,by=.id]
    
    if(return_ind_rate==TRUE){
    rate_par_df <- data.table(.id=1:length(rate_par),
                              ind_rate = rate_par)
    
    ds_df <- merge(ds_df,rate_par_df)
    }
    
    return(ds_df)
    
    #If the rate parameters have to be returned, modify the return structure to accommodate them
  } else{
    if(return_ind_rate==TRUE){
      return(list(tte=ds,ind_rate=rate_par))
    }  else{
      return(ds)
    }
  }
}

#' Calculate conditional multivariate normal values
#'
#' @param mu mean vector
#' @param Sigma covariance matrix
#' @param i index of the known parameter (1-based index)
#' @param xi known value of the i-th parameter
#' @param full_output boolean indicating whether to return the full list of parameters
#' 
#' @return List of length 2, one with new mu and other with covariance parameters
#'
#' @export
#' 
#' @details
#'Function to compute conditional multivariate normal values
#'
#' @examples
#' mu <- c(1, 2, 3)
#' Sigma <- matrix(c(0.2, 0.05, 0.1, 
#'                   0.05, 0.3, 0.05, 
#'                   0.1, 0.05, 0.4), nrow = 3)
#' 
#' i <- 1:2  # Index of the known parameter
#' xi <- c(1.2,2.3)  # Known value of the first parameter
#'
#'cond_mvn(mu, Sigma, i, xi,full_output = TRUE)

cond_mvn <- function(mu, Sigma, i, xi, full_output = FALSE) {
  
  if(length(mu)==length(i)){stop("Trying to condition on all parameters")}
  
  N <- length(mu)
  
  # Create index vectors
  all_indices <- 1:N
  not_i_indices <- all_indices[-i]
  
  # Extract necessary submatrices and subvectors
  mu_i <- mu[i]
  mu_not_i <- mu[not_i_indices]
  
  Sigma_ii <- Sigma[i, i]
  Sigma_i_not_i <- Sigma[i, not_i_indices]
  Sigma_not_i_i <- Sigma[not_i_indices, i]
  Sigma_not_i_not_i <- Sigma[not_i_indices, not_i_indices]
  
  # Compute the conditional mean and covariance
  mu_cond <- mu_not_i + Sigma_not_i_i %*% solve(Sigma_ii) %*% (xi - mu_i)
  Sigma_cond <- Sigma_not_i_not_i - Sigma_not_i_i %*% solve(Sigma_ii) %*% Sigma_i_not_i
  
  if (full_output) {
    # Create the full mean vector including the known value
    mu_full <- mu
    mu_full[not_i_indices] <- mu_cond
    mu_full[i] <- xi
    
    # Create the full mean vector including the known value
    Sigma_full <- Sigma
    Sigma_full[not_i_indices, not_i_indices] <- Sigma_cond
    Sigma_full[i, ] <- 0
    Sigma_full[, i] <- 0
    Sigma_full[i, i] <- 0  # Variance of the known value is zero
    
    # Return the full mean vector and conditional covariance matrix
    return(list(mean = mu_full, covariance = Sigma_full))
  } else {
    # Return only the conditional mean and covariance matrix
    return(list(mean = mu_cond, covariance = Sigma_cond))
  }
}

#' Calculate conditional dirichlet values
#'
#' @param alpha mean vector
#' @param i index of the known parameter (1-based index)
#' @param xi known value of the i-th parameter (should be >0)
#' @param full_output boolean indicating whether to return the full list of parameters
#' 
#' @return List of length 2, one with new mu and other with covariance parameters
#'
#' @export
#' 
#' @details
#'Function to compute conditional dirichlet values
#'
#' @examples
#' alpha <- c(2, 3, 4)
#' i <- 2  # Index of the known parameter
#' xi <- 0.5  # Known value of the second parameter
#' 
#' # Compute the conditional alpha parameters with full output
#' cond_dirichlet(alpha, i, xi, full_output = TRUE)

cond_dirichlet <- function(alpha, i, xi, full_output = FALSE) {
  
  if(length(alpha)==length(i)){stop("Trying to condition on all parameters")}
  if(sum(xi) > 1){stop("Sum of xi should be <= 1")}
  if(any(xi < 0)){stop("At least one xi is negative, and should be >=0")}
  
  # Create index vectors
  all_indices <- 1:length(alpha)
  not_i_indices <- all_indices[-i]
  
  # Remove the i-th alpha parameter
  alpha_not_i <- alpha[not_i_indices]
  
  # Adjust the remaining alpha parameters
  alpha_cond <- alpha_not_i * (1 - sum(xi)) / sum(alpha_not_i)
  
  if (full_output) {
    # Create the full alpha vector including the known value
    alpha_full <- numeric(length(alpha))
    alpha_full[not_i_indices] <- alpha_cond
    alpha_full[i] <- xi
    
    # Return the full alpha vector
    return(alpha_full)
  } else {
    # Return only the conditional alpha parameters
    return(alpha_cond)
  }
}


.resource_discrete_from_xptr <- function(xptr) {
  env <- new.env(parent = emptyenv())
  env$.ptr <- xptr
  
  env$size                 <- function()        discrete_resource_size_cpp(env$.ptr)
  env$queue_size           <- function()        discrete_resource_queue_size_cpp(env$.ptr)
  env$n_free               <- function()        discrete_resource_n_free_cpp(env$.ptr)
  env$patients_using       <- function()        discrete_resource_patients_using_cpp(env$.ptr)
  env$patients_using_times <- function()        discrete_resource_patients_using_times_cpp(env$.ptr)
  env$queue_start_times    <- function()        discrete_resource_queue_start_times_cpp(env$.ptr)
  env$queue_priorities     <- function()        discrete_resource_queue_priorities_cpp(env$.ptr)
  
  env$queue_info <- function(n = NULL) {
    if (is.null(n)) n <- env$queue_size()
    if (n <= 0) return(data.frame(patient_id = integer(0), priority = integer(0), start_time = numeric(0)))
    ids  <- env$next_patient_in_line(n)
    prio <- env$queue_priorities()[seq_len(length(ids))]
    st   <- env$queue_start_times()[seq_len(length(ids))]
    data.frame(patient_id = ids, priority = prio, start_time = st, stringsAsFactors = FALSE)
  }
  
  env$is_patient_in_queue <- function(patient_id) {
    if (!is.numeric(patient_id) || length(patient_id) != 1) stop("patient_id must be a single number")
    discrete_resource_is_patient_in_queue_cpp(env$.ptr, as.integer(patient_id))
  }
  env$is_patient_using <- function(patient_id) {
    if (!is.numeric(patient_id) || length(patient_id) != 1) stop("patient_id must be a single number")
    discrete_resource_is_patient_using_cpp(env$.ptr, as.integer(patient_id))
  }
  
  env$attempt_block <- function(patient_id = NULL, priority = 1L, start_time = NULL) {
    if (is.null(patient_id)) {
      patient_id <- tryCatch(get("i", envir = parent.frame(), inherits = TRUE),
                             error = function(e) stop("patient_id not provided and 'i' not found in parent frame"))
    }
    if (is.null(start_time)) {
      start_time <- tryCatch(get("curtime", envir = parent.frame(), inherits = TRUE),
                             error = function(e) stop("start_time not provided and 'curtime' not found in parent frame"))
    }
    if (!is.numeric(patient_id) || length(patient_id) != 1) stop("patient_id must be a single number")
    if (!is.numeric(priority)   || length(priority)   != 1) stop("priority must be a single number")
    if (!is.numeric(start_time) || length(start_time) != 1) stop("start_time must be a single number")
    discrete_resource_attempt_block_cpp(env$.ptr, as.integer(patient_id), as.integer(priority), as.numeric(start_time))
  }
  
  env$attempt_free <- function(patient_id = NULL, remove_all = FALSE) {
    if (is.null(patient_id)) {
      patient_id <- tryCatch(get("i", envir = parent.frame(), inherits = TRUE),
                             error = function(e) stop("patient_id not provided and 'i' not found in parent frame"))
    }
    if (!is.numeric(patient_id) || length(patient_id) != 1) stop("patient_id must be a single number")
    if (!is.logical(remove_all) || length(remove_all) != 1) stop("remove_all must be a single logical value")
    discrete_resource_attempt_free_cpp(env$.ptr, as.integer(patient_id), remove_all)
    invisible(NULL)
  }
  
  env$attempt_free_if_using <- function(patient_id = NULL, remove_all = FALSE) {
    if (is.null(patient_id)) {
      patient_id <- tryCatch(get("i", envir = parent.frame(), inherits = TRUE),
                             error = function(e) stop("patient_id not provided and 'i' not found in parent frame"))
    }
    if (!is.numeric(patient_id) || length(patient_id) != 1) stop("patient_id must be a single number")
    if (!is.logical(remove_all) || length(remove_all) != 1) stop("remove_all must be a single logical value")
    discrete_resource_attempt_free_if_using_cpp(env$.ptr, as.integer(patient_id), remove_all)
    invisible(NULL)
  }
  
  env$next_patient_in_line <- function(n = 1L) {
    if (!is.numeric(n) || length(n) != 1 || n < 1) stop("n must be a single positive integer")
    discrete_resource_next_patient_in_line_cpp(env$.ptr, as.integer(n))
  }
  
  env$modify_priority <- function(patient_id, new_priority) {
    if (!is.numeric(patient_id) || length(patient_id) != 1) stop("patient_id must be a single number")
    if (!is.numeric(new_priority) || length(new_priority) != 1) stop("new_priority must be a single number")
    discrete_resource_modify_priority_cpp(env$.ptr, as.integer(patient_id), as.integer(new_priority))
    invisible(NULL)
  }
  
  env$add_resource <- function(n_to_add) {
    if (!is.numeric(n_to_add) || length(n_to_add) != 1 || n_to_add < 1) stop("n_to_add must be a single positive integer")
    discrete_resource_add_resource_cpp(env$.ptr, as.integer(n_to_add))
    invisible(NULL)
  }
  
  env$remove_resource <- function(n_to_remove, current_time = NULL) {
    if (is.null(current_time)) {
      current_time <- tryCatch(get("curtime", envir = parent.frame(), inherits = TRUE),
                               error = function(e) stop("current_time not provided and 'curtime' not found in parent frame"))
    }
    if (!is.numeric(n_to_remove) || length(n_to_remove) != 1 || n_to_remove < 1) stop("n_to_remove must be a single positive integer")
    if (!is.numeric(current_time) || length(current_time) != 1) stop("current_time must be a single number")
    discrete_resource_remove_resource_cpp(env$.ptr, as.integer(n_to_remove), as.numeric(current_time))
    invisible(NULL)
  }
  
  class(env) <- "resource_discrete"
  env
}

# Exported: fast, always returns a list of wrapper envs (even for n = 1)
discrete_resource_clone_cpp <- function(x, n = 1) {
  if (!is.environment(x) || is.null(x$.ptr)) {
    stop("Expected a resource wrapper environment with field `.ptr`.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1 || n < 1) {
    stop("n must be a single integer >= 1", call. = FALSE)
  }
  xptrs <- discrete_resource_clone_xptrs_cpp(x, as.integer(n))
  lapply(xptrs, .resource_discrete_from_xptr)
}