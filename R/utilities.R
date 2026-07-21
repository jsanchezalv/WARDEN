# Utility functions for common HTA patterns --------------------------------

#' Compute effective hazard ratio with treatment waning
#'
#' Returns the effective hazard ratio at a given time on treatment,
#' accounting for waning of treatment effect.
#'
#' @param base_hr Numeric. Initial treatment HR (< 1 means benefit).
#' @param time_on_tx Numeric. Time the patient has been on treatment.
#' @param start Numeric. Time at which waning begins (HR = base_hr before
#'   this).
#' @param duration Numeric. Duration of the waning period (HR reaches 1 at
#'   start + duration).
#' @param type Character. Waning function: `"linear"` (default),
#'   `"exponential"`, or `"instant"`.
#'
#' @return Numeric scalar: effective HR at `time_on_tx`.
#' @export
waning_hr <- function(base_hr, time_on_tx, start, duration, type = "linear") {
  if (time_on_tx <= start) return(base_hr)
  if (time_on_tx >= start + duration || type == "instant") return(1.0)
  frac <- (time_on_tx - start) / duration
  switch(type,
    linear      = base_hr + frac * (1 - base_hr),
    exponential = base_hr * exp(frac * log(1 / base_hr)),
    stop('waning_hr(): unknown type "', type,
         '". Use "linear", "exponential", or "instant".',
         call. = FALSE)
  )
}

#' Draw a categorical index from probabilities
#'
#' Given a uniform random draw and a vector of probabilities, returns the
#' 1-based integer index of the selected category.
#'
#' @param luck Numeric between 0 and 1 (inclusive). A uniform random draw.
#' @param probs Numeric vector of probabilities. Will be normalized to sum
#'   to 1 internally. All values must be >= 0.
#'
#' @return Integer: 1-based index of the selected category.
#' @export
draw_categorical <- function(luck, probs) {
  if (any(probs < 0)) {
    stop("draw_categorical(): all probabilities must be >= 0.", call. = FALSE)
  }
  cdf <- cumsum(probs / sum(probs))
  which.max(luck <= cdf)
}

#' Quantile function for categorical distribution
#'
#' Alias for [draw_categorical()] using standard quantile-function naming.
#'
#' @param p Numeric between 0 and 1 (inclusive). Probability (same as `luck`
#'   in [draw_categorical()]).
#' @param probs Numeric vector of probabilities.
#'
#' @return Integer: 1-based index.
#' @export
qcategorical <- function(p, probs) {
  draw_categorical(luck = p, probs = probs)
}

#' Scale remaining time to an event by a hazard ratio
#'
#' Computes `curtime + (event_time - curtime) / hr`. Use inside event
#' reactions to shrink or expand the remaining time to an event.
#'
#' @param event_name Character. Name of the event whose time to scale.
#' @param hr Numeric greater than 0. Hazard ratio to apply. Values below 1
#'   extend time; values above 1 shorten it.
#'
#' @return Numeric: the new absolute event time.
#' @export
scale_remaining_time <- function(event_name, hr) {
  curtime     <- get("curtime",     envir = parent.frame(), inherits = TRUE)
  cur_evtlist <- get("cur_evtlist", envir = parent.frame(), inherits = TRUE)
  i           <- get("i",           envir = parent.frame(), inherits = TRUE)
  event_time  <- get_event_cpp(cur_evtlist, i, event_name)
  curtime + (event_time - curtime) / hr
}
