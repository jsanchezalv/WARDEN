# Model validation --------------------------------------------------------

#' Validate model inputs before running
#'
#' Performs static checks on the model inputs without running the simulation.
#' Catches common configuration errors early.
#'
#' @param arm_list Character vector of intervention names.
#' @param init_event_list List of initial events and event times (output of
#'   [add_tte()]).
#' @param evt_react_list List of event reactions (output of [add_reactevt()]).
#' @param util_ongoing_list Character vector of ongoing QALY variable names.
#' @param util_instant_list Character vector of instant QALY variable names.
#' @param util_cycle_list Character vector of cycle QALY variable names.
#' @param cost_ongoing_list Character vector of ongoing cost variable names.
#' @param cost_instant_list Character vector of instant cost variable names.
#' @param cost_cycle_list Character vector of cycle cost variable names.
#' @param other_ongoing_list Character vector of other ongoing variable names.
#' @param other_instant_list Character vector of other instant variable names.
#' @param npats Integer. Number of patients per simulation.
#' @param n_sim Integer. Number of simulations.
#' @param psa_bool Logical. Whether PSA mode is intended.
#' @param seed Numeric. Random seed.
#'
#' @return Invisibly `TRUE` if all checks pass. Prints a summary of checks.
#' @export
validate_model <- function(arm_list = c("int", "noint"),
                           init_event_list = NULL,
                           evt_react_list = NULL,
                           util_ongoing_list = NULL,
                           util_instant_list = NULL,
                           util_cycle_list = NULL,
                           cost_ongoing_list = NULL,
                           cost_instant_list = NULL,
                           cost_cycle_list = NULL,
                           other_ongoing_list = NULL,
                           other_instant_list = NULL,
                           npats = 500,
                           n_sim = 1,
                           psa_bool = NULL,
                           seed = NULL) {

  passed <- 0L
  warned <- 0L

  check_pass <- function(msg) {
    passed <<- passed + 1L
    cat("  [PASS]", msg, "\n")
  }

  check_fail <- function(msg) {
    stop("validate_model(): ", msg, call. = FALSE)
  }

  check_warn <- function(msg) {
    warned <<- warned + 1L
    cat("  [WARN]", msg, "\n")
  }

  cat("Validating model inputs...\n")

 # 1. arm_list
  if (!is.character(arm_list) || length(arm_list) < 1L) {
    check_fail("arm_list must be a non-empty character vector.")
  }
  if (anyDuplicated(arm_list) > 0L) {
    check_fail(paste0("arm_list contains duplicates: ",
                      paste(arm_list[duplicated(arm_list)], collapse = ", ")))
  }
  check_pass(paste0("arm_list valid (", length(arm_list), " arms)"))

  # 2. init_event_list coverage
  if (!is.null(init_event_list)) {
    defined_arms <- names(init_event_list)
    missing_arms <- setdiff(arm_list, defined_arms)
    if (length(missing_arms) > 0L) {
      check_fail(paste0("init_event_list is missing arm(s): ",
                        paste(missing_arms, collapse = ", ")))
    }
    check_pass("All arms have entries in init_event_list")

    # 3. Event/reaction cross-check
    if (!is.null(evt_react_list)) {
      declared_evts <- init_event_list[[1]]$evts
      reaction_evts <- names(evt_react_list)

      missing_reactions <- setdiff(declared_evts, reaction_evts)
      if (length(missing_reactions) > 0L) {
        check_fail(paste0("Events declared in init_event_list but missing from ",
                          "evt_react_list: ", paste(missing_reactions, collapse = ", ")))
      }
      check_pass("All declared events have matching reactions")

      orphan_reactions <- setdiff(reaction_evts, declared_evts)
      if (length(orphan_reactions) > 0L) {
        check_warn(paste0("Reactions defined for undeclared events: ",
                          paste(orphan_reactions, collapse = ", ")))
      }

      if (anyDuplicated(reaction_evts) > 0L) {
        check_fail(paste0("Duplicate reactions: ",
                          paste(reaction_evts[duplicated(reaction_evts)], collapse = ", ")))
      }
      check_pass("No duplicate reactions")
    }
  } else {
    check_warn("init_event_list is NULL (a 'start' event at t=0 will be auto-created)")
  }

  # 4. Seed overflow
  if (!is.null(seed)) {
    max_product <- as.double(n_sim) * 1007 * as.double(seed) +
      as.double(npats) * 349 * as.double(seed)
    if (max_product > .Machine$integer.max) {
      safe_max <- floor(.Machine$integer.max /
                          (as.double(n_sim) * 1007 + as.double(npats) * 349))
      check_fail(paste0("seed = ", seed, " will overflow. Max safe value: ", safe_max))
    }
    check_pass("Seed within safe range")
  }

  # 5. Dimension sanity
  if (!is.numeric(npats) || npats < 1L) {
    check_fail("npats must be >= 1")
  }
  if (!is.numeric(n_sim) || n_sim < 1L) {
    check_fail("n_sim must be >= 1")
  }
  check_pass(paste0("Dimensions valid (npats=", npats, ", n_sim=", n_sim, ")"))

  # 6. PSA warning
  if (n_sim > 1L && identical(psa_bool, FALSE)) {
    check_warn("n_sim > 1 but psa_bool = FALSE. Did you mean to enable PSA?")
  }

  cat("\nValidation complete: ", passed, " passed",
      if (warned > 0L) paste0(", ", warned, " warning(s)") else "", ".\n", sep = "")
  invisible(TRUE)
}
