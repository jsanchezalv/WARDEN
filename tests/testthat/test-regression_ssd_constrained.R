sm <- function(res, f) rowMeans(sapply(Filter(Negate(is.null), res), \(s) unname(s[[f]])))

# Shared constrained SSD model components
cssd_common_pt <- add_item(death = max(0.0000001, rnorm(n = 1, mean = 12, sd = 3)))
cssd_unique_pt <- add_item(
  fl.sick = 1, q_default = util.sick,
  c_default = cost.sick + if (arm == "int") { cost.int } else { 0 },
  success_blocking_bed = FALSE, had_to_queue = 0,
  time_in_queue = NA, time_start_queue = NA
)
cssd_events <- add_tte(arm = c("noint", "int"), evts = c("sick", "sicker", "death"), input = {
  sick   <- 0
  sicker <- draw_tte(1, dist = "exp", coef1 = coef_noint,
                     beta_tx = ifelse(arm == "int", HR_int, 1),
                     seed = random_seed_sicker_i[i])
})
cssd_reactions <-
  add_reactevt(name_evt = "sick", input = {
    shared_accumulator <- shared_accumulator$modify(shared_accumulator$value() + 1)
    value_accum <- shared_accumulator$value()
    beds_free <- beds$n_free(); time_in_queue <- NA
  }) |>
  add_reactevt(name_evt = "sicker", input = {
    success_blocking_bed <- beds$attempt_block()
    beds_free <- beds$n_free()
    if (!success_blocking_bed) {
      time_start_queue <- curtime
      modify_event(c(death = max(curtime, get_event("death") * 0.8)))
      had_to_queue <- 1
    } else {
      time_in_queue <- ifelse(had_to_queue == 1, curtime - time_start_queue, NA)
    }
    q_default <- util.sicker
    c_default <- cost.sicker + if (arm == "int") { cost.int } else { 0 }
    fl.sick <- 0
  }) |>
  add_reactevt(name_evt = "death", input = {
    beds$attempt_free()
    if (success_blocking_bed & beds$queue_size() > 0) {
      new_event(c(sicker = curtime), cur_evtlist, patient_id = beds$next_patient_in_line())
    }
    success_blocking_bed <- FALSE; time_in_queue <- NA
    beds_free <- beds$n_free(); q_default <- 0; c_default <- 0; curtime <- Inf
  })

cssd_sens_inputs <- add_item(
  indicators = if (sensitivity_bool) {
    create_indicators(sens, n_sensitivity * length(sensitivity_names), rep(1, 7L))
  } else { rep(1, 7L) }
)
cssd_dsa_common_all <- add_item(
  pick_val_v(
    base = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
    psa  = pick_psa(
      list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
      rep(1, 7),
      list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)),
      lapply(list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)), \(x) abs(x / 5))
    ),
    sens = list(
      base_value = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
      DSA_min    = list(0.6, 0.3, 1000, 5000,  800, log(0.1), 0.5),
      DSA_max    = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9)
    )[[sens_name_used]],
    psa_ind = psa_bool, sens_ind = sensitivity_bool,
    indicator = indicators,
    names_out = list("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int")
  )
) |>
  add_item(input = {
    random_seed_sicker_i = sample(1:1000, 1000, replace = FALSE)
    beds <- resource_discrete(1000); beds_free <- beds$n_free()
    shared_accumulator <- shared_input(0); value_accum <- shared_accumulator$value()
  })

test_that("SSD constrained regression - 650 beds (npats=200, seed=42)", {
  skip_on_cran()
  cssd_common_all_c <- add_item(input = {
    util.sick <- 0.8; util.sicker <- 0.5; cost.sick <- 3000; cost.sicker <- 7000; cost.int <- 1000
    coef_noint <- log(0.2); HR_int <- 0.8; drc <- 0.035; drq <- 0.035
    random_seed_sicker_i <- sample.int(100000, npats, replace = FALSE)
    beds <- resource_discrete(650); beds_free <- beds$n_free()
    shared_accumulator <- shared_input(0); value_accum <- shared_accumulator$value()
  })
  suppressMessages({
    res <- run_sim(
      npats = 200, n_sim = 1, psa_bool = FALSE, arm_list = c("int", "noint"),
      common_all_inputs = cssd_common_all_c, common_pt_inputs = cssd_common_pt,
      unique_pt_inputs = cssd_unique_pt, init_event_list = cssd_events,
      evt_react_list = cssd_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", constrained = TRUE, ipd = 1,
      input_out = c("beds_free","had_to_queue","time_in_queue","value_accum"), seed = 42
    )
    out <- res[[1]][[1]]
    tc  <- unname(out$total_costs); tq <- unname(out$total_qalys)
    dc  <- tc[1] - tc[2]; dq <- tq[1] - tq[2]
  })
  expect_equal(tc, c(58296.59031654, 50837.15629680), tolerance = 1e-4)
  expect_equal(tq, c(6.26533273, 6.09950094),         tolerance = 1e-4)
  expect_equal(dc, 7459.43401974,  tolerance = 1e-4)
  expect_equal(dq, 0.16583179,     tolerance = 1e-4)
  expect_equal(dc / dq, 44981.93038115, tolerance = 1e-3)
})

test_that("SSD constrained regression - unconstrained (npats=200, seed=42)", {
  skip_on_cran()
  cssd_common_all_c <- add_item(input = {
    util.sick <- 0.8; util.sicker <- 0.5; cost.sick <- 3000; cost.sicker <- 7000; cost.int <- 1000
    coef_noint <- log(0.2); HR_int <- 0.8; drc <- 0.035; drq <- 0.035
    random_seed_sicker_i <- sample.int(100000, npats, replace = FALSE)
    beds <- resource_discrete(650); beds_free <- beds$n_free()
    shared_accumulator <- shared_input(0); value_accum <- shared_accumulator$value()
  })
  suppressMessages({
    res <- run_sim(
      npats = 200, n_sim = 1, psa_bool = FALSE, arm_list = c("int", "noint"),
      common_all_inputs = cssd_common_all_c, common_pt_inputs = cssd_common_pt,
      unique_pt_inputs = cssd_unique_pt, init_event_list = cssd_events,
      evt_react_list = cssd_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", constrained = FALSE, ipd = 1,
      input_out = c("beds_free","had_to_queue","time_in_queue"), seed = 42
    )
    tc <- unname(res[[1]][[1]]$total_costs)
    tq <- unname(res[[1]][[1]]$total_qalys)
  })
  expect_equal(tc, c(58296.59031654, 50837.15629680), tolerance = 1e-4)
  expect_equal(tq, c(6.26533273, 6.09950094),         tolerance = 1e-4)
})

test_that("SSD constrained regression - unbinding 1000 beds (npats=200, seed=42)", {
  skip_on_cran()
  cssd_common_all_ub <- add_item(input = {
    util.sick <- 0.8; util.sicker <- 0.5; cost.sick <- 3000; cost.sicker <- 7000; cost.int <- 1000
    coef_noint <- log(0.2); HR_int <- 0.8; drc <- 0.035; drq <- 0.035
    random_seed_sicker_i <- sample.int(100000, npats, replace = FALSE)
    beds <- resource_discrete(1000); beds_free <- beds$n_free()
    shared_accumulator <- shared_input(0); value_accum <- shared_accumulator$value()
  })
  suppressMessages({
    res <- run_sim(
      npats = 200, n_sim = 1, psa_bool = FALSE, arm_list = c("int", "noint"),
      common_all_inputs = cssd_common_all_ub, common_pt_inputs = cssd_common_pt,
      unique_pt_inputs = cssd_unique_pt, init_event_list = cssd_events,
      evt_react_list = cssd_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", constrained = TRUE, ipd = 1,
      input_out = c("beds_free","had_to_queue","time_in_queue"), seed = 42
    )
    tc <- unname(res[[1]][[1]]$total_costs)
    tq <- unname(res[[1]][[1]]$total_qalys)
  })
  expect_equal(tc, c(58296.59031654, 50837.15629680), tolerance = 1e-4)
  expect_equal(tq, c(6.26533273, 6.09950094),         tolerance = 1e-4)
})

test_that("SSD constrained regression - DSA (npats=50, n_sim=1, seed=42)", {
  skip_on_cran()
  suppressMessages({
    res <- run_sim(
      npats = 50, n_sim = 1, psa_bool = FALSE, arm_list = c("int", "noint"),
      common_all_inputs = cssd_dsa_common_all, common_pt_inputs = cssd_common_pt,
      unique_pt_inputs = cssd_unique_pt, init_event_list = cssd_events,
      evt_react_list = cssd_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", sensitivity_inputs = cssd_sens_inputs,
      sensitivity_names = c("DSA_min","DSA_max"), sensitivity_bool = TRUE,
      n_sensitivity = 7L, constrained = TRUE,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),
      seed = 42
    )
    n_p <- 7L
  })
  # sens=1: base case
  expect_equal(unname(res[[1]][[1]]$total_costs),
               c(58723.77298567, 50996.99066339), tolerance = 1e-4)
  expect_equal(unname(res[[1]][[1]]$total_qalys),
               c(5.59193154, 5.53295731), tolerance = 1e-4)
  # sens=2: DSA_min param 1 (util.sick min)
  expect_equal(unname(res[[2]][[1]]$total_qalys),
               c(5.77100497, 5.47613386), tolerance = 1e-4)
  # sens=n_p+1: DSA_min last param (HR_int min)
  expect_equal(unname(res[[n_p + 1]][[1]]$total_qalys),
               c(6.69004338, 6.51312071), tolerance = 1e-4)
})

test_that("SSD constrained regression - probabilistic DSA (npats=50, n_sim=3, seed=42)", {
  skip_on_cran()
  suppressMessages({
    res <- run_sim(
      npats = 50, n_sim = 3, psa_bool = TRUE, arm_list = c("int", "noint"),
      common_all_inputs = cssd_dsa_common_all, common_pt_inputs = cssd_common_pt,
      unique_pt_inputs = cssd_unique_pt, init_event_list = cssd_events,
      evt_react_list = cssd_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", sensitivity_inputs = cssd_sens_inputs,
      sensitivity_names = c("DSA_min","DSA_max"), sensitivity_bool = TRUE,
      n_sensitivity = 7L, constrained = TRUE,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int","beds_free"),
      seed = 42
    )
  })
  expect_equal(sm(res[[1]], "total_costs"),
               c(53082.43004964, 45053.77933053), tolerance = 1e-4)
  expect_equal(sm(res[[1]], "total_qalys"),
               c(5.40829162, 5.32587271), tolerance = 1e-4)
  expect_equal(sm(res[[2]], "total_qalys"),
               c(6.57128610, 6.24504612), tolerance = 1e-4)
})

test_that("SSD constrained regression - PSA 65 beds (npats=50, n_sim=5, seed=42)", {
  skip_on_cran()
  cssd_common_all_psa <- add_item(input = {
    util.sick <- 0.8; util.sicker <- 0.5; cost.sick <- 3000; cost.sicker <- 7000; cost.int <- 1000
    coef_noint <- log(0.2); HR_int <- 0.8; drc <- 0.035; drq <- 0.035
    random_seed_sicker_i <- sample.int(100000, npats, replace = FALSE)
    beds <- resource_discrete(65); beds_free <- beds$n_free()
    shared_accumulator <- shared_input(0); value_accum <- shared_accumulator$value()
  })
  suppressMessages({
    res <- run_sim(
      npats = 50, n_sim = 5, psa_bool = TRUE, arm_list = c("int", "noint"),
      common_all_inputs = cssd_common_all_psa, common_pt_inputs = cssd_common_pt,
      unique_pt_inputs = cssd_unique_pt, init_event_list = cssd_events,
      evt_react_list = cssd_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", sensitivity_inputs = cssd_sens_inputs,
      sensitivity_bool = FALSE, n_sensitivity = 1, constrained = TRUE,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int","beds_free"),
      seed = 42
    )
  })
  expect_equal(sm(res[[1]], "total_costs"),
               c(59465.42474006, 52098.67606918), tolerance = 1e-4)
  expect_equal(sm(res[[1]], "total_qalys"),
               c(6.30041193, 6.11925998), tolerance = 1e-4)
})
