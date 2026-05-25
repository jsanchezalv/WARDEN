sm <- function(res, f) rowMeans(sapply(Filter(Negate(is.null), res), \(s) unname(s[[f]])))

# Shared inputs_selector model components — all data inlined in deferred expressions
is_sens_inputs <- add_item(iterator_sensitivity = sens_iterator(sens, n_sensitivity)) |>
  add_item(
    indicators = if (sensitivity_bool & sens_name_used %in% c("DSA_min","DSA_max")) {
      create_indicators(iterator_sensitivity, n_sensitivity * length(sensitivity_names),
                        rep(1, 7L))
    } else { rep(1, 7L) }
  )
is_simple <- add_item() |>
  add_item(
    pick_val_v(
      base = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
      psa  = pick_psa(
        list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
        as.list(rep(1, 7)),
        list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)),
        lapply(list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)), \(x) abs(x / 5))
      ),
      sens = list(
        base_value = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
        DSA_min    = list(0.6, 0.3, 1000, 5000,  800, log(0.1), 0.5),
        DSA_max    = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9),
        scenario_1 = list(0.6, 0.3, 1000, 5000,  800, log(0.1), 0.5),
        scenario_2 = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9)
      )[[sens_name_used]],
      psa_ind = psa_bool, sens_ind = sensitivity_bool,
      indicator = indicators,
      names_out = list("util.sick","util.sicker","cost.sick","cost.sicker",
                       "cost.int","coef_noint","HR_int"),
      indicator_psa = as.list(c(rep(1, 4), rep(0, 3)))
    )
  )
is_arm      <- add_item(q_default = util.sick, c_default = cost.sick + if (arm == "int") { cost.int } else { 0 })
is_events   <- add_tte(arm = c("noint","int"), evts = c("a1","b1"), input = { a1 <- 0; b1 <- 2 })
is_reactions <- add_reactevt(name_evt = "a1", input = {}) |>
  add_reactevt(name_evt = "b1", input = { q_default = 0; c_default = 0; curtime = Inf })

test_that("IS regression - deterministic (npats=5, seed=42)", {
  skip_on_cran()
  suppressMessages({
    res <- run_sim(
      npats = 5, n_sim = 1, psa_bool = FALSE, arm_list = c("int","noint"),
      common_all_inputs = is_simple, unique_pt_inputs = is_arm,
      init_event_list = is_events, evt_react_list = is_reactions,
      util_ongoing_list = "q_default", cost_ongoing_list = "c_default", ipd = 1,
      sensitivity_inputs = is_sens_inputs, sensitivity_names = NULL,
      sensitivity_bool = FALSE, n_sensitivity = 1, seed = 42
    )
    tc <- unname(res[[1]][[1]]$total_costs)
    tq <- unname(res[[1]][[1]]$total_qalys)
    dc <- tc[1] - tc[2]; dq <- tq[1] - tq[2]
  })
  expect_equal(tc, c(7768.12137341, 5826.09103006), tolerance = 1e-4)
  expect_equal(tq, c(1.55362427, 1.55362427),       tolerance = 1e-4)
  expect_equal(dc, 1942.03034335, tolerance = 1e-4)
  expect_equal(dq, 0,             tolerance = 1e-6)
  expect_true(is.infinite(dc / dq))
})

test_that("IS regression - PSA (npats=5, n_sim=2, seed=42)", {
  skip_on_cran()
  suppressMessages({
    res <- run_sim(
      npats = 5, n_sim = 2, psa_bool = TRUE, arm_list = c("int","noint"),
      common_all_inputs = is_simple, unique_pt_inputs = is_arm,
      init_event_list = is_events, evt_react_list = is_reactions,
      util_ongoing_list = "q_default", cost_ongoing_list = "c_default", ipd = 1,
      sensitivity_inputs = is_sens_inputs, sensitivity_names = NULL,
      sensitivity_bool = FALSE, n_sensitivity = 1, seed = 42
    )
  })
  expect_equal(sm(res[[1]], "total_costs"),
               c(6895.87622813, 4953.84588478), tolerance = 1e-4)
  expect_equal(sm(res[[1]], "total_qalys"),
               c(1.81975633, 1.81975633), tolerance = 1e-4)
})

test_that("IS regression - DSA (npats=5, n_sim=2, seed=42)", {
  skip_on_cran()
  suppressMessages({
    res <- run_sim(
      npats = 5, n_sim = 2, psa_bool = TRUE, arm_list = c("int","noint"),
      common_all_inputs = is_simple, unique_pt_inputs = is_arm,
      init_event_list = is_events, evt_react_list = is_reactions,
      util_ongoing_list = "q_default", cost_ongoing_list = "c_default", ipd = 1,
      sensitivity_inputs = is_sens_inputs,
      sensitivity_names = c("DSA_min","DSA_max"), sensitivity_bool = TRUE,
      n_sensitivity = 7L,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker",
                    "cost.int","coef_noint","HR_int"),
      seed = 42
    )
  })
  expect_equal(sm(res[[1]], "total_costs"),
               c(6895.87622813, 4953.84588478), tolerance = 1e-4)
  expect_equal(sm(res[[1]], "total_qalys"),
               c(1.16521821, 1.16521821), tolerance = 1e-4)
  expect_equal(sm(res[[2]], "total_qalys"),
               c(1.81975633, 1.81975633), tolerance = 1e-4)
  expect_equal(sm(res[[8]], "total_qalys"),
               c(1.74782731, 1.74782731), tolerance = 1e-4)
})

test_that("IS regression - scenario analysis (npats=5, n_sim=2, seed=42)", {
  skip_on_cran()
  suppressMessages({
    res <- run_sim(
      npats = 5, n_sim = 2, psa_bool = TRUE, arm_list = c("int","noint"),
      common_all_inputs = is_simple, unique_pt_inputs = is_arm,
      init_event_list = is_events, evt_react_list = is_reactions,
      util_ongoing_list = "q_default", cost_ongoing_list = "c_default", ipd = 1,
      sensitivity_inputs = is_sens_inputs,
      sensitivity_names = c("scenario_1","scenario_2"), sensitivity_bool = TRUE,
      n_sensitivity = 1,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker",
                    "cost.int","coef_noint","HR_int"),
      seed = 42
    )
  })
  expect_equal(sm(res[[1]], "total_costs"),
               c(3495.65461804, 1942.03034335), tolerance = 1e-4)
  expect_equal(sm(res[[1]], "total_qalys"),
               c(1.16521821, 1.16521821), tolerance = 1e-4)
  expect_equal(sm(res[[2]], "total_costs"),
               c(13594.21240347, 9710.15171676), tolerance = 1e-4)
  expect_equal(sm(res[[2]], "total_qalys"),
               c(1.74782731, 1.74782731), tolerance = 1e-4)
})

test_that("IS regression - split params (npats=5, n_sim=2, psa_bool=FALSE, seed=42)", {
  skip_on_cran()
  is_sens_split <- add_item(iterator_sensitivity = sens_iterator(sens, n_sensitivity)) |>
    add_item(
      indicators = if (sensitivity_bool & sens_name_used %in% c("DSA_min","DSA_max")) {
        create_indicators(iterator_sensitivity, n_sensitivity * length(sensitivity_names),
                          rep(1, 7L))
      } else { rep(1, 7L) }
    ) |>
    add_item(
      indicators_pat = if (sensitivity_bool & sens_name_used %in% c("DSA_min","DSA_max")) {
        create_indicators(iterator_sensitivity, n_sensitivity * length(sensitivity_names),
                          rep(1, 2L), 7L)
      } else { rep(1, 2L) }
    )
  is_pat <- add_item() |>
    add_item(pick_val_v(
      base = list(60, 1),
      psa  = pick_psa(list("rnorm","rbinom"), as.list(rep(1, 2)),
                      list(60, 1), list(10, 0.5)),
      sens = list(
        base_value = list(60, 1),
        DSA_min    = list(30, 0),
        DSA_max    = list(80, 1),
        scenario_1 = list(55, 1),
        scenario_2 = list(45, 0)
      )[[sens_name_used]],
      psa_ind = psa_bool, sens_ind = sensitivity_bool,
      indicator = indicators_pat, names_out = list("age","sex"),
      indicator_psa = as.list(rep(1, 2))
    ))
  suppressMessages({
    res <- run_sim(
      npats = 5, n_sim = 2, psa_bool = FALSE, arm_list = c("int","noint"),
      common_all_inputs = is_simple, unique_pt_inputs = is_arm,
      common_pt_inputs = is_pat, init_event_list = is_events,
      evt_react_list = is_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", ipd = 1,
      sensitivity_inputs = is_sens_split,
      sensitivity_names = c("DSA_min","DSA_max"), sensitivity_bool = TRUE,
      n_sensitivity = 9L,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker",
                    "cost.int","coef_noint","HR_int","age","sex"),
      seed = 42
    )
  })
  expect_equal(unname(res[[1]][[1]]$total_costs),
               c(7768.12137341, 5826.09103006), tolerance = 1e-4)
  expect_equal(unname(res[[1]][[1]]$total_qalys),
               c(1.16521821, 1.16521821), tolerance = 1e-4)
  expect_equal(unname(res[[2]][[1]]$total_qalys),
               c(1.55362427, 1.55362427), tolerance = 1e-4)
  expect_equal(unname(res[[8]][[1]]$total_qalys),
               c(1.55362427, 1.55362427), tolerance = 1e-4)
})

test_that("IS regression - covaried params (npats=5, n_sim=2, psa_bool=FALSE, seed=42)", {
  skip_on_cran()
  is_sens_cov <- add_item(iterator_sensitivity = sens_iterator(sens, n_sensitivity))
  is_simple_cov <- add_item() |>
    add_item(pick_val_v(
      base = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
      psa  = pick_psa(
        list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
        as.list(rep(1, 7)),
        list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)),
        lapply(list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)), \(x) abs(x / 5))
      ),
      sens = list(
        base_value = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
        DSA_min    = list(0.6, 0.3, 1000, 5000,  800, log(0.1), 0.5),
        DSA_max    = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9),
        scenario_1 = list(0.6, 0.3, 1000, 5000,  800, log(0.1), 0.5),
        scenario_2 = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9)
      )[[sens_name_used]],
      psa_ind = psa_bool, sens_ind = sensitivity_bool,
      indicator = as.list(c(1L, 1L, 2L, 2L, 2L, 3L, 4L)),
      sens_iterator = iterator_sensitivity, indicator_sens_binary = FALSE,
      names_out = list("util.sick","util.sicker","cost.sick","cost.sicker",
                       "cost.int","coef_noint","HR_int"),
      indicator_psa = as.list(c(rep(1, 4), rep(0, 3))),
      distributions = list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
      covariances = lapply(list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)), \(x) abs(x / 5))
    ))
  is_pat_cov <- add_item() |>
    add_item(pick_val_v(
      base = list(60, 1),
      psa  = pick_psa(list("rnorm","rbinom"), as.list(rep(1, 2)), list(60, 1), list(10, 0.5)),
      sens = list(
        base_value = list(60, 1),
        DSA_min    = list(30, 0),
        DSA_max    = list(80, 1),
        scenario_1 = list(55, 1),
        scenario_2 = list(45, 0)
      )[[sens_name_used]],
      psa_ind = psa_bool, sens_ind = sensitivity_bool,
      indicator = list(5L, 5L), sens_iterator = iterator_sensitivity,
      indicator_sens_binary = FALSE, names_out = list("age","sex"),
      indicator_psa = as.list(rep(1, 2)),
      distributions = list("rnorm","rbinom"), covariances = list(10, 0.5)
    ))
  suppressMessages({
    res <- run_sim(
      npats = 5, n_sim = 2, psa_bool = FALSE, arm_list = c("int","noint"),
      common_all_inputs = is_simple_cov, unique_pt_inputs = is_arm,
      common_pt_inputs = is_pat_cov, init_event_list = is_events,
      evt_react_list = is_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", ipd = 1,
      sensitivity_inputs = is_sens_cov,
      sensitivity_names = c("DSA_min","DSA_max"), sensitivity_bool = TRUE,
      n_sensitivity = 5L,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker",
                    "cost.int","coef_noint","HR_int","age","sex"),
      seed = 42
    )
  })
  expect_equal(unname(res[[1]][[1]]$total_costs),
               c(7768.12137341, 5826.09103006), tolerance = 1e-4)
  expect_equal(unname(res[[1]][[1]]$total_qalys),
               c(1.16521821, 1.16521821), tolerance = 1e-4)
  expect_equal(unname(res[[2]][[1]]$total_costs),
               c(3495.65461804, 1942.03034335), tolerance = 1e-4)
  expect_equal(unname(res[[3]][[1]]$total_qalys),
               c(1.55362427, 1.55362427), tolerance = 1e-4)
})

test_that("IS regression - vector params (npats=5, n_sim=2, psa_bool=FALSE, seed=42)", {
  skip_on_cran()
  is_sens_cov <- add_item(iterator_sensitivity = sens_iterator(sens, n_sensitivity))
  is_simple_cov <- add_item() |>
    add_item(pick_val_v(
      base = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
      psa  = pick_psa(
        list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
        as.list(rep(1, 7)),
        list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)),
        lapply(list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)), \(x) abs(x / 5))
      ),
      sens = list(
        base_value = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
        DSA_min    = list(0.6, 0.3, 1000, 5000,  800, log(0.1), 0.5),
        DSA_max    = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9),
        scenario_1 = list(0.6, 0.3, 1000, 5000,  800, log(0.1), 0.5),
        scenario_2 = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9)
      )[[sens_name_used]],
      psa_ind = psa_bool, sens_ind = sensitivity_bool,
      indicator = as.list(c(1L, 1L, 2L, 2L, 2L, 3L, 4L)),
      sens_iterator = iterator_sensitivity, indicator_sens_binary = FALSE,
      names_out = list("util.sick","util.sicker","cost.sick","cost.sicker",
                       "cost.int","coef_noint","HR_int"),
      indicator_psa = as.list(c(rep(1, 4), rep(0, 3))),
      distributions = list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
      covariances = lapply(list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)), \(x) abs(x / 5))
    ))
  is_pat_vec <- add_item() |>
    add_item(pick_val_v(
      base = list(60, 1, c(10, 20)),
      psa  = pick_psa(list("rnorm","rbinom","mvrnorm"), as.list(rep(1, 3)),
                      list(60, 1, c(10, 20)), list(10, 0.5, matrix(c(2, 1, 4, 1), 2, 2))),
      sens = list(
        base_value = list(60, 1, c(10, 20)),
        DSA_min    = list(30, 0, c(5, 10)),
        DSA_max    = list(80, 1, c(15, 25)),
        scenario_1 = list(55, 1, c(12, 21)),
        scenario_2 = list(45, 0, c(16, 10))
      )[[sens_name_used]],
      psa_ind = psa_bool, sens_ind = sensitivity_bool,
      indicator = list(5L, 5L, c(6L, 6L)), sens_iterator = iterator_sensitivity,
      indicator_sens_binary = FALSE, names_out = list("age","sex","v_state"),
      indicator_psa = list(1, 1, c(1, 0)),
      distributions = list("rnorm","rbinom","mvrnorm"),
      covariances = list(10, 0.5, matrix(c(2, 1, 4, 1), 2, 2))
    ))
  suppressMessages({
    res <- run_sim(
      npats = 5, n_sim = 2, psa_bool = FALSE, arm_list = c("int","noint"),
      common_all_inputs = is_simple_cov, unique_pt_inputs = is_arm,
      common_pt_inputs = is_pat_vec, init_event_list = is_events,
      evt_react_list = is_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", ipd = 1,
      sensitivity_inputs = is_sens_cov,
      sensitivity_names = c("DSA_min","DSA_max"), sensitivity_bool = TRUE,
      n_sensitivity = 6L,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker",
                    "cost.int","coef_noint","HR_int","age","sex"),
      seed = 42
    )
  })
  expect_equal(unname(res[[1]][[1]]$total_costs),
               c(7768.12137341, 5826.09103006), tolerance = 1e-4)
  expect_equal(unname(res[[1]][[1]]$total_qalys),
               c(1.16521821, 1.16521821), tolerance = 1e-4)
  expect_equal(unname(res[[2]][[1]]$total_costs),
               c(3495.65461804, 1942.03034335), tolerance = 1e-4)
  expect_equal(unname(res[[3]][[1]]$total_qalys),
               c(1.55362427, 1.55362427), tolerance = 1e-4)
})
