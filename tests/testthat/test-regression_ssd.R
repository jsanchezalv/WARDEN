sm <- function(res, f) rowMeans(sapply(Filter(Negate(is.null), res), \(s) unname(s[[f]])))

# Shared SSD model components — only references WARDEN simulation variables
ssd_common_all <- add_item(input = {
  util.sick   <- 0.8; util.sicker <- 0.5
  cost.sick   <- 3000; cost.sicker <- 7000; cost.int <- 1000
  coef_noint  <- log(0.2); coef_death <- log(0.05)
  HR_int <- 0.7; drc <- 0.035; drq <- 0.035
  random_seed_sicker_i <- sample.int(100000, npats, replace = FALSE)
  random_seed_death_i  <- runif(npats)
})
ssd_common_pt  <- add_item(death = qexp(random_seed_death_i[i], exp(coef_death)))
ssd_unique_pt  <- add_item(
  fl.sick = 1, q_default = util.sick,
  c_default = cost.sick + if (arm == "int") { cost.int } else { 0 }
)
ssd_events <- add_tte(arm = c("noint", "int"), evts = c("sick", "sicker", "death"), input = {
  sick   <- 0
  sicker <- draw_tte(1, dist = "exp", coef1 = coef_noint,
                     beta_tx = ifelse(arm == "int", HR_int, 1),
                     seed = random_seed_sicker_i[i])
})
ssd_reactions <-
  add_reactevt(name_evt = "sick",   input = {}) |>
  add_reactevt(name_evt = "sicker", input = {
    q_default <- util.sicker
    c_default <- cost.sicker + if (arm == "int") { cost.int } else { 0 }
    fl.sick <- 0
    modify_event(c("death" = max(curtime + (get_event("death") - curtime) * 0.8, curtime)))
  }) |>
  add_reactevt(name_evt = "death", input = { q_default <- 0; c_default <- 0; curtime <- Inf })

# DSA inputs — all data inlined so no file-level variables are referenced inside deferred exprs
ssd_sens_inputs <- add_item(
  indicators = if (sensitivity_bool) {
    create_indicators(sens, n_sensitivity * length(sensitivity_names), rep(1, 8L))
  } else { rep(1, 8L) }
)
ssd_dsa_common_all <- add_item(
  pick_val_v(
    base = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.05), 0.7),
    psa  = pick_psa(
      list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rnorm","rlnorm"),
      rep(1, 8),
      list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.05), log(0.7)),
      lapply(list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.05), log(0.7)), \(x) abs(x / 5))
    ),
    sens = list(
      base_value = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.05), 0.7),
      DSA_min    = list(0.6, 0.3, 1000, 5000,  800, log(0.1), log(0.1),  0.5),
      DSA_max    = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), log(0.03), 0.9)
    )[[sens_name_used]],
    psa_ind = psa_bool, sens_ind = sensitivity_bool,
    indicator = indicators,
    names_out = list("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","coef_death","HR_int")
  )
) |>
  add_item(
    random_seed_sicker_i = sample(1:1000, 1000, replace = FALSE),
    random_seed_death_i  = runif(npats)
  )

test_that("SSD regression - deterministic (npats=200, seed=42)", {
  skip_on_cran()
  suppressMessages({
    res <- run_sim(
      npats = 200, n_sim = 1, psa_bool = FALSE, arm_list = c("int", "noint"),
      common_all_inputs = ssd_common_all, common_pt_inputs = ssd_common_pt,
      unique_pt_inputs = ssd_unique_pt, init_event_list = ssd_events,
      evt_react_list = ssd_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", ipd = 1, seed = 42
    )
    out <- res[[1]][[1]]
    tc  <- unname(out$total_costs); tq <- unname(out$total_qalys)
    dc  <- tc[1] - tc[2]; dq <- tq[1] - tq[2]
  })
  expect_equal(tc, c(69095.87045811, 61054.45126504), tolerance = 1e-4)
  expect_equal(tq, c(6.80157013, 6.47829882),         tolerance = 1e-4)
  expect_equal(dc, 8041.41919307,  tolerance = 1e-4)
  expect_equal(dq, 0.32327131,     tolerance = 1e-4)
  expect_equal(dc / dq, 24875.15,  tolerance = 1e-3)
})

test_that("SSD regression - DSA (npats=50, n_sim=1, seed=42)", {
  skip_on_cran()
  suppressMessages({
    res <- run_sim(
      npats = 50, n_sim = 1, psa_bool = FALSE, arm_list = c("int", "noint"),
      common_all_inputs = ssd_dsa_common_all, common_pt_inputs = ssd_common_pt,
      unique_pt_inputs = ssd_unique_pt, init_event_list = ssd_events,
      evt_react_list = ssd_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", sensitivity_inputs = ssd_sens_inputs,
      sensitivity_names = c("DSA_min", "DSA_max"), sensitivity_bool = TRUE,
      n_sensitivity = 8L,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker",
                    "cost.int","coef_noint","coef_death","HR_int"),
      seed = 42
    )
  })
  expect_equal(unname(res[[1]][[1]]$total_costs),
               c(66966.47293112, 59165.70798435), tolerance = 1e-4)
  expect_equal(unname(res[[1]][[1]]$total_qalys),
               c(5.90984544, 5.74230143), tolerance = 1e-4)
  # sens=2: DSA_min param 1 (util.sick min)
  expect_equal(unname(res[[2]][[1]]$total_qalys),
               c(5.71377581, 5.17494850), tolerance = 1e-4)
  # sens=9: DSA_min last param (HR_int min)
  expect_equal(unname(res[[9]][[1]]$total_qalys),
               c(6.89524023, 6.52846854), tolerance = 1e-4)
})

test_that("SSD regression - probabilistic DSA (npats=50, n_sim=3, seed=42)", {
  skip_on_cran()
  suppressMessages({
    res <- run_sim(
      npats = 50, n_sim = 3, psa_bool = TRUE, arm_list = c("int", "noint"),
      common_all_inputs = ssd_dsa_common_all, common_pt_inputs = ssd_common_pt,
      unique_pt_inputs = ssd_unique_pt, init_event_list = ssd_events,
      evt_react_list = ssd_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", sensitivity_inputs = ssd_sens_inputs,
      sensitivity_names = c("DSA_min", "DSA_max"), sensitivity_bool = TRUE,
      n_sensitivity = 8L,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker",
                    "cost.int","coef_noint","coef_death","HR_int"),
      seed = 42
    )
  })
  expect_equal(sm(res[[1]], "total_costs"),
               c(48425.92455638, 41798.88917101), tolerance = 1e-4)
  expect_equal(sm(res[[1]], "total_qalys"),
               c(4.90225962, 4.68184622), tolerance = 1e-4)
  expect_equal(sm(res[[2]], "total_qalys"),
               c(6.00973835, 5.35265384), tolerance = 1e-4)
})

test_that("SSD regression - simple PSA (npats=50, n_sim=5, seed=42)", {
  skip_on_cran()
  suppressMessages({
    res <- run_sim(
      npats = 50, n_sim = 5, psa_bool = TRUE, arm_list = c("int", "noint"),
      common_all_inputs = ssd_dsa_common_all, common_pt_inputs = ssd_common_pt,
      unique_pt_inputs = ssd_unique_pt, init_event_list = ssd_events,
      evt_react_list = ssd_reactions, util_ongoing_list = "q_default",
      cost_ongoing_list = "c_default", sensitivity_inputs = ssd_sens_inputs,
      sensitivity_bool = FALSE, n_sensitivity = 1L,
      input_out = c("util.sick","util.sicker","cost.sick","cost.sicker",
                    "cost.int","coef_noint","coef_death","HR_int"),
      seed = 42
    )
  })
  expect_equal(sm(res[[1]], "total_costs"),
               c(50772.00292763, 44617.76712723), tolerance = 1e-4)
  expect_equal(sm(res[[1]], "total_qalys"),
               c(5.85673061, 5.43802396), tolerance = 1e-4)
})
