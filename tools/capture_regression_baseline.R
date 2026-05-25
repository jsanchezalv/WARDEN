# tools/capture_regression_baseline.R
# Run once to capture baseline regression values, then hardcode them in the test files.
# Usage: "/c/Program Files/R/R-4.2.3/bin/Rscript.exe" -e "source('tools/capture_regression_baseline.R')"

devtools::load_all(quiet = TRUE)
library(purrr)

fv <- function(x, d = 8) {
  paste0("c(", paste(formatC(x, digits = d, format = "f"), collapse = ", "), ")")
}
fs <- function(x, d = 8) formatC(x, digits = d, format = "f")

# helper to get sim means
sim_mean <- function(res_list, field) {
  vals <- sapply(Filter(Negate(is.null), res_list), function(s) unname(s[[field]]))
  rowMeans(vals)
}

# ==============================================================
# SHARED SSD BASE MODEL COMPONENTS (reused across det/DSA/PSA)
# ==============================================================
ssd_common_all <- add_item(input = {
  util.sick   <- 0.8
  util.sicker <- 0.5
  cost.sick   <- 3000
  cost.sicker <- 7000
  cost.int    <- 1000
  coef_noint  <- log(0.2)
  coef_death  <- log(0.05)
  HR_int      <- 0.7
  drc         <- 0.035
  drq         <- 0.035
  random_seed_sicker_i <- sample.int(100000, npats, replace = FALSE)
  random_seed_death_i  <- runif(npats)
})

ssd_common_pt <- add_item(death = qexp(random_seed_death_i[i], exp(coef_death)))

ssd_unique_pt <- add_item(
  fl.sick   = 1,
  q_default = util.sick,
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
    fl.sick   <- 0
    modify_event(c("death" = max(curtime + (get_event("death") - curtime) * 0.8, curtime)))
  }) |>
  add_reactevt(name_evt = "death", input = {
    q_default <- 0
    c_default <- 0
    curtime   <- Inf
  })

util_ongoing <- "q_default"
cost_ongoing <- "c_default"

# ==============================================================
# 1. SSD DETERMINISTIC
# ==============================================================
cat("\n=== 1. SSD DETERMINISTIC (npats=200, seed=42) ===\n")

res <- run_sim(
  npats = 200, n_sim = 1, psa_bool = FALSE,
  arm_list = c("int", "noint"),
  common_all_inputs = ssd_common_all,
  common_pt_inputs  = ssd_common_pt,
  unique_pt_inputs  = ssd_unique_pt,
  init_event_list   = ssd_events,
  evt_react_list    = ssd_reactions,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  ipd = 1, seed = 42
)
out <- res[[1]][[1]]
det <- summary_results_det(out)
cat(sprintf("total_costs:         %s\n", fv(unname(out$total_costs))))
cat(sprintf("total_costs_undisc:  %s\n", fv(unname(out$total_costs_undisc))))
cat(sprintf("total_qalys:         %s\n", fv(unname(out$total_qalys))))
cat(sprintf("total_qalys_undisc:  %s\n", fv(unname(out$total_qalys_undisc))))
cat(sprintf("dcosts[2]:           %s\n", fs(det$dcosts[2])))
cat(sprintf("dqalys[2]:           %s\n", fs(det$dqalys[2])))
cat(sprintf("ICUR[2]:             %s\n", fs(det$ICUR[2])))

# ==============================================================
# 2. SSD DSA
# ==============================================================
cat("\n=== 2. SSD DSA (npats=50, n_sim=1, seed=42) ===\n")

list_par <- list(
  parameter_name = list("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","coef_death","HR_int"),
  base_value     = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.05), 0.7),
  DSA_min        = list(0.6, 0.3, 1000, 5000,  800, log(0.1), log(0.1),  0.5),
  DSA_max        = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), log(0.03), 0.9),
  PSA_dist       = list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rnorm","rlnorm"),
  a              = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.05), log(0.7)),
  b              = lapply(list(0.8,0.5,3000,7000,1000,log(0.2),log(0.05),log(0.7)), \(x) abs(x/5)),
  scenario_1     = list(0.6, 0.3, 1000, 5000,  800, log(0.1), log(0.03), 0.5),
  scenario_2     = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), log(0.2),  0.9)
)

ssd_sens_inputs <- add_item(
  indicators = if (sensitivity_bool) {
    create_indicators(sens, n_sensitivity * length(sensitivity_names), rep(1, length(list_par[[1]])))
  } else {
    rep(1, length(list_par[[1]]))
  }
)

ssd_dsa_common_all <- add_item(
  pick_val_v(
    base     = list_par[["base_value"]],
    psa      = pick_psa(list_par[["PSA_dist"]], rep(1, length(list_par[["PSA_dist"]])), list_par[["a"]], list_par[["b"]]),
    sens     = list_par[[sens_name_used]],
    psa_ind  = psa_bool,
    sens_ind = sensitivity_bool,
    indicator  = indicators,
    names_out  = list_par[["parameter_name"]]
  )
) |>
  add_item(
    random_seed_sicker_i = sample(1:1000, 1000, replace = FALSE),
    random_seed_death_i  <- runif(npats)
  )

res_dsa <- run_sim(
  npats = 50, n_sim = 1, psa_bool = FALSE,
  arm_list = c("int", "noint"),
  common_all_inputs = ssd_dsa_common_all,
  common_pt_inputs  = ssd_common_pt,
  unique_pt_inputs  = ssd_unique_pt,
  init_event_list   = ssd_events,
  evt_react_list    = ssd_reactions,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  sensitivity_inputs = ssd_sens_inputs,
  sensitivity_names  = c("DSA_min", "DSA_max"),
  sensitivity_bool   = TRUE,
  n_sensitivity      = length(list_par[[1]]),
  input_out          = unlist(list_par[["parameter_name"]]),
  seed = 42
)

# base case (sens index 1) and two specific DSA iterations
for (si in c(1, 2, length(list_par[[1]]) + 1)) {
  o <- res_dsa[[si]][[1]]
  cat(sprintf("sens=%d total_costs: %s  total_qalys: %s\n",
              si, fv(unname(o$total_costs)), fv(unname(o$total_qalys))))
}

# ==============================================================
# 3. SSD PROBABILISTIC DSA
# ==============================================================
cat("\n=== 3. SSD PROB DSA (npats=50, n_sim=3, seed=42) ===\n")

res_pdsa <- run_sim(
  npats = 50, n_sim = 3, psa_bool = TRUE,
  arm_list = c("int", "noint"),
  common_all_inputs = ssd_dsa_common_all,
  common_pt_inputs  = ssd_common_pt,
  unique_pt_inputs  = ssd_unique_pt,
  init_event_list   = ssd_events,
  evt_react_list    = ssd_reactions,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  sensitivity_inputs = ssd_sens_inputs,
  sensitivity_names  = c("DSA_min", "DSA_max"),
  sensitivity_bool   = TRUE,
  n_sensitivity      = length(list_par[[1]]),
  input_out          = unlist(list_par[["parameter_name"]]),
  seed = 42
)

for (si in c(1, 2)) {
  cat(sprintf("sens=%d mean_costs: %s  mean_qalys: %s\n",
              si,
              fv(sim_mean(res_pdsa[[si]], "total_costs")),
              fv(sim_mean(res_pdsa[[si]], "total_qalys"))))
}

# ==============================================================
# 4. SSD SIMPLE PSA
# ==============================================================
cat("\n=== 4. SSD SIMPLE PSA (npats=50, n_sim=5, seed=42) ===\n")

res_psa <- run_sim(
  npats = 50, n_sim = 5, psa_bool = TRUE,
  arm_list = c("int", "noint"),
  common_all_inputs = ssd_dsa_common_all,
  common_pt_inputs  = ssd_common_pt,
  unique_pt_inputs  = ssd_unique_pt,
  init_event_list   = ssd_events,
  evt_react_list    = ssd_reactions,
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  sensitivity_inputs = ssd_sens_inputs,
  sensitivity_bool   = FALSE,
  n_sensitivity      = 1,
  input_out          = unlist(list_par[["parameter_name"]]),
  seed = 42
)

cat(sprintf("mean_costs:  %s\n", fv(sim_mean(res_psa[[1]], "total_costs"))))
cat(sprintf("mean_qalys:  %s\n", fv(sim_mean(res_psa[[1]], "total_qalys"))))

# ==============================================================
# 5. eBC DETERMINISTIC
# ==============================================================
cat("\n=== 5. eBC DETERMINISTIC (npats=200, seed=42) ===\n")

df_util <- data.frame(
  name  = c("util.idfs.ontx","util.idfs.offtx","util.remission","util.recurrence","util.mbc.progression.mbc","util.mbc.pps"),
  value = c(0.75, 0.8, 0.9, 0.7, 0.6, 0.5),
  se    = rep(0.02, 6)
)
df_cost <- data.frame(
  name  = c("cost.idfs.tx","cost.recurrence","cost.mbc.tx","cost.tx.beva","cost.idfs.txnoint","cost.idfs","cost.mbc.progression.mbc","cost.mbc.pps","cost.2ndline","cost.ae"),
  value = c(40000, 5000, 3000, 10000, 30000, 10000, 20000, 30000, 20000, 1000)
)
df_cost$se <- df_cost$value / 5

ebc_common_all <- add_item(input = {
  pick_val_v(base = df_util$value, psa = MASS::mvrnorm(1, df_util$value, diag(df_util$se^2)),
             sens = df_util$value, psa_ind = psa_bool, sens_ind = sensitivity_bool,
             indicator = rep(0, nrow(df_util)), names_out = df_util$name, deploy_env = TRUE)
  pick_val_v(base = df_cost$value, psa = rgamma_mse(1, df_cost$value, df_cost$se),
             sens = df_cost$value, psa_ind = psa_bool, sens_ind = sensitivity_bool,
             indicator = rep(0, nrow(df_cost)), names_out = df_cost$name, deploy_env = TRUE)
})

ebc_common_pt <- add_item(input = {
  sex_pt   <- ifelse(rbinom(1, 1, p = 0.01), "male", "female")
  nat.os.s <- rcond_gompertz(1,
                             shape = if (sex_pt == "male") { 0.102 } else { 0.115 },
                             rate  = if (sex_pt == "male") { 0.000016 } else { 0.0000041 },
                             lower_bound = 50)
  fl.remission <- rbinom(1, 1, 0.8)
})

ebc_unique_pt <- add_item(input = {
  fl.idfs.ontx            <- 1; fl.idfs <- 1; fl.mbcs.ontx <- 1
  fl.mbcs.progression.mbc <- 1; fl.tx.beva <- 1; fl.mbcs <- 0
  fl.mbcs_2ndline <- 0; fl.recurrence <- 0
  q_default <- if (fl.idfs == 1) {
    util.idfs.ontx * fl.idfs.ontx + (1 - fl.idfs.ontx) * (1 - fl.idfs.ontx)
  } else if (fl.idfs == 0 & fl.mbcs == 0) {
    util.remission * fl.remission + fl.recurrence * util.recurrence
  } else if (fl.mbcs == 1) {
    util.mbc.progression.mbc * fl.mbcs.progression.mbc + (1 - fl.mbcs.progression.mbc) * util.mbc.pps
  }
  c_default <- if (arm == "noint") {
    cost.idfs.txnoint * fl.idfs.ontx + cost.idfs
  } else {
    cost.idfs.tx * fl.idfs.ontx + cost.tx.beva * fl.tx.beva + cost.idfs
  }
  c_ae <- 0
  rnd_stream_ae  <- random_stream(100)
  rnd_stream_mbc <- random_stream(100)
})

ebc_events <-
  add_tte(arm = "int",
          evts = c("start","ttot","ttot.beva","progression.mbc","os","idfs","ttot.early","remission","recurrence","start.early.mbc","ae","2ndline_mbc"),
          other_inp = c("os.early","os.mbc"),
          input = {
            start <- 0
            idfs  <- draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2))
            ttot.early <- min(draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2)), idfs)
            ttot.beva  <- draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2))
            os.early   <- draw_tte(1, "lnorm", coef1 = 3, coef2 = log(0.2))
            if (fl.remission) {
              recurrence <- idfs + draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2))
              remission  <- idfs
              if (min(os.early, nat.os.s) > recurrence) {
                os.mbc        <- draw_tte(1, "lnorm", coef1 = 0.8, coef2 = log(0.2)) + idfs + recurrence
                progression.mbc <- draw_tte(1, "lnorm", coef1 = 0.5, coef2 = log(0.2)) + idfs + recurrence
                ttot <- draw_tte(1, "lnorm", coef1 = 0.5, coef2 = log(0.2)) + idfs + recurrence
              }
            } else {
              start.early.mbc <- draw_tte(1, "lnorm", coef1 = 2.3, coef2 = log(0.2))
              idfs <- ifelse(start.early.mbc < idfs, start.early.mbc, idfs)
              ttot.early <- min(ifelse(start.early.mbc < idfs, start.early.mbc, idfs), ttot.early)
              os.mbc <- draw_tte(1, "lnorm", coef1 = 0.8, coef2 = log(0.2)) + start.early.mbc
              progression.mbc <- draw_tte(1, "lnorm", coef1 = 0.5, coef2 = log(0.2)) + start.early.mbc
              ttot <- draw_tte(1, "lnorm", coef1 = 0.5, coef2 = log(0.2)) + start.early.mbc
            }
            os <- min(os.mbc, os.early, nat.os.s)
          }) |>
  add_tte(arm = "noint",
          evts = c("start","ttot","ttot.beva","progression.mbc","os","idfs","ttot.early","remission","recurrence","start.early.mbc"),
          other_inp = c("os.early","os.mbc"),
          input = {
            start <- 0
            idfs  <- draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2), beta_tx = 1.2)
            ttot.early <- min(draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2), beta_tx = 1.2), idfs)
            os.early   <- draw_tte(1, "lnorm", coef1 = 3, coef2 = log(0.2), beta_tx = 1.2)
            if (fl.remission) {
              recurrence <- idfs + draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2))
              remission  <- idfs
              if (min(os.early, nat.os.s) > recurrence) {
                os.mbc <- draw_tte(1, "lnorm", coef1 = 0.8, coef2 = log(0.2)) + idfs + recurrence
                progression.mbc <- draw_tte(1, "lnorm", coef1 = 0.5, coef2 = log(0.2)) + idfs + recurrence
                ttot <- draw_tte(1, "lnorm", coef1 = 0.5, coef2 = log(0.2)) + idfs + recurrence
              }
            } else {
              start.early.mbc <- draw_tte(1, "lnorm", coef1 = 2.3, coef2 = log(0.2))
              idfs <- ifelse(start.early.mbc < idfs, start.early.mbc, idfs)
              ttot.early <- min(ifelse(start.early.mbc < idfs, start.early.mbc, idfs), ttot.early)
              os.mbc <- draw_tte(1, "lnorm", coef1 = 0.8, coef2 = log(0.2)) + start.early.mbc
              progression.mbc <- draw_tte(1, "lnorm", coef1 = 0.5, coef2 = log(0.2)) + start.early.mbc
              ttot <- draw_tte(1, "lnorm", coef1 = 0.5, coef2 = log(0.2)) + start.early.mbc
            }
            os <- min(os.mbc, os.early, nat.os.s)
          })

.q_expr <- quote(
  if (fl.idfs==1) {
    util.idfs.ontx * fl.idfs.ontx + (1-fl.idfs.ontx)*(1-fl.idfs.ontx)
  } else if (fl.idfs==0 & fl.mbcs==0) {
    util.remission * fl.remission + fl.recurrence*util.recurrence
  } else if (fl.mbcs==1) {
    util.mbc.progression.mbc * fl.mbcs.progression.mbc + (1-fl.mbcs.progression.mbc)*util.mbc.pps
  }
)
.c_mbc_expr <- quote(
  cost.mbc.tx * fl.mbcs.ontx + cost.mbc.progression.mbc * fl.mbcs.progression.mbc +
    cost.mbc.pps * (1-fl.mbcs.progression.mbc) + cost.2ndline*fl.mbcs_2ndline
)

ebc_reactions <-
  add_reactevt("start",        input = {}) |>
  add_reactevt("ttot",         input = { q_default <- eval(.q_expr); c_default <- eval(.c_mbc_expr); fl.mbcs.ontx <- 0 }) |>
  add_reactevt("ttot.beva",    input = { q_default <- eval(.q_expr); c_default <- eval(.c_mbc_expr); fl.tx.beva <- 0 }) |>
  add_reactevt("progression.mbc", input = {
    q_default <- eval(.q_expr); c_default <- eval(.c_mbc_expr)
    fl.mbcs.progression.mbc <- 0; fl.mbcs_2ndline <- 1
    new_event(c("2ndline_mbc" = curtime + qexp(rnd_stream_mbc$draw_n(), 0.08)/12))
  }) |>
  add_reactevt("idfs",         input = {
    q_default = eval(.q_expr)
    c_default <- if (arm=="noint") { cost.idfs.txnoint*fl.idfs.ontx + cost.idfs } else { cost.idfs.tx*fl.idfs.ontx + cost.tx.beva*fl.tx.beva + cost.idfs }
    fl.idfs <- 0
  }) |>
  add_reactevt("ttot.early",   input = {
    q_default <- eval(.q_expr)
    c_default <- if (arm=="noint") { cost.idfs.txnoint*fl.idfs.ontx + cost.idfs } else { cost.idfs.tx*fl.idfs.ontx + cost.tx.beva*fl.tx.beva + cost.idfs }
    fl.idfs.ontx <- 0; fl.tx.beva <- 0
    n_ae <- qpois(rnd_stream_ae$draw_n(), lambda = 0.25*(curtime - prevtime))
    if (n_ae > 0) { new_event(c("ae" = curtime + 0.0001)); n_ae <- n_ae - 1 }
  }) |>
  add_reactevt("remission",    input = { q_default <- eval(.q_expr); c_default <- cost.recurrence * fl.recurrence; fl.remission <- 1 }) |>
  add_reactevt("recurrence",   input = { q_default <- eval(.q_expr); c_default <- cost.recurrence * fl.recurrence; fl.recurrence <- 1; fl.remission <- 0; fl.mbcs <- 1; fl.mbcs.progression.mbc <- 1 }) |>
  add_reactevt("start.early.mbc", input = { q_default <- eval(.q_expr); c_default <- cost.recurrence * fl.recurrence; fl.mbcs <- 1; fl.mbcs.progression.mbc <- 1 }) |>
  add_reactevt("2ndline_mbc",  input = {
    q_default <- eval(.q_expr); c_default <- eval(.c_mbc_expr); fl.mbcs_2ndline <- 0
    n_ae <- qpois(rnd_stream_ae$draw_n(), lambda = 0.25*(curtime - prevtime))
    if (n_ae > 0) { new_event(c("ae" = curtime + 0.0001)); n_ae <- n_ae - 1 }
  }) |>
  add_reactevt("ae", input = {
    if (n_ae > 0) { new_event(c("ae" = curtime)); n_ae <- n_ae - 1 }
    q_default = eval(.q_expr); c_default <- eval(.c_mbc_expr); c_ae <- cost.ae
    modify_event(c("os" = max(get_event("os") - 0.125, curtime + 0.0001)))
  }) |>
  add_reactevt("os", input = {
    q_default <- eval(.q_expr); c_default <- eval(.c_mbc_expr)
    fl.tx.beva <- 0; fl.mbcs.ontx <- 0; fl.idfs <- 0; fl.mbcs <- 0; curtime <- Inf
  })

res_ebc <- run_sim(
  npats = 200, n_sim = 1, psa_bool = FALSE,
  arm_list = c("int", "noint"),
  common_all_inputs = ebc_common_all,
  common_pt_inputs  = ebc_common_pt,
  unique_pt_inputs  = ebc_unique_pt,
  init_event_list   = ebc_events,
  evt_react_list    = ebc_reactions,
  util_ongoing_list = "q_default",
  cost_ongoing_list = "c_default",
  cost_instant_list = "c_ae",
  input_out = c("os.early","os.mbc","nat.os.s","sex_pt"),
  seed = 42
)

out_ebc <- res_ebc[[1]][[1]]
det_ebc <- summary_results_det(out_ebc)
cat(sprintf("total_costs:         %s\n", fv(unname(out_ebc$total_costs))))
cat(sprintf("total_costs_undisc:  %s\n", fv(unname(out_ebc$total_costs_undisc))))
cat(sprintf("total_qalys:         %s\n", fv(unname(out_ebc$total_qalys))))
cat(sprintf("total_qalys_undisc:  %s\n", fv(unname(out_ebc$total_qalys_undisc))))
cat(sprintf("dcosts[2]:           %s\n", fs(det_ebc$dcosts[2])))
cat(sprintf("dqalys[2]:           %s\n", fs(det_ebc$dqalys[2])))
cat(sprintf("ICUR[2]:             %s\n", fs(det_ebc$ICUR[2])))

# ==============================================================
# 6. SSD CONSTRAINED
# ==============================================================
cat("\n=== 6. SSD CONSTRAINED (npats=200, seed=42) ===\n")

cssd_common_all <- add_item(input = {
  util.sick   <- 0.8; util.sicker <- 0.5
  cost.sick   <- 3000; cost.sicker <- 7000; cost.int <- 1000
  coef_noint  <- log(0.2); HR_int <- 0.8
  drc <- 0.035; drq <- 0.035
  random_seed_sicker_i <- sample.int(100000, npats, replace = FALSE)
  beds              <- resource_discrete(650)
  beds_free         <- beds$n_free()
  shared_accumulator <- shared_input(0)
  value_accum       <- shared_accumulator$value()
})

cssd_common_pt <- add_item(death = max(0.0000001, rnorm(n = 1, mean = 12, sd = 3)))

cssd_unique_pt <- add_item(
  fl.sick = 1,
  q_default = util.sick,
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
  add_reactevt("sick", input = {
    shared_accumulator <- shared_accumulator$modify(shared_accumulator$value() + 1)
    value_accum <- shared_accumulator$value()
    beds_free <- beds$n_free(); time_in_queue <- NA
  }) |>
  add_reactevt("sicker", input = {
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
  add_reactevt("death", input = {
    beds$attempt_free()
    if (success_blocking_bed & beds$queue_size() > 0) {
      new_event(c(sicker = curtime), cur_evtlist, patient_id = beds$next_patient_in_line())
    }
    success_blocking_bed <- FALSE; time_in_queue <- NA
    beds_free <- beds$n_free(); q_default <- 0; c_default <- 0; curtime <- Inf
  })

res_cssd <- run_sim(
  npats = 200, n_sim = 1, psa_bool = FALSE,
  arm_list = c("int", "noint"),
  common_all_inputs = cssd_common_all,
  common_pt_inputs  = cssd_common_pt,
  unique_pt_inputs  = cssd_unique_pt,
  init_event_list   = cssd_events,
  evt_react_list    = cssd_reactions,
  util_ongoing_list = "q_default",
  cost_ongoing_list = "c_default",
  constrained = TRUE, ipd = 1,
  input_out = c("beds_free","had_to_queue","time_in_queue","value_accum"),
  seed = 42
)
out_c <- res_cssd[[1]][[1]]
det_c <- summary_results_det(out_c)
cat(sprintf("CONSTRAINED  total_costs: %s  total_qalys: %s\n", fv(unname(out_c$total_costs)), fv(unname(out_c$total_qalys))))
cat(sprintf("  total_costs_undisc: %s\n", fv(unname(out_c$total_costs_undisc))))
cat(sprintf("  total_qalys_undisc: %s\n", fv(unname(out_c$total_qalys_undisc))))
cat(sprintf("  dcosts[2]: %s  dqalys[2]: %s  ICUR[2]: %s\n", fs(det_c$dcosts[2]), fs(det_c$dqalys[2]), fs(det_c$ICUR[2])))

# Unconstrained
res_unc <- run_sim(
  npats = 200, n_sim = 1, psa_bool = FALSE,
  arm_list = c("int", "noint"),
  common_all_inputs = cssd_common_all,
  common_pt_inputs  = cssd_common_pt,
  unique_pt_inputs  = cssd_unique_pt,
  init_event_list   = cssd_events,
  evt_react_list    = cssd_reactions,
  util_ongoing_list = "q_default",
  cost_ongoing_list = "c_default",
  constrained = FALSE, ipd = 1,
  input_out = c("beds_free","had_to_queue","time_in_queue"),
  seed = 42
)
out_u <- res_unc[[1]][[1]]
det_u <- summary_results_det(out_u)
cat(sprintf("UNCONSTRAINED total_costs: %s  total_qalys: %s\n", fv(unname(out_u$total_costs)), fv(unname(out_u$total_qalys))))
cat(sprintf("  dcosts[2]: %s  dqalys[2]: %s  ICUR[2]: %s\n", fs(det_u$dcosts[2]), fs(det_u$dqalys[2]), fs(det_u$ICUR[2])))

# Unbinding constraint (1000 beds)
cssd_common_all_ub <- add_item(input = {
  util.sick <- 0.8; util.sicker <- 0.5
  cost.sick <- 3000; cost.sicker <- 7000; cost.int <- 1000
  coef_noint <- log(0.2); HR_int <- 0.8
  drc <- 0.035; drq <- 0.035
  random_seed_sicker_i <- sample.int(100000, npats, replace = FALSE)
  beds <- resource_discrete(1000)
  beds_free <- beds$n_free()
  shared_accumulator <- shared_input(0)
  value_accum <- shared_accumulator$value()
})
res_ub <- run_sim(
  npats = 200, n_sim = 1, psa_bool = FALSE,
  arm_list = c("int", "noint"),
  common_all_inputs = cssd_common_all_ub,
  common_pt_inputs  = cssd_common_pt,
  unique_pt_inputs  = cssd_unique_pt,
  init_event_list   = cssd_events,
  evt_react_list    = cssd_reactions,
  util_ongoing_list = "q_default",
  cost_ongoing_list = "c_default",
  constrained = TRUE, ipd = 1,
  input_out = c("beds_free","had_to_queue","time_in_queue"),
  seed = 42
)
out_ub <- res_ub[[1]][[1]]
det_ub <- summary_results_det(out_ub)
cat(sprintf("UNBINDING    total_costs: %s  total_qalys: %s\n", fv(unname(out_ub$total_costs)), fv(unname(out_ub$total_qalys))))
cat(sprintf("  dcosts[2]: %s  dqalys[2]: %s  ICUR[2]: %s\n", fs(det_ub$dcosts[2]), fs(det_ub$dqalys[2]), fs(det_ub$ICUR[2])))

# Constrained DSA
cssd_list_par <- list(
  parameter_name = list("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),
  base_value = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
  DSA_min    = list(0.6, 0.3, 1000, 5000,  800, log(0.1), 0.5),
  DSA_max    = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9),
  PSA_dist   = list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
  a          = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)),
  b          = lapply(list(0.8,0.5,3000,7000,1000,log(0.2),log(0.8)), \(x) abs(x/5))
)

cssd_sens_inputs <- add_item(
  indicators = if (sensitivity_bool) {
    create_indicators(sens, n_sensitivity * length(sensitivity_names), rep(1, length(cssd_list_par[[1]])))
  } else {
    rep(1, length(cssd_list_par[[1]]))
  }
)

cssd_dsa_common_all <- add_item(
  pick_val_v(
    base     = cssd_list_par[["base_value"]],
    psa      = pick_psa(cssd_list_par[["PSA_dist"]], rep(1, length(cssd_list_par[["PSA_dist"]])), cssd_list_par[["a"]], cssd_list_par[["b"]]),
    sens     = cssd_list_par[[sens_name_used]],
    psa_ind  = psa_bool, sens_ind = sensitivity_bool,
    indicator = indicators, names_out = cssd_list_par[["parameter_name"]]
  )
) |>
  add_item(input = {
    random_seed_sicker_i = sample(1:1000, 1000, replace = FALSE)
    beds <- resource_discrete(1000)
    beds_free <- beds$n_free()
    shared_accumulator <- shared_input(0)
    value_accum <- shared_accumulator$value()
  })

cat("\n=== 7. SSD CONSTRAINED DSA (npats=50, n_sim=1, seed=42) ===\n")
res_cdsa <- run_sim(
  npats = 50, n_sim = 1, psa_bool = FALSE,
  arm_list = c("int", "noint"),
  common_all_inputs = cssd_dsa_common_all,
  common_pt_inputs  = cssd_common_pt,
  unique_pt_inputs  = cssd_unique_pt,
  init_event_list   = cssd_events,
  evt_react_list    = cssd_reactions,
  util_ongoing_list = "q_default",
  cost_ongoing_list = "c_default",
  sensitivity_inputs = cssd_sens_inputs,
  sensitivity_names  = c("DSA_min","DSA_max"),
  sensitivity_bool   = TRUE,
  n_sensitivity      = length(cssd_list_par[[1]]),
  constrained = TRUE,
  input_out = unlist(cssd_list_par[["parameter_name"]]),
  seed = 42
)
for (si in c(1, 2, length(cssd_list_par[[1]]) + 1)) {
  o <- res_cdsa[[si]][[1]]
  cat(sprintf("sens=%d total_costs: %s  total_qalys: %s\n",
              si, fv(unname(o$total_costs)), fv(unname(o$total_qalys))))
}

cat("\n=== 8. SSD CONSTRAINED PROB DSA (npats=50, n_sim=3, seed=42) ===\n")
res_cpdsa <- run_sim(
  npats = 50, n_sim = 3, psa_bool = TRUE,
  arm_list = c("int", "noint"),
  common_all_inputs = cssd_dsa_common_all,
  common_pt_inputs  = cssd_common_pt,
  unique_pt_inputs  = cssd_unique_pt,
  init_event_list   = cssd_events,
  evt_react_list    = cssd_reactions,
  util_ongoing_list = "q_default",
  cost_ongoing_list = "c_default",
  sensitivity_inputs = cssd_sens_inputs,
  sensitivity_names  = c("DSA_min","DSA_max"),
  sensitivity_bool   = TRUE,
  n_sensitivity      = length(cssd_list_par[[1]]),
  constrained = TRUE,
  input_out = c(unlist(cssd_list_par[["parameter_name"]]), "beds_free"),
  seed = 42
)
for (si in c(1, 2)) {
  cat(sprintf("sens=%d mean_costs: %s  mean_qalys: %s\n",
              si,
              fv(sim_mean(res_cpdsa[[si]], "total_costs")),
              fv(sim_mean(res_cpdsa[[si]], "total_qalys"))))
}

cat("\n=== 9. SSD CONSTRAINED PSA (npats=50, n_sim=5, seed=42) ===\n")
cssd_common_all_psa <- add_item(input = {
  util.sick <- 0.8; util.sicker <- 0.5
  cost.sick <- 3000; cost.sicker <- 7000; cost.int <- 1000
  coef_noint <- log(0.2); HR_int <- 0.8
  drc <- 0.035; drq <- 0.035
  random_seed_sicker_i <- sample.int(100000, npats, replace = FALSE)
  beds <- resource_discrete(65)
  beds_free <- beds$n_free()
  shared_accumulator <- shared_input(0)
  value_accum <- shared_accumulator$value()
})
res_cpsa <- run_sim(
  npats = 50, n_sim = 5, psa_bool = TRUE,
  arm_list = c("int", "noint"),
  common_all_inputs = cssd_common_all_psa,
  common_pt_inputs  = cssd_common_pt,
  unique_pt_inputs  = cssd_unique_pt,
  init_event_list   = cssd_events,
  evt_react_list    = cssd_reactions,
  util_ongoing_list = "q_default",
  cost_ongoing_list = "c_default",
  sensitivity_inputs = cssd_sens_inputs,
  sensitivity_bool   = FALSE, n_sensitivity = 1,
  constrained = TRUE,
  input_out = c(unlist(cssd_list_par[["parameter_name"]]), "beds_free"),
  seed = 42
)
cat(sprintf("mean_costs: %s\nmean_qalys: %s\n",
            fv(sim_mean(res_cpsa[[1]], "total_costs")),
            fv(sim_mean(res_cpsa[[1]], "total_qalys"))))

# ==============================================================
# 10-16. INPUTS SELECTOR
# ==============================================================
cat("\n=== 10. INPUTS_SELECTOR MODELS ===\n")

l_inputs <- list(
  parameter_name = list("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),
  base_value  = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
  PSA_dist    = list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
  a           = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)),
  b           = lapply(list(0.8,0.5,3000,7000,1000,log(0.2),log(0.8)), \(x) abs(x/5)),
  n           = as.list(rep(1, 7)),
  DSA_min     = list(0.6, 0.3, 1000, 5000,  800, log(0.1), 0.5),
  DSA_max     = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9),
  scenario_1  = list(0.6, 0.3, 1000, 5000,  800, log(0.1), 0.5),
  scenario_2  = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9),
  psa_indicators = as.list(c(rep(1,4), rep(0,3)))
)

is_sens_inputs <- add_item(
  iterator_sensitivity = sens_iterator(sens, n_sensitivity)
) |>
  add_item(
    indicators = if (sensitivity_bool & sens_name_used %in% c("DSA_min","DSA_max")) {
      create_indicators(iterator_sensitivity, n_sensitivity * length(sensitivity_names), rep(1, length(l_inputs[[1]])))
    } else { rep(1, length(l_inputs[[1]])) }
  )

is_simple <- add_item() |>
  add_item(
    pick_val_v(
      base = l_inputs[["base_value"]],
      psa  = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]], l_inputs[["a"]], l_inputs[["b"]]),
      sens = l_inputs[[sens_name_used]],
      psa_ind = psa_bool, sens_ind = sensitivity_bool,
      indicator = indicators, names_out = l_inputs[["parameter_name"]],
      indicator_psa = l_inputs[["psa_indicators"]]
    )
  )

is_arm <- add_item(
  q_default = util.sick,
  c_default = cost.sick + if (arm == "int") { cost.int } else { 0 }
)

is_events <- add_tte(arm = c("noint","int"), evts = c("a1","b1"), input = { a1 <- 0; b1 <- 2 })
is_reactions <- add_reactevt("a1", input = {}) |>
  add_reactevt("b1", input = { q_default = 0; c_default = 0; curtime = Inf })

# Deterministic
res_is_det <- run_sim(npats=5, n_sim=1, psa_bool=FALSE, arm_list=c("int","noint"),
  common_all_inputs=is_simple, unique_pt_inputs=is_arm,
  init_event_list=is_events, evt_react_list=is_reactions,
  util_ongoing_list="q_default", cost_ongoing_list="c_default", ipd=1, seed=42)
cat(sprintf("IS DET   total_costs: %s  total_qalys: %s\n",
            fv(unname(res_is_det[[1]][[1]]$total_costs)),
            fv(unname(res_is_det[[1]][[1]]$total_qalys))))

# PSA only
res_is_psa <- run_sim(npats=5, n_sim=2, psa_bool=TRUE, arm_list=c("int","noint"),
  common_all_inputs=is_simple, unique_pt_inputs=is_arm,
  init_event_list=is_events, evt_react_list=is_reactions,
  util_ongoing_list="q_default", cost_ongoing_list="c_default", ipd=1,
  sensitivity_inputs=is_sens_inputs, sensitivity_names=NULL,
  sensitivity_bool=FALSE, n_sensitivity=1, seed=42)
cat(sprintf("IS PSA   mean_costs: %s  mean_qalys: %s\n",
            fv(sim_mean(res_is_psa[[1]], "total_costs")),
            fv(sim_mean(res_is_psa[[1]], "total_qalys"))))

# DSA
res_is_dsa <- run_sim(npats=5, n_sim=2, psa_bool=TRUE, arm_list=c("int","noint"),
  common_all_inputs=is_simple, unique_pt_inputs=is_arm,
  init_event_list=is_events, evt_react_list=is_reactions,
  util_ongoing_list="q_default", cost_ongoing_list="c_default", ipd=1,
  sensitivity_inputs=is_sens_inputs,
  sensitivity_names=c("DSA_min","DSA_max"), sensitivity_bool=TRUE,
  n_sensitivity=length(l_inputs[[1]]),
  input_out=unlist(l_inputs[["parameter_name"]]), seed=42)
for (si in c(1, 2, length(l_inputs[[1]])+1)) {
  cat(sprintf("IS DSA sens=%d mean_costs: %s  mean_qalys: %s\n",
              si,
              fv(sim_mean(res_is_dsa[[si]], "total_costs")),
              fv(sim_mean(res_is_dsa[[si]], "total_qalys"))))
}
# also print mean parameter values for sens=2 (util.sick should be 0.6)
merged2 <- bind_rows(map(res_is_dsa[[2]], "merged_df"))
cat(sprintf("IS DSA sens=2 mean(util.sick)=%.4f  mean(cost.sick)=%.2f\n",
            mean(merged2$util.sick, na.rm=TRUE),
            mean(merged2$cost.sick, na.rm=TRUE)))

# Scenario
res_is_scen <- run_sim(npats=5, n_sim=2, psa_bool=TRUE, arm_list=c("int","noint"),
  common_all_inputs=is_simple, unique_pt_inputs=is_arm,
  init_event_list=is_events, evt_react_list=is_reactions,
  util_ongoing_list="q_default", cost_ongoing_list="c_default", ipd=1,
  sensitivity_inputs=is_sens_inputs,
  sensitivity_names=c("scenario_1","scenario_2"), sensitivity_bool=TRUE,
  n_sensitivity=1,
  input_out=unlist(l_inputs[["parameter_name"]]), seed=42)
for (si in c(1, 2, 3)) {
  cat(sprintf("IS SCEN sens=%d mean_costs: %s  mean_qalys: %s\n",
              si,
              fv(sim_mean(res_is_scen[[si]], "total_costs")),
              fv(sim_mean(res_is_scen[[si]], "total_qalys"))))
}

# Params split across levels (DSA)
l_inputs_pat <- list(
  parameter_name = list("age","sex"),
  base_value = list(60, 1),
  PSA_dist   = list("rnorm","rbinom"),
  a = list(60, 1), b = list(10, 0.5), n = as.list(rep(1,2)),
  DSA_min = list(30, 0), DSA_max = list(80, 1),
  scenario_1 = list(55, 1), scenario_2 = list(45, 0),
  psa_indicators = as.list(rep(1,2))
)
is_sens_split <- add_item(iterator_sensitivity = sens_iterator(sens, n_sensitivity)) |>
  add_item(
    indicators = if (sensitivity_bool & sens_name_used %in% c("DSA_min","DSA_max")) {
      create_indicators(iterator_sensitivity, n_sensitivity*length(sensitivity_names), rep(1,length(l_inputs[[1]])))
    } else { rep(1,length(l_inputs[[1]])) }
  ) |>
  add_item(
    indicators_pat = if (sensitivity_bool & sens_name_used %in% c("DSA_min","DSA_max")) {
      create_indicators(iterator_sensitivity, n_sensitivity*length(sensitivity_names),
                        rep(1,length(l_inputs_pat[[1]])), length(l_inputs[[1]]))
    } else { rep(1,length(l_inputs_pat[[1]])) }
  )
is_pat <- add_item() |>
  add_item(pick_val_v(
    base=l_inputs_pat[["base_value"]],
    psa=pick_psa(l_inputs_pat[["PSA_dist"]],l_inputs_pat[["n"]],l_inputs_pat[["a"]],l_inputs_pat[["b"]]),
    sens=l_inputs_pat[[sens_name_used]], psa_ind=psa_bool, sens_ind=sensitivity_bool,
    indicator=indicators_pat, names_out=l_inputs_pat[["parameter_name"]],
    indicator_psa=l_inputs_pat[["psa_indicators"]]
  ))
res_is_split <- run_sim(npats=5, n_sim=2, psa_bool=FALSE, arm_list=c("int","noint"),
  common_all_inputs=is_simple, unique_pt_inputs=is_arm, common_pt_inputs=is_pat,
  init_event_list=is_events, evt_react_list=is_reactions,
  util_ongoing_list="q_default", cost_ongoing_list="c_default", ipd=1,
  sensitivity_inputs=is_sens_split,
  sensitivity_names=c("DSA_min","DSA_max"), sensitivity_bool=TRUE,
  n_sensitivity=length(l_inputs[[1]])+length(l_inputs_pat[[1]]),
  input_out=c(unlist(l_inputs[["parameter_name"]]),unlist(l_inputs_pat[["parameter_name"]])),
  seed=42)
for (si in c(1, 2, length(l_inputs[[1]])+1)) {
  cat(sprintf("IS SPLIT sens=%d mean_costs: %s  mean_qalys: %s\n",
              si,
              fv(sim_mean(res_is_split[[si]], "total_costs")),
              fv(sim_mean(res_is_split[[si]], "total_qalys"))))
}
merged_s2 <- bind_rows(map(res_is_split[[2]], "merged_df"))
cat(sprintf("IS SPLIT sens=2 mean(util.sick)=%.4f  mean(age)=%.2f\n",
            mean(merged_s2$util.sick, na.rm=TRUE), mean(merged_s2$age, na.rm=TRUE)))

# Covaried parameters
l_inputs_cov <- list(
  parameter_name = list("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),
  base_value  = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), 0.8),
  PSA_dist    = list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
  a           = list(0.8, 0.5, 3000, 7000, 1000, log(0.2), log(0.8)),
  b           = lapply(list(0.8,0.5,3000,7000,1000,log(0.2),log(0.8)), \(x) abs(x/5)),
  n           = as.list(rep(1,7)),
  DSA_min     = list(0.6, 0.3, 1000, 5000, 800, log(0.1), 0.5),
  DSA_max     = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9),
  scenario_1  = list(0.6, 0.3, 1000, 5000, 800, log(0.1), 0.5),
  scenario_2  = list(0.9, 0.7, 5000, 9000, 2000, log(0.4), 0.9),
  psa_indicators = as.list(c(rep(1,4),rep(0,3))),
  dsa_indicators = as.list(c(1,1,2,2,2,3,4))
)
l_inputs_pat_cov <- list(
  parameter_name = list("age","sex"),
  base_value = list(60,1),
  PSA_dist = list("rnorm","rbinom"),
  a=list(60,1), b=list(10,0.5), n=as.list(rep(1,2)),
  DSA_min=list(30,0), DSA_max=list(80,1),
  scenario_1=list(55,1), scenario_2=list(45,0),
  psa_indicators=as.list(rep(1,2)),
  dsa_indicators=list(5,5)
)
is_sens_cov <- add_item(iterator_sensitivity = sens_iterator(sens, n_sensitivity))
is_simple_cov <- add_item() |>
  add_item(pick_val_v(
    base=l_inputs_cov[["base_value"]],
    psa=pick_psa(l_inputs_cov[["PSA_dist"]],l_inputs_cov[["n"]],l_inputs_cov[["a"]],l_inputs_cov[["b"]]),
    sens=l_inputs_cov[[sens_name_used]], psa_ind=psa_bool, sens_ind=sensitivity_bool,
    indicator=l_inputs_cov[["dsa_indicators"]], sens_iterator=iterator_sensitivity,
    indicator_sens_binary=FALSE, names_out=l_inputs_cov[["parameter_name"]],
    indicator_psa=l_inputs_cov[["psa_indicators"]], distributions=l_inputs_cov[["PSA_dist"]], covariances=l_inputs_cov[["b"]]
  ))
is_pat_cov <- add_item() |>
  add_item(pick_val_v(
    base=l_inputs_pat_cov[["base_value"]],
    psa=pick_psa(l_inputs_pat_cov[["PSA_dist"]],l_inputs_pat_cov[["n"]],l_inputs_pat_cov[["a"]],l_inputs_pat_cov[["b"]]),
    sens=l_inputs_pat_cov[[sens_name_used]], psa_ind=psa_bool, sens_ind=sensitivity_bool,
    indicator=l_inputs_pat_cov[["dsa_indicators"]], sens_iterator=iterator_sensitivity,
    indicator_sens_binary=FALSE, names_out=l_inputs_pat_cov[["parameter_name"]],
    indicator_psa=l_inputs_pat_cov[["psa_indicators"]], distributions=l_inputs_pat_cov[["PSA_dist"]], covariances=l_inputs_pat_cov[["b"]]
  ))
n_cov_sens <- length(unique(unlist(l_inputs_cov[["dsa_indicators"]]))) + length(unique(unlist(l_inputs_pat_cov[["dsa_indicators"]])))
res_is_cov <- run_sim(npats=5, n_sim=2, psa_bool=FALSE, arm_list=c("int","noint"),
  common_all_inputs=is_simple_cov, unique_pt_inputs=is_arm, common_pt_inputs=is_pat_cov,
  init_event_list=is_events, evt_react_list=is_reactions,
  util_ongoing_list="q_default", cost_ongoing_list="c_default", ipd=1,
  sensitivity_inputs=is_sens_cov,
  sensitivity_names=c("DSA_min","DSA_max"), sensitivity_bool=TRUE,
  n_sensitivity=n_cov_sens,
  input_out=c(unlist(l_inputs_cov[["parameter_name"]]),unlist(l_inputs_pat_cov[["parameter_name"]])),
  seed=42)
for (si in c(1, 2, 3)) {
  cat(sprintf("IS COV sens=%d mean_costs: %s  mean_qalys: %s\n",
              si,
              fv(sim_mean(res_is_cov[[si]], "total_costs")),
              fv(sim_mean(res_is_cov[[si]], "total_qalys"))))
}
merged_c2 <- bind_rows(map(res_is_cov[[2]], "merged_df"))
cat(sprintf("IS COV sens=2 mean(util.sick)=%.4f  mean(util.sicker)=%.4f\n",
            mean(merged_c2$util.sick, na.rm=TRUE), mean(merged_c2$util.sicker, na.rm=TRUE)))

# Vector parameters
l_inputs_pat_vec <- list(
  parameter_name = list("age","sex","v_state"),
  base_value = list(60,1,c(10,20)),
  PSA_dist = list("rnorm","rbinom","mvrnorm"),
  a=list(60,1,c(10,20)), b=list(10,0.5,matrix(c(2,1,4,1),2,2)),
  n=as.list(rep(1,3)),
  DSA_min=list(30,0,c(5,10)), DSA_max=list(80,1,c(15,25)),
  scenario_1=list(55,1,c(12,21)), scenario_2=list(45,0,c(16,10)),
  psa_indicators=list(1,1,c(1,0)),
  dsa_indicators=list(5,5,c(6,6))
)
n_vec_sens <- length(unique(unlist(l_inputs_cov[["dsa_indicators"]]))) + length(unique(unlist(l_inputs_pat_vec[["dsa_indicators"]])))
is_pat_vec <- add_item() |>
  add_item(pick_val_v(
    base=l_inputs_pat_vec[["base_value"]],
    psa=pick_psa(l_inputs_pat_vec[["PSA_dist"]],l_inputs_pat_vec[["n"]],l_inputs_pat_vec[["a"]],l_inputs_pat_vec[["b"]]),
    sens=l_inputs_pat_vec[[sens_name_used]], psa_ind=psa_bool, sens_ind=sensitivity_bool,
    indicator=l_inputs_pat_vec[["dsa_indicators"]], sens_iterator=iterator_sensitivity,
    indicator_sens_binary=FALSE, names_out=l_inputs_pat_vec[["parameter_name"]],
    indicator_psa=l_inputs_pat_vec[["psa_indicators"]], distributions=l_inputs_pat_vec[["PSA_dist"]], covariances=l_inputs_pat_vec[["b"]]
  ))
res_is_vec <- run_sim(npats=5, n_sim=2, psa_bool=FALSE, arm_list=c("int","noint"),
  common_all_inputs=is_simple_cov, unique_pt_inputs=is_arm, common_pt_inputs=is_pat_vec,
  init_event_list=is_events, evt_react_list=is_reactions,
  util_ongoing_list="q_default", cost_ongoing_list="c_default", ipd=1,
  sensitivity_inputs=is_sens_cov,
  sensitivity_names=c("DSA_min","DSA_max"), sensitivity_bool=TRUE,
  n_sensitivity=n_vec_sens,
  input_out=c(unlist(l_inputs_cov[["parameter_name"]]),unlist(l_inputs_pat_vec[["parameter_name"]][1:2])),
  seed=42)
for (si in c(1, 2, 3)) {
  cat(sprintf("IS VEC sens=%d mean_costs: %s  mean_qalys: %s\n",
              si,
              fv(sim_mean(res_is_vec[[si]], "total_costs")),
              fv(sim_mean(res_is_vec[[si]], "total_qalys"))))
}
# check v_state values for base case (sens=1, sim=1)
v_state_s1 <- res_is_vec[[1]][[1]]$v_state
cat(sprintf("IS VEC sens=1 v_state int: %s  noint: %s\n",
            fv(colMeans(v_state_s1$int)), fv(colMeans(v_state_s1$noint))))
# check v_state for DSA_min sens=6 (v_state should be c(5,10))
v_state_s6 <- res_is_vec[[6]][[1]]$v_state
cat(sprintf("IS VEC sens=6 v_state int: %s  noint: %s\n",
            fv(colMeans(v_state_s6$int)), fv(colMeans(v_state_s6$noint))))

cat("\n=== DONE ===\n")
