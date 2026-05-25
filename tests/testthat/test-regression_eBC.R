# Early Breast Cancer (eBC) regression test — deterministic only

test_that("eBC regression - deterministic (npats=200, seed=42)", {
  skip_on_cran()

  ebc_common_all <- add_item(input = {
    pick_val_v(
      base = c(0.75, 0.8, 0.9, 0.7, 0.6, 0.5),
      psa  = MASS::mvrnorm(1, c(0.75, 0.8, 0.9, 0.7, 0.6, 0.5), diag(rep(0.0004, 6))),
      sens = c(0.75, 0.8, 0.9, 0.7, 0.6, 0.5),
      psa_ind = psa_bool, sens_ind = sensitivity_bool,
      indicator = rep(0, 6L),
      names_out = c("util.idfs.ontx","util.idfs.offtx","util.remission",
                    "util.recurrence","util.mbc.progression.mbc","util.mbc.pps"),
      deploy_env = TRUE
    )
    pick_val_v(
      base = c(40000, 5000, 3000, 10000, 30000, 10000, 20000, 30000, 20000, 1000),
      psa  = rgamma_mse(1, c(40000, 5000, 3000, 10000, 30000, 10000, 20000, 30000, 20000, 1000),
                           c(8000, 1000, 600, 2000, 6000, 2000, 4000, 6000, 4000, 200)),
      sens = c(40000, 5000, 3000, 10000, 30000, 10000, 20000, 30000, 20000, 1000),
      psa_ind = psa_bool, sens_ind = sensitivity_bool,
      indicator = rep(0, 10L),
      names_out = c("cost.idfs.tx","cost.recurrence","cost.mbc.tx","cost.tx.beva",
                    "cost.idfs.txnoint","cost.idfs","cost.mbc.progression.mbc",
                    "cost.mbc.pps","cost.2ndline","cost.ae"),
      deploy_env = TRUE
    )
    .q <- quote(
      if (fl.idfs==1) {
        util.idfs.ontx*fl.idfs.ontx+(1-fl.idfs.ontx)*(1-fl.idfs.ontx)
      } else if (fl.idfs==0&fl.mbcs==0) {
        util.remission*fl.remission+fl.recurrence*util.recurrence
      } else if (fl.mbcs==1) {
        util.mbc.progression.mbc*fl.mbcs.progression.mbc+
          (1-fl.mbcs.progression.mbc)*util.mbc.pps
      }
    )
    .cmbc <- quote(
      cost.mbc.tx*fl.mbcs.ontx+cost.mbc.progression.mbc*fl.mbcs.progression.mbc+
        cost.mbc.pps*(1-fl.mbcs.progression.mbc)+cost.2ndline*fl.mbcs_2ndline
    )
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
    fl.idfs.ontx <- 1; fl.idfs <- 1; fl.mbcs.ontx <- 1
    fl.mbcs.progression.mbc <- 1; fl.tx.beva <- 1; fl.mbcs <- 0
    fl.mbcs_2ndline <- 0; fl.recurrence <- 0
    q_default <- if (fl.idfs == 1) {
      util.idfs.ontx * fl.idfs.ontx + (1 - fl.idfs.ontx) * (1 - fl.idfs.ontx)
    } else if (fl.idfs == 0 & fl.mbcs == 0) {
      util.remission * fl.remission + fl.recurrence * util.recurrence
    } else if (fl.mbcs == 1) {
      util.mbc.progression.mbc * fl.mbcs.progression.mbc +
        (1 - fl.mbcs.progression.mbc) * util.mbc.pps
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
            evts = c("start","ttot","ttot.beva","progression.mbc","os","idfs",
                     "ttot.early","remission","recurrence","start.early.mbc","ae","2ndline_mbc"),
            other_inp = c("os.early","os.mbc"), input = {
              start <- 0
              idfs  <- draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2))
              ttot.early <- min(draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2)), idfs)
              ttot.beva  <- draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2))
              os.early   <- draw_tte(1, "lnorm", coef1 = 3, coef2 = log(0.2))
              if (fl.remission) {
                recurrence <- idfs + draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2))
                remission  <- idfs
                if (min(os.early, nat.os.s) > recurrence) {
                  os.mbc <- draw_tte(1,"lnorm",coef1=0.8,coef2=log(0.2))+idfs+recurrence
                  progression.mbc <- draw_tte(1,"lnorm",coef1=0.5,coef2=log(0.2))+idfs+recurrence
                  ttot <- draw_tte(1,"lnorm",coef1=0.5,coef2=log(0.2))+idfs+recurrence
                }
              } else {
                start.early.mbc <- draw_tte(1, "lnorm", coef1 = 2.3, coef2 = log(0.2))
                idfs <- ifelse(start.early.mbc < idfs, start.early.mbc, idfs)
                ttot.early <- min(ifelse(start.early.mbc < idfs, start.early.mbc, idfs), ttot.early)
                os.mbc <- draw_tte(1,"lnorm",coef1=0.8,coef2=log(0.2))+start.early.mbc
                progression.mbc <- draw_tte(1,"lnorm",coef1=0.5,coef2=log(0.2))+start.early.mbc
                ttot <- draw_tte(1,"lnorm",coef1=0.5,coef2=log(0.2))+start.early.mbc
              }
              os <- min(os.mbc, os.early, nat.os.s)
            }) |>
    add_tte(arm = "noint",
            evts = c("start","ttot","ttot.beva","progression.mbc","os","idfs",
                     "ttot.early","remission","recurrence","start.early.mbc"),
            other_inp = c("os.early","os.mbc"), input = {
              start <- 0
              idfs  <- draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2), beta_tx = 1.2)
              ttot.early <- min(draw_tte(1,"lnorm",coef1=2,coef2=log(0.2),beta_tx=1.2), idfs)
              os.early   <- draw_tte(1, "lnorm", coef1 = 3, coef2 = log(0.2), beta_tx = 1.2)
              if (fl.remission) {
                recurrence <- idfs + draw_tte(1, "lnorm", coef1 = 2, coef2 = log(0.2))
                remission  <- idfs
                if (min(os.early, nat.os.s) > recurrence) {
                  os.mbc <- draw_tte(1,"lnorm",coef1=0.8,coef2=log(0.2))+idfs+recurrence
                  progression.mbc <- draw_tte(1,"lnorm",coef1=0.5,coef2=log(0.2))+idfs+recurrence
                  ttot <- draw_tte(1,"lnorm",coef1=0.5,coef2=log(0.2))+idfs+recurrence
                }
              } else {
                start.early.mbc <- draw_tte(1, "lnorm", coef1 = 2.3, coef2 = log(0.2))
                idfs <- ifelse(start.early.mbc < idfs, start.early.mbc, idfs)
                ttot.early <- min(ifelse(start.early.mbc < idfs, start.early.mbc, idfs), ttot.early)
                os.mbc <- draw_tte(1,"lnorm",coef1=0.8,coef2=log(0.2))+start.early.mbc
                progression.mbc <- draw_tte(1,"lnorm",coef1=0.5,coef2=log(0.2))+start.early.mbc
                ttot <- draw_tte(1,"lnorm",coef1=0.5,coef2=log(0.2))+start.early.mbc
              }
              os <- min(os.mbc, os.early, nat.os.s)
            })

  ebc_reactions <-
    add_reactevt(name_evt = "start",       input = {}) |>
    add_reactevt(name_evt = "ttot",        input = { q_default <- eval(.q); c_default <- eval(.cmbc); fl.mbcs.ontx <- 0 }) |>
    add_reactevt(name_evt = "ttot.beva",   input = { q_default <- eval(.q); c_default <- eval(.cmbc); fl.tx.beva <- 0 }) |>
    add_reactevt(name_evt = "progression.mbc", input = {
      q_default <- eval(.q); c_default <- eval(.cmbc)
      fl.mbcs.progression.mbc <- 0; fl.mbcs_2ndline <- 1
      new_event(c("2ndline_mbc" = curtime + qexp(rnd_stream_mbc$draw_n(), 0.08) / 12))
    }) |>
    add_reactevt(name_evt = "idfs", input = {
      q_default = eval(.q)
      c_default <- if (arm == "noint") {
        cost.idfs.txnoint * fl.idfs.ontx + cost.idfs
      } else { cost.idfs.tx * fl.idfs.ontx + cost.tx.beva * fl.tx.beva + cost.idfs }
      fl.idfs <- 0
    }) |>
    add_reactevt(name_evt = "ttot.early", input = {
      q_default <- eval(.q)
      c_default <- if (arm == "noint") {
        cost.idfs.txnoint * fl.idfs.ontx + cost.idfs
      } else { cost.idfs.tx * fl.idfs.ontx + cost.tx.beva * fl.tx.beva + cost.idfs }
      fl.idfs.ontx <- 0; fl.tx.beva <- 0
      n_ae <- qpois(rnd_stream_ae$draw_n(), lambda = 0.25 * (curtime - prevtime))
      if (n_ae > 0) { new_event(c("ae" = curtime + 0.0001)); n_ae <- n_ae - 1 }
    }) |>
    add_reactevt(name_evt = "remission",      input = { q_default <- eval(.q); c_default <- cost.recurrence*fl.recurrence; fl.remission <- 1 }) |>
    add_reactevt(name_evt = "recurrence",     input = { q_default <- eval(.q); c_default <- cost.recurrence*fl.recurrence; fl.recurrence <- 1; fl.remission <- 0; fl.mbcs <- 1; fl.mbcs.progression.mbc <- 1 }) |>
    add_reactevt(name_evt = "start.early.mbc",input = { q_default <- eval(.q); c_default <- cost.recurrence*fl.recurrence; fl.mbcs <- 1; fl.mbcs.progression.mbc <- 1 }) |>
    add_reactevt(name_evt = "2ndline_mbc", input = {
      q_default <- eval(.q); c_default <- eval(.cmbc); fl.mbcs_2ndline <- 0
      n_ae <- qpois(rnd_stream_ae$draw_n(), lambda = 0.25 * (curtime - prevtime))
      if (n_ae > 0) { new_event(c("ae" = curtime + 0.0001)); n_ae <- n_ae - 1 }
    }) |>
    add_reactevt(name_evt = "ae", input = {
      if (n_ae > 0) { new_event(c("ae" = curtime)); n_ae <- n_ae - 1 }
      q_default = eval(.q); c_default <- eval(.cmbc); c_ae <- cost.ae
      modify_event(c("os" = max(get_event("os") - 0.125, curtime + 0.0001)))
    }) |>
    add_reactevt(name_evt = "os", input = {
      q_default <- eval(.q); c_default <- eval(.cmbc)
      fl.tx.beva <- 0; fl.mbcs.ontx <- 0; fl.idfs <- 0; fl.mbcs <- 0; curtime <- Inf
    })

  suppressMessages({
    res <- run_sim(
      npats = 200, n_sim = 1, psa_bool = FALSE, arm_list = c("int", "noint"),
      common_all_inputs = ebc_common_all, common_pt_inputs = ebc_common_pt,
      unique_pt_inputs  = ebc_unique_pt,  init_event_list   = ebc_events,
      evt_react_list    = ebc_reactions,
      util_ongoing_list = "q_default", cost_ongoing_list = "c_default",
      cost_instant_list = "c_ae",
      input_out = c("os.early","os.mbc","nat.os.s","sex_pt"), seed = 42
    )
    out <- res[[1]][[1]]
    tc  <- unname(out$total_costs)
    tq  <- unname(out$total_qalys)
    dc  <- tc[1] - tc[2]; dq <- tq[1] - tq[2]
  })

  expect_equal(tc, c(424780.04076015, 238112.40527392), tolerance = 1e-4)
  expect_equal(tq, c(11.08750972, 9.99000132),          tolerance = 1e-4)
  expect_equal(dc, 186667.63548623, tolerance = 1e-4)
  expect_equal(dq, 1.09750840,      tolerance = 1e-4)
  expect_equal(dc / dq, 170083.10402456, tolerance = 1e-3)
})
