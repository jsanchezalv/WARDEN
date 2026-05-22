suppressMessages({
test_that("Test minimal model runs with some basic settings", {

# No discount rate, no qalys nor costs --------------------------------------------------------

  
  common_all_inputs <-add_item(input = {
    drc         <- 0 
    drq         <- 0
  })

    init_event_list <- 
      add_tte(arm=c("aa","bb"), evts = c("start","death") ,input={
        start <- 0
        death <- 1
      })
    evt_react_list <-
      add_reactevt(name_evt = "start",
                   input = {}) %>%
      add_reactevt(name_evt = "death",
                   input = {
                     curtime   <- Inf
                   }) 
    
  
    
    results <- run_sim(  
      npats=1,                              
      n_sim=1,                                  
      psa_bool = FALSE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,          
      ipd = 1
    )
  
    expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
    expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
    expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
    expect_no_error(summary_results_det(results[[1]][[1]]))
    expect_no_error(summary_results_sim(results[[1]]))
    expect_no_error(summary_results_sens(results))
    
# No discount rate, qalys and costs = 0 --------------------------------------------------------
    
    common_all_inputs <-add_item(input = {
      drc         <- 0 
      drq         <- 0
      q_default <- 0
      c_default <- 0
    })
    util_ongoing <- "q_default"
    cost_ongoing <- "c_default"
    results <- run_sim(  
      npats=1,                              
      n_sim=1,                                  
      psa_bool = FALSE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,          
      util_ongoing_list = util_ongoing,
      cost_ongoing_list = cost_ongoing,
      ipd = 1
    )
    
    expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
    expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
    expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
    expect_no_error(summary_results_det(results[[1]][[1]]))
    expect_no_error(summary_results_sim(results[[1]]))
    expect_no_error(summary_results_sens(results))
    
    # No discount rate, qalys and costs != 0 --------------------------------------------------------
    
    common_all_inputs <-add_item(input = {
      drc         <- 0 
      drq         <- 0
      q_default <- 0.8
      c_default <- 10
    })
    
    results <- run_sim(  
      npats=10,                              
      n_sim=1,                                  
      psa_bool = FALSE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,          
      util_ongoing_list = util_ongoing,
      cost_ongoing_list = cost_ongoing,
      ipd = 1
    )
    expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
    expect_equal(unname(results[[1]][[1]]$total_qalys), c(0.8,0.8))
    expect_equal(unname(results[[1]][[1]]$total_costs), c(10,10))
    expect_no_error(summary_results_det(results[[1]][[1]]))
    expect_no_error(summary_results_sim(results[[1]]))
    expect_no_error(summary_results_sens(results))
    
    # Inf discount rate, qalys and costs != 0 --------------------------------------------------------
    
    common_all_inputs <-add_item(input = {
      drc         <- Inf 
      drq         <- Inf
      q_default <- 0.8
      c_default <- 10
    })
    
    results <- run_sim(  
      npats=10,                              
      n_sim=1,                                  
      psa_bool = FALSE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,          
      util_ongoing_list = util_ongoing,
      cost_ongoing_list = cost_ongoing,
      ipd = 1
    )
    expect_equal(unname(results[[1]][[1]]$total_lys), c(0,0))
    expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
    expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
    expect_no_error(summary_results_det(results[[1]][[1]]))
    expect_no_error(summary_results_sim(results[[1]]))
    expect_no_error(summary_results_sens(results))
    
    # Death is last event but does not set curtime to Inf --------------------------------------------------------
    
    common_all_inputs <-add_item(input = {
      drc         <- 0 
      drq         <- 0
    })
    
    evt_react_list <-
      add_reactevt(name_evt = "start",
                   input = {}) %>%
      add_reactevt(name_evt = "death",
                   input = {}) 
    
    results <- run_sim(  
      npats=10,                              
      n_sim=1,                                  
      psa_bool = FALSE,                         
      arm_list = c("aa", "bb"),             # intervention list
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,
      ipd = 1
    )
    expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
    expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
    expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
    expect_no_error(summary_results_det(results[[1]][[1]]))
    expect_no_error(summary_results_sim(results[[1]]))
    expect_no_error(summary_results_sens(results))

# Can export other inputs -------------------------------------------------

    
    results <- run_sim(  
      npats=10,                              
      n_sim=1,                                  
      psa_bool = FALSE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,
      ipd = 1, input_out = c("drc","drq")
    )
    expect_equal(unname(results[[1]][[1]]$drc), c(0,0))
    expect_equal(unname(results[[1]][[1]]$drq), c(0,0))
    expect_no_error(summary_results_det(results[[1]][[1]]))
    expect_no_error(summary_results_sim(results[[1]]))
    expect_no_error(summary_results_sens(results))

# Export works for length > 1 --------------------------------------------

    common_all_inputs <-add_item(input = {
      drc         <- 0 
      drq         <- 0
      vec <- list(c(0,5),c(0,5,6,7))
    })
    
    
    results <- run_sim(  
      npats=10,                              
      n_sim=1,                                  
      psa_bool = FALSE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,
      ipd = 1, input_out = c("drc","drq","vec")
    )
    expect_equal(unname(results[[1]][[1]]$vec), list(list(c(0,5),c(0,5,6,7)),list(c(0,5),c(0,5,6,7))))
    expect_equal(unname(results[[1]][[1]]$extradata_raw[[1]]$vec), list(c(0,5),c(0,5,6,7)))
    expect_no_error(summary_results_det(results[[1]][[1]]))
    expect_no_error(summary_results_sim(results[[1]]))
    expect_no_error(summary_results_sens(results))


    
    
# Many options set as TRUE ------------------------------------------------

    results <- run_sim(  
      npats=10,                              
      n_sim=1,                                  
      psa_bool = TRUE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,
      ipd = 1,
      input_out = c("drc","drq","vec"),
      timed_freq = TRUE, 
      accum_backwards = TRUE,
      continue_on_error = TRUE,
      seed = 2
    )
    expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
    expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
    expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
    expect_no_error(summary_results_det(results[[1]][[1]]))
    expect_no_error(summary_results_sim(results[[1]]))
    expect_no_error(summary_results_sens(results))
    

# Check continue on error works and inputs can change across simulations -------------------------------------------


    common_all_inputs <-add_item(input = {
      drc         <- 0 
      drq         <- ifelse(simulation == 3, Inf,0)
      err_flag <- if(simulation == 2){stop("ERROR")}else{NULL}
    })
    
    
    results <- run_sim(  
      npats=10,                              
      n_sim=3,                                  
      psa_bool = TRUE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,
      ipd = 2,
      input_out = c("drc","drq"),
      timed_freq = TRUE, 
      accum_backwards = TRUE,
      continue_on_error = TRUE,
      seed = 2
    )
    
    expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
    expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
    expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
    expect_equal(unname(results[[1]][[1]]$timed_outputs$total_lys$aa), c(0,1))
    expect_equal(results[[1]][[2]], NULL)
    expect_equal(unname(results[[1]][[3]]$total_lys), c(0,0))
    expect_equal(unname(results[[1]][[3]]$total_qalys), c(0,0))
    expect_equal(unname(results[[1]][[3]]$total_costs), c(0,0))
    expect_no_error(summary_results_det(results[[1]][[1]]))
    expect_no_error(summary_results_sim(results[[1]]))
    expect_no_error(summary_results_sens(results))
    

# several analyses --------------------------------------------------------
    
    
    common_all_inputs <-add_item(input = {
      drc         <- 0 
      drq         <- ifelse(simulation == 3, Inf,0)
      err_flag <- if(simulation == 2){stop("ERROR")}else{NULL}
    })
    
    
    results <- run_sim(  
      npats=10,                              
      n_sim=3,                                  
      psa_bool = TRUE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,
      ipd = 2,
      input_out = c("drc","drq"),
      timed_freq = TRUE, 
      accum_backwards = TRUE,
      continue_on_error = TRUE,
      seed = 2,
      n_sensitivity = 2,
      sensitivity_bool = TRUE,
      sensitivity_names = c("A","B")
    )
    
    expect_equal(unname(results[[4]][[1]]$total_lys), c(1,1))
    expect_equal(unname(results[[4]][[1]]$total_qalys), c(0,0))
    expect_equal(unname(results[[4]][[1]]$total_costs), c(0,0))
    expect_equal(unname(results[[4]][[1]]$timed_outputs$total_lys$aa), c(0,1))
    expect_equal(results[[4]][[2]], NULL)
    expect_equal(unname(results[[4]][[3]]$total_lys), c(0,0))
    expect_equal(unname(results[[4]][[3]]$total_qalys), c(0,0))
    expect_equal(unname(results[[4]][[3]]$total_costs), c(0,0))
    expect_no_error(summary_results_det(results[[4]][[1]]))
    expect_no_error(summary_results_sim(results[[4]]))
    expect_no_error(summary_results_sens(results))
    
    
    
    # Over time outputs make sense --------------------------------------------
    common_all_inputs <-add_item(input = {
      drc         <- 0 
      drq         <- 0
      q_default <- 0.8
      
      cost_a <- 1
      cost_a_cycle_starttime <- 0
      cost_a_cycle_l <- 0.125
      cost_a_max_cycles <- 100
    })
    
    init_event_list <- 
      add_tte(arm=c("aa","bb"), evts = c("start","cycle","death") ,input={
        start <- 0
        death <- 10.1
        cycle <- 1
      })
    evt_react_list <-
      add_reactevt(name_evt = "start",
                   input = {
                     cost_insta <- 100
                     
                   }) %>%
      add_reactevt(name_evt = "cycle",
                   input = {
                     if(curtime<9){
                       new_event(list(cycle = curtime + 1.1))
                     }
                     cost_insta <- 200
                     
                   }) %>%
      add_reactevt(name_evt = "death",
                   input = {
                   }) 
    
    results <- run_sim(  
      npats=10,                              
      n_sim=1,                                  
      psa_bool = TRUE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,
      util_ongoing_list = 'q_default',
      cost_instant_list = "cost_insta",
      cost_cycle_list = "cost_a",
      ipd = 1,
      timed_freq = 0.25,
    )
    
    expect_equal(results[[1]][[1]]$timed_outputs$total_lys$aa, c(seq(0,10.1, by = 0.25),10.1))
    expect_equal(results[[1]][[1]]$timed_outputs$timepoints, c(seq(0,10.1, by = 0.25),10.1))
    expect_equal(results[[1]][[1]]$timed_outputs$cost_insta$aa,
                 c(100, 100, 100, 100, 300, 300, 300, 300, 300, 500, 500, 500, 
                   500, 700, 700, 700, 700, 700, 900, 900, 900, 900, 1100, 1100, 
                   1100, 1100, 1300, 1300, 1300, 1300, 1300, 1500, 1500, 1500, 1500, 
                   1700, 1700, 1700, 1700, 1700, 1900, 1900))
    expect_equal(results[[1]][[1]]$timed_outputs$cost_a$aa,
                 c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 
                   32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 
                   64, 66, 68, 70, 72, 74, 76, 78, 80, 81))
    expect_equal(results[[1]][[1]]$cost_a[["aa"]],max(results[[1]][[1]]$timed_outputs$cost_a$aa))
    expect_equal(results[[1]][[1]]$total_lys[["aa"]],max(results[[1]][[1]]$timed_outputs$total_lys$aa))
    expect_equal(results[[1]][[1]]$total_qalys[["aa"]],max(results[[1]][[1]]$timed_outputs$total_qalys$aa))
    
    
    
    # Over time outputs make sense, backwards --------------------------------------------
    common_all_inputs <-add_item(input = {
      drc         <- 0 
      drq         <- 0
      q_default <- 0.8
      
      cost_a <- 0
      cost_a_cycle_starttime <- 0
      cost_a_cycle_l <- 0.125
      cost_a_max_cycles <- 100
    })
    
    init_event_list <- 
      add_tte(arm=c("aa","bb"), evts = c("start","cycle","death") ,input={
        start <- 0
        death <- 10.1
        cycle <- 1
      })
    evt_react_list <-
      add_reactevt(name_evt = "start",
                   input = {
                     cost_insta <- 100
                     
                   }) %>%
      add_reactevt(name_evt = "cycle",
                   input = {
                     if(curtime<9){
                       new_event(list(cycle = curtime + 1.1))
                     }
                     cost_insta <- 200
                     cost_a <- 1
                     
                   }) %>%
      add_reactevt(name_evt = "death",
                   input = {
                   }) 
    
    results <- run_sim(  
      npats=10,                              
      n_sim=1,                                  
      psa_bool = TRUE,                         
      arm_list = c("aa", "bb"),             
      common_all_inputs = common_all_inputs,    
      init_event_list = init_event_list,        
      evt_react_list = evt_react_list,
      util_ongoing_list = 'q_default',
      cost_instant_list = "cost_insta",
      cost_cycle_list = "cost_a",
      ipd = 1,
      timed_freq = 0.25, accum_backwards = TRUE
    )
    
    expect_equal(results[[1]][[1]]$timed_outputs$total_lys$aa, c(seq(0,10.1, by = 0.25),10.1))
    expect_equal(results[[1]][[1]]$timed_outputs$timepoints, c(seq(0,10.1, by = 0.25),10.1))
    expect_equal(results[[1]][[1]]$timed_outputs$cost_insta$aa,
                 c(100, 100, 100, 100, 300, 300, 300, 300, 300, 500, 500, 500, 
                   500, 700, 700, 700, 700, 700, 900, 900, 900, 900, 1100, 1100, 
                   1100, 1100, 1300, 1300, 1300, 1300, 1300, 1500, 1500, 1500, 1500, 
                   1700, 1700, 1700, 1700, 1700, 1900, 1900))
    expect_equal(results[[1]][[1]]$timed_outputs$cost_a$aa,
                 c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 
                   32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 
                   64, 66, 68, 70, 72, 74, 76, 78, 80, 81))
    expect_equal(results[[1]][[1]]$cost_a[["aa"]],max(results[[1]][[1]]$timed_outputs$cost_a$aa))
    expect_equal(results[[1]][[1]]$total_lys[["aa"]],max(results[[1]][[1]]$timed_outputs$total_lys$aa))
    expect_equal(results[[1]][[1]]$total_qalys[["aa"]],max(results[[1]][[1]]$timed_outputs$total_qalys$aa))
    
    
})


test_that("Test everything but with constrained = TRUE", {
  
  
  # No discount rate, no qalys nor costs --------------------------------------------------------
  
  
  common_all_inputs <-add_item(input = {
    drc         <- 0 
    drq         <- 0
  })
  
  init_event_list <- 
    add_tte(arm=c("aa","bb"), evts = c("start","death") ,input={
      start <- 0
      death <- 1
    })
  evt_react_list <-
    add_reactevt(name_evt = "start",
                 input = {}) %>%
    add_reactevt(name_evt = "death",
                 input = {
                   curtime   <- Inf
                 }) 
  
  
  
  results <- run_sim(  
    npats=1,                              
    n_sim=1,                                  
    psa_bool = FALSE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,          
    ipd = 1,
    constrained = TRUE
  )
  
  expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
  expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
  expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
  expect_no_error(summary_results_det(results[[1]][[1]]))
  expect_no_error(summary_results_sim(results[[1]]))
  expect_no_error(summary_results_sens(results))
  
  # No discount rate, qalys and costs = 0 --------------------------------------------------------
  
  common_all_inputs <-add_item(input = {
    drc         <- 0 
    drq         <- 0
    q_default <- 0
    c_default <- 0
  })
  util_ongoing <- "q_default"
  cost_ongoing <- "c_default"
  results <- run_sim(  
    npats=1,                              
    n_sim=1,                                  
    psa_bool = FALSE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,          
    util_ongoing_list = util_ongoing,
    cost_ongoing_list = cost_ongoing,
    ipd = 1,
    constrained = TRUE
  )
  
  expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
  expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
  expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
  expect_no_error(summary_results_det(results[[1]][[1]]))
  expect_no_error(summary_results_sim(results[[1]]))
  expect_no_error(summary_results_sens(results))
  
  # No discount rate, qalys and costs != 0 --------------------------------------------------------
  
  common_all_inputs <-add_item(input = {
    drc         <- 0 
    drq         <- 0
    q_default <- 0.8
    c_default <- 10
  })
  
  results <- run_sim(  
    npats=10,                              
    n_sim=1,                                  
    psa_bool = FALSE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,          
    util_ongoing_list = util_ongoing,
    cost_ongoing_list = cost_ongoing,
    ipd = 1,
    constrained = TRUE
  )
  expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
  expect_equal(unname(results[[1]][[1]]$total_qalys), c(0.8,0.8))
  expect_equal(unname(results[[1]][[1]]$total_costs), c(10,10))
  expect_no_error(summary_results_det(results[[1]][[1]]))
  expect_no_error(summary_results_sim(results[[1]]))
  expect_no_error(summary_results_sens(results))
  
  # Inf discount rate, qalys and costs != 0 --------------------------------------------------------
  
  common_all_inputs <-add_item(input = {
    drc         <- Inf 
    drq         <- Inf
    q_default <- 0.8
    c_default <- 10
  })
  
  results <- run_sim(  
    npats=10,                              
    n_sim=1,                                  
    psa_bool = FALSE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,          
    util_ongoing_list = util_ongoing,
    cost_ongoing_list = cost_ongoing,
    ipd = 1,
    constrained = TRUE
  )
  expect_equal(unname(results[[1]][[1]]$total_lys), c(0,0))
  expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
  expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
  expect_no_error(summary_results_det(results[[1]][[1]]))
  expect_no_error(summary_results_sim(results[[1]]))
  expect_no_error(summary_results_sens(results))
  
  # Death is last event but does not set curtime to Inf --------------------------------------------------------
  
  common_all_inputs <-add_item(input = {
    drc         <- 0 
    drq         <- 0
  })
  
  evt_react_list <-
    add_reactevt(name_evt = "start",
                 input = {}) %>%
    add_reactevt(name_evt = "death",
                 input = {}) 
  
  results <- run_sim(  
    npats=10,                              
    n_sim=1,                                  
    psa_bool = FALSE,                         
    arm_list = c("aa", "bb"),             # intervention list
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,
    ipd = 1,
    constrained = TRUE
  )
  expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
  expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
  expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
  expect_no_error(summary_results_det(results[[1]][[1]]))
  expect_no_error(summary_results_sim(results[[1]]))
  expect_no_error(summary_results_sens(results))
  
  # Can export other inputs -------------------------------------------------
  
  
  results <- run_sim(  
    npats=10,                              
    n_sim=1,                                  
    psa_bool = FALSE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,
    ipd = 1, input_out = c("drc","drq"),
    constrained = TRUE
  )
  expect_equal(unname(results[[1]][[1]]$drc), c(0,0))
  expect_equal(unname(results[[1]][[1]]$drq), c(0,0))
  expect_no_error(summary_results_det(results[[1]][[1]]))
  expect_no_error(summary_results_sim(results[[1]]))
  expect_no_error(summary_results_sens(results))
  
  # Export works for length > 1 --------------------------------------------
  
  common_all_inputs <-add_item(input = {
    drc         <- 0 
    drq         <- 0
    vec <- list(c(0,5),c(0,5,6,7))
  })
  
  
  results <- run_sim(  
    npats=10,                              
    n_sim=1,                                  
    psa_bool = FALSE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,
    ipd = 1, input_out = c("drc","drq","vec"),
    constrained = TRUE
  )
  expect_equal(unname(results[[1]][[1]]$vec), list(list(c(0,5),c(0,5,6,7)),list(c(0,5),c(0,5,6,7))))
  expect_equal(unname(results[[1]][[1]]$extradata_raw[[1]]$vec), list(c(0,5),c(0,5,6,7)))
  expect_no_error(summary_results_det(results[[1]][[1]]))
  expect_no_error(summary_results_sim(results[[1]]))
  expect_no_error(summary_results_sens(results))
  
  
  
  
  # Many options set as TRUE ------------------------------------------------
  
  results <- run_sim(  
    npats=10,                              
    n_sim=1,                                  
    psa_bool = TRUE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,
    ipd = 1,
    input_out = c("drc","drq","vec"),
    timed_freq = TRUE, 
    accum_backwards = TRUE,
    continue_on_error = TRUE,
    seed = 2,
    constrained = TRUE
  )
  expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
  expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
  expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
  expect_no_error(summary_results_det(results[[1]][[1]]))
  expect_no_error(summary_results_sim(results[[1]]))
  expect_no_error(summary_results_sens(results))
  
  
  # Check continue on error works and inputs can change across simulations -------------------------------------------
  
  
  common_all_inputs <-add_item(input = {
    drc         <- 0 
    drq         <- ifelse(simulation == 3, Inf,0)
    err_flag <- if(simulation == 2){stop("ERROR")}else{NULL}
  })
  
  
  results <- run_sim(  
    npats=10,                              
    n_sim=3,                                  
    psa_bool = TRUE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,
    ipd = 2,
    input_out = c("drc","drq"),
    timed_freq = TRUE, 
    accum_backwards = TRUE,
    continue_on_error = TRUE,
    seed = 2,
    constrained = TRUE
  )
  
  expect_equal(unname(results[[1]][[1]]$total_lys), c(1,1))
  expect_equal(unname(results[[1]][[1]]$total_qalys), c(0,0))
  expect_equal(unname(results[[1]][[1]]$total_costs), c(0,0))
  expect_equal(unname(results[[1]][[1]]$timed_outputs$total_lys$aa), c(0,1))
  expect_equal(results[[1]][[2]], NULL)
  expect_equal(unname(results[[1]][[3]]$total_lys), c(0,0))
  expect_equal(unname(results[[1]][[3]]$total_qalys), c(0,0))
  expect_equal(unname(results[[1]][[3]]$total_costs), c(0,0))
  expect_no_error(summary_results_det(results[[1]][[1]]))
  expect_no_error(summary_results_sim(results[[1]]))
  expect_no_error(summary_results_sens(results))
  
  
  # several analyses --------------------------------------------------------
  
  
  common_all_inputs <-add_item(input = {
    drc         <- 0 
    drq         <- ifelse(simulation == 3, Inf,0)
    err_flag <- if(simulation == 2){stop("ERROR")}else{NULL}
  })
  
  
  results <- run_sim(  
    npats=10,                              
    n_sim=3,                                  
    psa_bool = TRUE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,
    ipd = 2,
    input_out = c("drc","drq"),
    timed_freq = TRUE, 
    accum_backwards = TRUE,
    continue_on_error = TRUE,
    seed = 2,
    n_sensitivity = 2,
    sensitivity_bool = TRUE,
    sensitivity_names = c("A","B"),
    constrained = TRUE
  )
  
  expect_equal(unname(results[[4]][[1]]$total_lys), c(1,1))
  expect_equal(unname(results[[4]][[1]]$total_qalys), c(0,0))
  expect_equal(unname(results[[4]][[1]]$total_costs), c(0,0))
  expect_equal(unname(results[[4]][[1]]$timed_outputs$total_lys$aa), c(0,1))
  expect_equal(results[[4]][[2]], NULL)
  expect_equal(unname(results[[4]][[3]]$total_lys), c(0,0))
  expect_equal(unname(results[[4]][[3]]$total_qalys), c(0,0))
  expect_equal(unname(results[[4]][[3]]$total_costs), c(0,0))
  expect_no_error(summary_results_det(results[[4]][[1]]))
  expect_no_error(summary_results_sim(results[[4]]))
  expect_no_error(summary_results_sens(results))
  
  
  
  # Over time outputs make sense --------------------------------------------
  common_all_inputs <-add_item(input = {
    drc         <- 0 
    drq         <- 0
    q_default <- 0.8
    
    cost_a <- 1
    cost_a_cycle_starttime <- 0
    cost_a_cycle_l <- 0.125
    cost_a_max_cycles <- 100
  })
  
  init_event_list <- 
    add_tte(arm=c("aa","bb"), evts = c("start","cycle","death") ,input={
      start <- 0
      death <- 10.1
      cycle <- 1
    })
  evt_react_list <-
    add_reactevt(name_evt = "start",
                 input = {
                   cost_insta <- 100
                   
                 }) %>%
    add_reactevt(name_evt = "cycle",
                 input = {
                   if(curtime<9){
                     new_event(list(cycle = curtime + 1.1))
                   }
                   cost_insta <- 200
                   
                 }) %>%
    add_reactevt(name_evt = "death",
                 input = {
                 }) 
  
  results <- run_sim(  
    npats=10,                              
    n_sim=1,                                  
    psa_bool = TRUE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,
    util_ongoing_list = 'q_default',
    cost_instant_list = "cost_insta",
    cost_cycle_list = "cost_a",
    ipd = 1,
    timed_freq = 0.25,
    constrained = TRUE
  )
  
  expect_equal(results[[1]][[1]]$timed_outputs$total_lys$aa, c(seq(0,10.1, by = 0.25),10.1))
  expect_equal(results[[1]][[1]]$timed_outputs$timepoints, c(seq(0,10.1, by = 0.25),10.1))
  expect_equal(results[[1]][[1]]$timed_outputs$cost_insta$aa,
               c(100, 100, 100, 100, 300, 300, 300, 300, 300, 500, 500, 500, 
                 500, 700, 700, 700, 700, 700, 900, 900, 900, 900, 1100, 1100, 
                 1100, 1100, 1300, 1300, 1300, 1300, 1300, 1500, 1500, 1500, 1500, 
                 1700, 1700, 1700, 1700, 1700, 1900, 1900))
  expect_equal(results[[1]][[1]]$timed_outputs$cost_a$aa,
               c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 
                 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 
                 64, 66, 68, 70, 72, 74, 76, 78, 80, 81))
  expect_equal(results[[1]][[1]]$cost_a[["aa"]],max(results[[1]][[1]]$timed_outputs$cost_a$aa))
  expect_equal(results[[1]][[1]]$total_lys[["aa"]],max(results[[1]][[1]]$timed_outputs$total_lys$aa))
  expect_equal(results[[1]][[1]]$total_qalys[["aa"]],max(results[[1]][[1]]$timed_outputs$total_qalys$aa))
  
  
  
  # Over time outputs make sense, backwards --------------------------------------------
  common_all_inputs <-add_item(input = {
    drc         <- 0 
    drq         <- 0
    q_default <- 0.8
    
    cost_a <- 0
    cost_a_cycle_starttime <- 0
    cost_a_cycle_l <- 0.125
    cost_a_max_cycles <- 100
  })
  
  init_event_list <- 
    add_tte(arm=c("aa","bb"), evts = c("start","cycle","death") ,input={
      start <- 0
      death <- 10.1
      cycle <- 1
    })
  evt_react_list <-
    add_reactevt(name_evt = "start",
                 input = {
                   cost_insta <- 100
                   
                 }) %>%
    add_reactevt(name_evt = "cycle",
                 input = {
                   if(curtime<9){
                     new_event(list(cycle = curtime + 1.1))
                   }
                   cost_insta <- 200
                   cost_a <- 1
                   
                 }) %>%
    add_reactevt(name_evt = "death",
                 input = {
                 }) 
  
  results <- run_sim(  
    npats=10,                              
    n_sim=1,                                  
    psa_bool = TRUE,                         
    arm_list = c("aa", "bb"),             
    common_all_inputs = common_all_inputs,    
    init_event_list = init_event_list,        
    evt_react_list = evt_react_list,
    util_ongoing_list = 'q_default',
    cost_instant_list = "cost_insta",
    cost_cycle_list = "cost_a",
    ipd = 1,
    timed_freq = 0.25, accum_backwards = TRUE,
    constrained = TRUE
  )
  
  expect_equal(results[[1]][[1]]$timed_outputs$total_lys$aa, c(seq(0,10.1, by = 0.25),10.1))
  expect_equal(results[[1]][[1]]$timed_outputs$timepoints, c(seq(0,10.1, by = 0.25),10.1))
  expect_equal(results[[1]][[1]]$timed_outputs$cost_insta$aa,
               c(100, 100, 100, 100, 300, 300, 300, 300, 300, 500, 500, 500, 
                 500, 700, 700, 700, 700, 700, 900, 900, 900, 900, 1100, 1100, 
                 1100, 1100, 1300, 1300, 1300, 1300, 1300, 1500, 1500, 1500, 1500, 
                 1700, 1700, 1700, 1700, 1700, 1900, 1900))
  expect_equal(results[[1]][[1]]$timed_outputs$cost_a$aa,
               c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 
                 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 
                 64, 66, 68, 70, 72, 74, 76, 78, 80, 81))
  expect_equal(results[[1]][[1]]$cost_a[["aa"]],max(results[[1]][[1]]$timed_outputs$cost_a$aa))
  expect_equal(results[[1]][[1]]$total_lys[["aa"]],max(results[[1]][[1]]$timed_outputs$total_lys$aa))
  expect_equal(results[[1]][[1]]$total_qalys[["aa"]],max(results[[1]][[1]]$timed_outputs$total_qalys$aa))
  
  

# shared resources and inputs ---------------------------------------------------

  common_all_inputs <-add_item(input = {
    util.sick   <- 0.8
    util.sicker <- 0.5
    cost.sick   <- 3000
    cost.sicker <- 7000
    cost.int    <- 1000
    coef_noint  <- log(0.2)
    HR_int      <- 0.8
    drc         <- 0.035 #different values than what's assumed by default
    drq         <- 0.035
    random_seed_sicker_i <- sample.int(100000,npats,replace = FALSE)
    beds <- resource_discrete(6)
    beds_free <- beds$n_free()
    shared_accumulator <- shared_input(0)
    value_accum <- shared_accumulator$value()
  })  
  
  
  #Put objects here that do not change as we loop through treatments for a patient
  common_pt_inputs <- add_item(death= max(0.0000001,rnorm(n=1, mean=12, sd=3))) 
  
  #Put objects here that change as we loop through treatments for each patient (e.g. events can affect fl.tx, but events do not affect nat.os.s)
  unique_pt_inputs <- add_item(fl.sick = 1,
                               q_default = util.sick,
                               c_default = cost.sick + if(arm=="int"){cost.int}else{0},
                               success_blocking_bed = FALSE,
                               had_to_queue = 0,
                               time_in_queue = NA,
                               time_start_queue = NA) 
  
  
  init_event_list <- 
    add_tte(arm=c("noint","int"), evts = c("sick","sicker","death") ,input={
      sick <- 0
      sicker <- draw_tte(1,dist="exp", coef1=coef_noint, beta_tx = ifelse(arm=="int",HR_int,1), seed = random_seed_sicker_i[i]) #this way the value would be the same if it wasn't for the HR, effectively "cloning" patients luck
      
    })
  
  evt_react_list <-
        add_reactevt(name_evt = "sick",
                     input = {
                       shared_accumulator <- shared_accumulator$modify(shared_accumulator$value() + 1)
                       value_accum <- shared_accumulator$value()
                       beds_free <- beds$n_free()
                       time_in_queue <- NA
                     }) %>%
        add_reactevt(name_evt = "sicker",
                     input = {
                       success_blocking_bed <- beds$attempt_block()
                       beds_free <- beds$n_free()
                       if(!success_blocking_bed){
                         time_start_queue <- curtime
                         modify_event(c(death = max(curtime,get_event("death") * 0.8)))
                         had_to_queue <- 1
                       }else{
                         time_in_queue <- ifelse(had_to_queue == 1, curtime - time_start_queue,NA)
                       }
                       q_default <- util.sicker
                       c_default <- cost.sicker + if(arm=="int"){cost.int}else{0}
                       fl.sick   <- 0 
                     }) %>%
        add_reactevt(name_evt = "death",
                     input = {
                       beds$attempt_free() #remove from using or from the queue
                       if(success_blocking_bed & beds$queue_size() > 0){
                         new_event(c(sicker = curtime),
                                   cur_evtlist,
                                   patient_id = beds$next_patient_in_line())
                       }
                       success_blocking_bed <- FALSE
                       time_in_queue <- NA
                       beds_free <- beds$n_free()
                       q_default <- 0
                       c_default <- 0
                       curtime   <- Inf
                     }) 
      
  util_ongoing <- "q_default"
  cost_ongoing <- "c_default"
        
  results <- run_sim(  
          npats=10,                               # number of patients to be simulated
          n_sim=1,                                  # number of simulations to run
          psa_bool = FALSE,                         # use PSA or not. If n_sim > 1 and psa_bool = FALSE, then difference in outcomes is due to sampling (number of pats simulated)  
          arm_list = c("int", "noint"),             # intervention list
          common_all_inputs = common_all_inputs,    # inputs common that do not change within a simulation
          common_pt_inputs = common_pt_inputs,      # inputs that change within a simulation but are not affected by the intervention
          unique_pt_inputs = unique_pt_inputs,
          init_event_list = init_event_list,        # initial event list
          evt_react_list = evt_react_list,          # reaction of events
          util_ongoing_list = util_ongoing,
          cost_ongoing_list = cost_ongoing,
          constrained = TRUE,
          ipd = 1,
          timed_freq = 0.25,
          input_out = c("beds_free","had_to_queue","time_in_queue","value_accum")
  )
  
      
  expect_true(max(results[[1]][[1]]$merged_df$time_in_queue,na.rm=TRUE)>0)
  expect_true(results[[1]][[1]]$value_accum[[1]]>1)
  expect_true(results[[1]][[1]]$had_to_queue[[1]]>0)
  expect_equal(results[[1]][[1]]$timed_outputs$beds_free[[1]],
               c(0, 6, 5.9, 5.9, 5.7, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 5.4, 
                 5.4, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 3.9, 3.9, 3.9, 
                 3.9, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3, 2.5, 
                 2.5, 2.5, 2.5, 2.8, 2.8, 2.8, 2.8, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 
                 2, 2, 2, 1.9, 1.9, 1.9, 1.9, 1.9, 2.2, 2.1),
               tolerance = 0.01)
  
})

})

    