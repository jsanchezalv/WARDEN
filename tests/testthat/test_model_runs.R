test_that("Test minimal model runs with some basic settings", {
  

# No discount rate, no qalys nor costs --------------------------------------------------------

  
  common_all_inputs <-add_item2(input = {
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
    
    common_all_inputs <-add_item2(input = {
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
    
    common_all_inputs <-add_item2(input = {
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
    
    common_all_inputs <-add_item2(input = {
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
    
    common_all_inputs <-add_item2(input = {
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

    common_all_inputs <-add_item2(input = {
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


    common_all_inputs <-add_item2(input = {
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
    
    
    common_all_inputs <-add_item2(input = {
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
    
})
    