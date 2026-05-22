test_that("Vector discounting equal to non-vector for single-length elements", {
  #Ongoing
  expect_equal(disc_ongoing_v(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500),
               suppressWarnings(disc_ongoing(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500)))
  
  #Instant
  expect_equal(disc_instant_v(lcldr=0.035, lclcurtime=3, lclval=2500),
               suppressWarnings(disc_instant(lcldr=0.035, lclcurtime=3, lclval=2500)))

 })


test_that("Discounting works with no discounting", {
  #Ongoing
  expect_equal(disc_ongoing_v(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500),
               2500)
  expect_equal(suppressWarnings(disc_ongoing(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500)),
               2500)
  
  #Instant
  expect_equal(disc_instant_v(lcldr=0, lclcurtime=2, lclval=2500),
               2500)
  expect_equal(suppressWarnings(disc_instant(lcldr=0, lclcurtime=2, lclval=2500)),
               2500)
  
  #Cycle
  expect_equal(disc_cycle_v(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500, cyclelength=1/12, starttime=0),
               12*2500)
  expect_equal(suppressWarnings(disc_cycle(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500, cyclelength=1/12, starttime=0)),
               12*2500)
})



test_that("Discounting works with odd numbers", {
  #Ongoing
  expect_equal(disc_ongoing_v(lcldr=0,lclprvtime=0, lclcurtime=Inf, lclval=2500),
               Inf)
  expect_equal(suppressWarnings(disc_ongoing(lcldr=0,lclprvtime=0, lclcurtime=Inf, lclval=2500)),
               Inf)
  expect_equal(disc_ongoing_v(lcldr=0.035,lclprvtime=0, lclcurtime=0, lclval=2500),
               0)
  expect_equal(suppressWarnings(disc_ongoing(lcldr=0.035,lclprvtime=0, lclcurtime=0, lclval=2500)),
               0)
  expect_equal(disc_ongoing_v(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=2500),
               0)
  expect_equal(suppressWarnings(disc_ongoing(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=2500)),
               0)
  #Inf*0 gives NaN, while the element-wise function just check wehterh prevtime and curtime are equal
  expect_equal(disc_ongoing_v(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=Inf),
               NaN)
  expect_equal(suppressWarnings(disc_ongoing(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=Inf)),
               NaN) 
  
  #Instant
  expect_equal(disc_instant_v(lcldr=0, lclcurtime=Inf, lclval=2500),
               2500)
  expect_equal(suppressWarnings(disc_instant(lcldr=0, lclcurtime=Inf, lclval=2500)),
               2500)
  expect_equal(disc_instant_v(lcldr=0.035, lclcurtime=0, lclval=2500),
               2500)
  expect_equal(suppressWarnings(disc_instant(lcldr=0.035, lclcurtime=0, lclval=2500)),
               2500)
  expect_equal(disc_instant_v(lcldr=0.035, lclcurtime=5, lclval=2500),
               2104.93292)
  expect_equal(suppressWarnings(disc_instant(lcldr=0.035, lclcurtime=5, lclval=2500)),
               2104.93292)
  #Inf*0 gives NaN
  expect_equal(disc_instant_v(lcldr=0.035, lclcurtime=5, lclval=Inf),
               Inf)
  expect_equal(suppressWarnings(disc_instant(lcldr=0.035, lclcurtime=5, lclval=Inf)),
               Inf) 
  expect_equal(disc_instant_v(lcldr=5, lclcurtime=5, lclval=2500),
               0.32150206)
  expect_equal(suppressWarnings(disc_instant(lcldr=5, lclcurtime=5, lclval=2500)),
               0.32150206) 
  
  #Cycle
  expect_equal(disc_cycle_v(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500, cyclelength=1/12, starttime=1),
               12*2500)
  expect_equal(disc_cycle_v(lcldr=0,lclprvtime=0, lclcurtime=0, lclval=2500, cyclelength=1/12, starttime=0),
               2500)
  expect_equal(disc_cycle_v(lcldr=0.035,lclprvtime=0, lclcurtime=0, lclval=2500, cyclelength=2, starttime=0),
               2500)
  # #Inf*0 gives NaN
  # expect_equal(disc_cycle_v(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=Inf, cyclelength=1/12, starttime=0),
  #              NaN)
})



test_that("Vectorial discounting working as expected with vectors", {
  #Ongoing
  expect_equal(disc_ongoing_v(lcldr=0.035,lclprvtime=c(0.5,0.5,0.5), lclcurtime=c(3,3,3), lclval=c(0,1000,Inf)),
               c(0,2354.66015,Inf))
  
  #Instant
  expect_equal(disc_instant_v(lcldr=0.035, lclcurtime=c(3,3,3), lclval=c(0,1000,Inf)),
               c(0,901.9427,Inf))
  
  #Cycle
  expect_equal(disc_cycle_v(lcldr=0.035, lclcurtime=c(3,3,3), lclval=c(0,1000,Inf),lclprvtime=c(0.5,0.5,0.5), cyclelength=c(1/12,1/12,1/12),starttime=c(0,0,0)),
               c(0,28296.443022132232727,Inf))  
  
  expect_equal(disc_cycle_v(0.035,0,1,1,1000,0,5),1000)
  expect_equal(disc_cycle_v(0.035,1,1,1,1000,0,5),966.1835748792273)
  expect_equal(disc_cycle_v(0.035,0,1,0,1000,0,5),1000)
  expect_equal(disc_cycle_v(0.035,0,1,2,1000,0,5),1966.1835748792273)
})

test_that("Create indicators works correctly",{
  expect_equal(
    create_indicators(
      2,
      10,
      c(1,1)
    ),
    c(0,1)
  )
  
  expect_equal(
    create_indicators(
      2,
      10,
      c(1,1),
      5
    ),
    c(0,0)
  )
  
  expect_equal(
    create_indicators(
      6,
      10,
      c(1,1),
      5
    ),
    c(1,0)
  )
  
  
  expect_equal(
    create_indicators(
      9,
      10,
      c(1,1),
      5
    ),
    c(0,0)
  )
  
  expect_error(
    create_indicators(
      12,
      10,
      c(1,1),
      5
    )
  )
  
  expect_error(
    create_indicators(
      9,
      10,
      rep(2,20),
      5
    )
  )
  
  expect_error(
    create_indicators(
      9,
      10,
      rep(2,3),
      20
    )
  )
})
  



test_that("Pick values vectorized work correctly",{
  expect_equal(
    pick_val_v(base = list(2,3,c(2, 3, 4)),
               psa =sapply(1:3,
                           function(x) eval(call(
                             c("rnorm","rnorm","rdirichlet")[[x]],
                             1,
                             c(2,3,list(c(2, 3, 4)))[[x]],
                             c(0.1,0.1,NULL)[[x]]
                             ))),
               sens = list(4,5,c(0.4,0.8,0.1)),
               psa_ind = FALSE,
               sens_ind = TRUE,
               indicator=list(1,2,c(3,4,5)),
               names_out=c("util","util2","dirichlet_vector") ,
               indicator_sens_binary = FALSE,
               sens_iterator = 5,
               distributions = list("rnorm","rnorm","rdirichlet"),
               covariances = list(0.1,0.1,NULL), deploy = FALSE ),
    list(util = 2, util2 = 3, dirichlet_vector = c(0.36, 0.54, 0.1))
    
  )
  
  expect_equal(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = FALSE,
      sens_ind = FALSE,
      indicator=c(1,0), deploy = FALSE
    ),
    list(0,0)
  )
  
  expect_equal(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = FALSE,
      sens_ind = TRUE,
      indicator=c(1,0), deploy = FALSE
    ),
    list(2,0)
  )
  
  expect_equal(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = FALSE,
      sens_ind = TRUE,
      indicator=c(0,1), deploy = FALSE
    ),
    list(0,3)
  )
  
  expect_equal(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = FALSE,
      sens_ind = TRUE,
      indicator=c(1,1), deploy = FALSE
    ),
    list(2,3)
  )
  
  expect_equal(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = TRUE,
      sens_ind = TRUE,
      indicator=c(1,1), deploy = FALSE
    ),
    list(2,3)
  )
  
  expect_error(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = TRUE,
      sens_ind = TRUE,
      indicator=c(1,5), deploy = FALSE
    )
  )
  
  expect_error(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = 5,
      sens_ind = TRUE,
      indicator=c(1,1), deploy = FALSE
    )
  )
  
  expect_error(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = TRUE,
      sens_ind = 3,
      indicator=c(1,1), deploy = FALSE
    )
  )
  
  expect_equal(
    pick_val_v(
      base = c(0,0),
      psa ={c(draw_tte(1,'norm',0,0.1,seed=1),draw_tte(1,'norm',0,0.1,seed=2))},
      sens = c(2,3),
      psa_ind = TRUE,
      sens_ind = TRUE,
      indicator=c(0,0), deploy = FALSE
    ),
    {list(draw_tte(1,'norm',0,0.1,seed=1),draw_tte(1,'norm',0,0.1,seed=2))}
  )
  
  expect_equal(
    pick_val_v(
      base = c(0,0),
      psa ={c(draw_tte(1,'norm',0,0.1,seed=1),draw_tte(1,'norm',0,0.1,seed=2))},
      sens = c(2,3),
      psa_ind = TRUE,
      sens_ind = TRUE,
      indicator=c(1,0), deploy = FALSE
    ),
    {list(2,draw_tte(1,'norm',0,0.1,seed=2))}
  )
  

})


test_that("Conditional Multivariate normal works as expected",{
  expect_equal(
    cond_mvn(mu = c(1, 2, 3),
                    Sigma = matrix(c(0.2, 0.05, 0.1, 
                                     0.05, 0.3, 0.05, 
                                     0.1, 0.05, 0.4), nrow = 3),
                    i = 1:2,
                    xi = c(1.2,2.3),
                    full_output = TRUE
                    )$mean,
    c(1.2, 2.3, 3.1217391)
  )
  
  expect_equal(
    cond_mvn(mu = c(1, 2, 3),
                    Sigma = matrix(c(0.2, 0.05, 0.1, 
                                     0.05, 0.3, 0.05, 
                                     0.1, 0.05, 0.4), nrow = 3),
                    i = 1:2,
                    xi = c(1.2,2.3),
                    full_output = TRUE
    )$covariance,
    structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0.347826086956522), dim = c(3L, 
                                                                    3L))
  )
  
  expect_error(
    cond_mvn(mu = c(1, 2, 3),
                    Sigma = matrix(c(0, 0, 0, 
                                     0, 0, 0, 
                                     0, 0, 0), nrow = 3),
                    i = 1:2,
                    xi = c(1.2,2.3),
                    full_output = TRUE
    )
  )
  
  expect_error(
    cond_mvn(mu = c(1, 2, 3),
                    Sigma = matrix(c(0, 0, 0, 
                                     0, 0, 0, 
                                     0, 0, 0), nrow = 3),
                    i = 5,
                    xi = c(1.2,2.3),
                    full_output = TRUE
    )
  )
  
})


test_that("Model Reactions Interactivity summary can be created",{
  expr <- substitute({
    
    a <- sum(5+7)
    
    ggplot()
    
    data.frame(x=1,b=2)
    
    list(b=5)
    
    a <- list(s=7)
    
    
    j <- 6
    if(TRUE){modify_event(list(j=5))}
    
    l <- 9
    
    afsa=ifelse(TRUE,"asda",NULL)
    

      o_exn = o_exn + 1
      
      a = NULL
      
      b = if(a){"CZ"}else{"AW"}
      
      rnd_prob_exn_sev = runif(1)
      
      exn_sev = rnd_prob_exn_sev <= p_sev
      
      o_exn_mod = o_exn_mod + if(exn_sev) { 0 } else { 1 }
      
      o_exn_sev = o_exn_sev + if(exn_sev) { 1 } else { 0 }
      
      o_rec_time_without_exn = (o_exn == 0) * 1
      
      o_rec_time_without_exn_sev = (o_exn_sev == 0) * 1
      
      o_c_exn = if(exn_sev) { c_sev } else { c_mod }
  
      o_other_c_exn_mod = if(exn_sev) { 0 } else { c_mod }
      
      o_other_c_exn_sev = if(exn_sev) { c_sev } else { 0 }
      
      o_qloss_exn = -if(exn_sev) { q_sev } else { q_mod }
      
      o_other_qloss_exn_mod = -if(exn_sev) { 0 } else { q_mod }
    
      o_other_qloss_exn_sev = -if(exn_sev) { q_sev } else { 0 }
      
      o_qloss_cg_exn = -if(exn_sev) { q_cg_sev } else { q_cg_mod }
    
      o_other_qloss_cg_exn_mod = -if(exn_sev) { 0 } else { q_cg_mod }
      
      o_other_qloss_cg_exn_sev = -if(exn_sev) { q_cg_sev } else { 0 }
      
      o_q = utility
      
      o_other_q_gold1 = if(gold == 1) { utility } else { 0 }
      
      o_other_q_gold2 = if(gold == 2) { utility } else { 0 }
      
      o_other_q_gold3 = if(gold == 3) { utility } else { 0 }
      
      o_other_q_gold4 = if(gold == 4) { utility } else { 0 }
      
      o_other_q_on_dup = if(on_dup) { utility } else { 0 }
    
      n_exn = n_exn + 1
      
      n_exn_mod = n_exn_mod + (1 - exn_sev)
      
      n_exn_sev = n_exn_sev + exn_sev
      
      u_adj_exn_lt = u_adj_exn_lt + if(exn_sev) { u_adj_sev_lt } else { u_adj_mod_lt }
      
      utility = u_gold - u_adj_exn_lt - u_mace_lt
      
      o_rec_utility = utility
      
      rnd_exn = runif(1)
    
    
    if(a==1){
      a=list(6+b)
      
      modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
    } else{
      modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
      if(a>6){
        a=8
      }
      
    }
    
    
    if (sel_resp_incl == 1 & on_dup == 1) {
      
      modify_event(list(e_response = curtime, z = 6))
      
    }
    
  })
  
  expect_length(ast_as_list(expr),43) 
  
  expect_type(ast_as_list(expr),"list")
  
  expect_equal(class(extract_elements_from_list(ast_as_list(expr))),"data.frame")
  
  expect_length(extract_elements_from_list(ast_as_list(expr)),4) #4 columns
  
  expect_equal(nrow(extract_elements_from_list(ast_as_list(expr))),43) #43 items/events changed, including assignments
  
  a <- add_reactevt(name_evt="example",
                    input={
                      a <- 5
                      w=5
                    })
  
  
  expect_equal(nrow(extract_from_reactions(a)),2) #2 items/events changed, including assignments
  
  expect_equal(extract_from_reactions(a),
               data.table(event = c("example","example"),
                          name = c("a","w"),
                          type = c("item","item"),
                          conditional_flag = c(FALSE,FALSE),
                          definition = c("5","5"))
               ) 
  
  
  
})


test_that("add_tte works as expected", {
  initial_data <- list()
  arm <- "control"
  evts <- c("start", "end")
  
  result <- add_tte(.data = initial_data, arm = arm, evts = evts, input = {
    start <- 0
    end <- 100
  })
  
  expect_true("control" %in% names(result))
  expect_equal(result$control$evts, evts)
})

test_that("modify_event modifies events correctly", {
  input_list_arm <- environment()
  input_list_arm$qaly_default_instant = 100
  input_list_arm$debug = FALSE
  input_list_arm$cur_evtlist = FALSE
  input_list_arm$accum_backwards = FALSE
  event_queue <- queue_create(c("ae","new_event","nat.death","nonexistent"))
  new_event(c(ae=5,nat.death=10),event_queue,1)

  # Modify an existing event
  modify_event(list(ae = 10), create_if_missing = TRUE, event_queue,1)
  expect_equal(get_event("ae",event_queue,1),10)
  
  # Create new event if not exists
  modify_event(list(new_event = 50), create_if_missing = TRUE, event_queue,1)
  expect_equal(get_event("new_event",event_queue,1), 50)
  
  # Ignore non-existent event
  expect_no_condition(modify_event(list(nonexistent = 20), create_if_missing = FALSE, event_queue,1))
  expect_error(get_event("nonexistent",event_queue,1))
})

test_that("new_event adds new events correctly", {
  input_list_arm <- environment()
  input_list_arm$qaly_default_instant = 100
  input_list_arm$debug = FALSE
  input_list_arm$cur_evtlist = FALSE
  input_list_arm$accum_backwards = FALSE
  input_list_arm$cur_evtlist = c()
  event_queue <- queue_create(c("ae","not_numeric"))
  new_event(c(ae=5),event_queue,1)
  
  expect_equal(get_event("ae",event_queue,1), 5)
  
  expect_error(new_event(c("not_numeric" = "five")))
})


test_that("replicate_profiles works correctly", {
  profiles <- data.frame(id = 1:10, age = rnorm(10, 60, 5))
  
  # Test replication with replacement
  set.seed(42)
  result <- replicate_profiles(profiles, replications = 20, replacement = TRUE)
  expect_equal(nrow(result), 20)
  expect_true(all(result$id %in% profiles$id))
  
  # Test replication without replacement
  set.seed(42)
  result_no_replacement <- replicate_profiles(profiles, replications = 10, replacement = FALSE)
  expect_equal(nrow(result_no_replacement), 10)
  expect_equal(sort(result_no_replacement$id), sort(profiles$id))
  
})


test_that("add_reactevt adds reactions correctly", {
  # Create an empty data list
  data_list <- list()
  
  # Add a reaction
  result <- add_reactevt(.data = data_list, name_evt = "start", input = { curtime <- Inf })
  expect_true("start" %in% names(result))

  # Test error handling for invalid event name
  expect_error(add_reactevt(name_evt = c("evt1", "evt2"), input = {}), 
               "name_evt argument in add_reactevt should be a single string with at least 2 characters")
})


# Luck_adj ----------------------------------------------------------------


test_that("luck_adj adjusts luck correctly", {
  # Test single values
  adj <- luck_adj(prevsurv = 0.8, cursurv = 0.6, luck = 0.9, condq = TRUE)
  expect_true(adj > 0 & adj < 1)
  
  # Test vectorized adjustment
  adj_vec <- luck_adj(prevsurv = c(0.8, 0), cursurv = c(0.6, 0.5), luck = c(0.9, 0.8), condq = TRUE)
  expect_equal(length(adj_vec), 2)
  expect_equal(adj_vec[2], 0.8)
  
  # Test conditional adjustment
  adj_cond <- luck_adj(prevsurv = 0.8, cursurv = 0.6, luck = 0.9, condq = FALSE)
  expect_true(adj_cond > 0 & adj_cond < 1)
})


# Random Streams ----------------------------------------------------------


test_that("random_stream initializes correctly", {
  random_stream_instance <- random_stream(10)
  
  expect_type(random_stream_instance, "environment")
  expect_equal(length(random_stream_instance$stream), 10, info = "The stream should initialize with 10 elements")
})

test_that("draw_n draws correct number of elements", {
  random_stream_instance <- random_stream(10)
  
  drawn_numbers <- random_stream_instance$draw_n(3)
  
  expect_equal(length(drawn_numbers), 3, info = "draw_n should draw 3 elements")
  expect_equal(length(random_stream_instance$stream), 7, info = "The stream should have 7 elements left after drawing 3")
})

test_that("draw_n with larger n regenerates stream", {
  random_stream_instance <- random_stream(5)
  
  expect_warning(drawn_numbers_large <- random_stream_instance$draw_n(8),
                 "Stream is smaller than the number of numbers drawn",
                 info = "Should warn when trying to draw more elements than available")
  
  expect_equal(length(drawn_numbers_large), 8, info = "After regenerating, draw_n should return 8 elements")
  expect_equal(length(random_stream_instance$stream), 0, info = "After drawing all elements, the stream should be empty")
})

test_that("generate_stream changes stream length", {
  random_stream_instance <- random_stream(5)
  
  random_stream_instance$generate_stream(20)
  
  expect_equal(length(random_stream_instance$stream), 20, info = "generate_stream should change stream to the correct size")
})


# Time covariate tte ------------------------------------------------------
library(flexsurv)

test_that("qtimecov returns numeric scalar within bounds", {
  param_fun_factory <- function(p0, p1, p2, p3) {
    function(t) p0 + p1 * t + p2 * t^2 + p3 * (floor(t) + 1)
  }
  
  rate_exp <- param_fun_factory(0.1, 0, 0, 0)
  set.seed(1)
  tte <- qtimecov(runif(1), a_fun = rate_exp, dist = "exp")
  expect_type(tte$tte, "double")
  expect_length(tte, 2)
  expect_gt(tte$tte, 0)
  expect_lt(tte$tte, 100)
})

test_that("qtimecov works for all supported distributions", {
  set.seed(1)
  param_fun_factory <- function(p0, p1, p2, p3) {
    function(.time) p0 + p1 * .time + p2 * .time^2 + p3 * (floor(.time) + 1)
  }
  
  # 1. Exponential
  rate_exp <- param_fun_factory(0.1, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = rate_exp, dist = "exp"))
  expect_equal(qtimecov(luck = 0.5,a_fun = rate_exp,dist = "exp", dt = 0.001
  )$tte,qexp(0.5,0.1), tolerance = 0.01)
  
  # 2. Gamma
  shape <- param_fun_factory(2, 0, 0, 0)
  rate <- param_fun_factory(0.2, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = shape, b_fun = rate, dist = "gamma"))
  expect_equal(qtimecov(luck = 0.5,a_fun = shape, b_fun = rate, dist = "gamma", dt = 0.001
  )$tte,qgamma(0.5,2,0.2), tolerance = 0.01)
  
  # 3. Lognormal
  meanlog <- param_fun_factory(log(10) - 0.5^2 / 2, 0, 0, 0)
  sdlog <- param_fun_factory(0.5, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = meanlog, b_fun = sdlog, dist = "lnorm"))
  expect_equal(qtimecov(0.5, a_fun = meanlog, b_fun = sdlog, dist = "lnorm",dt=0.01)$tte,
               qlnorm(0.5,log(10) - 0.5^2 / 2,0.5), tolerance = 0.01)
  
  # 4. Normal
  mean <- param_fun_factory(10, 0, 0, 0)
  sd <- param_fun_factory(2, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = mean, b_fun = sd, dist = "norm"))
  expect_equal(qtimecov(0.5, a_fun = mean, b_fun = sd, dist = "norm",dt=0.01)$tte,
               qnorm(0.5,10,2), tolerance = 0.01)
  
  # 5. Weibull
  shape <- param_fun_factory(2, 0, 0, 0)
  scale <- param_fun_factory(10, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = shape, b_fun = scale, dist = "weibull"))
  expect_equal(qtimecov(0.5, a_fun = shape, b_fun = scale, dist = "weibull",dt=0.01)$tte,
               qweibull(0.5,2,10), tolerance = 0.01)
  
  # 5.1 WeibullPH
  shape <- param_fun_factory(2, 0, 0, 0)
  scale <- param_fun_factory(0.01, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = shape, b_fun = scale, dist = "weibullPH"))
  expect_equal(qtimecov(0.5, a_fun = shape, b_fun = scale, dist = "weibullPH",dt=0.01)$tte,
               flexsurv::qweibullPH(0.5,2,0.01), tolerance = 0.01)
  
  # 6. Loglogistic
  shape <- param_fun_factory(2.5, 0, 0, 0)
  scale <- param_fun_factory(7.6, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = shape, b_fun = scale, dist = "llogis"))
  expect_equal(qtimecov(0.5, a_fun = shape, b_fun = scale, dist = "llogis",dt=0.01)$tte,
               flexsurv::qllogis(0.5,2.5,7.6), tolerance = 0.01)
  
  # 7. Gompertz
  shape <- param_fun_factory(0.01, 0, 0, 0)
  rate <- param_fun_factory(0.091, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = shape, b_fun = rate, dist = "gompertz"))
  expect_equal(qtimecov(0.5, a_fun = shape, b_fun = rate, dist = "gompertz",dt=0.01)$tte,
               qgompertz(0.5,0.01,0.091), tolerance = 0.01)
  
  rate_exp <- function(t) 0.1
  init_luck <- 0.95
  expect_equal(qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001)$tte,{
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001,max_time = 10)
    
    qtimecov(luck = a$luck,a_fun = rate_exp,dist = "exp", dt = 0.001, start_time=a$tte)$tte},
    tolerance = 0.01)
  
  init_luck <- 0.99
  expect_equal(qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001)$tte,{
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001,max_time = 5)
    
    qtimecov(luck = a$luck,a_fun = rate_exp,dist = "exp", dt = 0.001, start_time=a$tte)$tte},
    tolerance = 0.01)
  
  
  init_luck <- 0.3
  expect_equal(qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001)$tte,{
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001,max_time = 1)
    
    qtimecov(luck = a$luck,a_fun = rate_exp,dist = "exp", dt = 0.001, start_time=a$tte)$tte},
    tolerance = 0.01)
  
  
  
  init_luck <- 0.3
  expect_equal(qtimecov(luck = init_luck,a_fun = rate_exp,b_fun = function(t) 90000 ,dist = "weibull", dt = 0.001)$tte,{
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,b_fun = function(t) 90000,dist = "weibull", dt = 0.001,max_time = 1)
    
    qtimecov(luck = a$luck,a_fun = rate_exp,b_fun = function(t) 90000,dist = "weibull", dt = 0.001, start_time=a$tte)$tte},
    tolerance = 0.01)
  
  rate_exp <- function(.time) 0.1
  rate_exp2 <- function(.time) 0.2
  time_change <- 10
  init_luck <- 0.95
  
  expect_equal({
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005, max_time = time_change)
    qtimecov(luck = a$luck,a_fun = rate_exp2,dist = "exp", dt = 0.005, start_time=a$tte)$tte
    
  },{
    new_luck <- luck_adj(prevsurv = 1 - pexp(q=time_change,rate_exp(1)),
                         cursurv = 1 - pexp(q=time_change,rate_exp2(1)),
                         luck = init_luck,
                         condq = FALSE) #time 10 change
    qexp(new_luck,rate_exp2(1))
  }, tolerance = 0.01)
  
  

# time varying and an event -----------------------------------------------

  rate_exp <- function(.time) 0.1 + 0.01*.time * 0.00001*.time^2
  rate_exp2 <- function(.time) 0.2 + 0.02*.time
  time_change <- 8
  init_luck <- 0.95
  
  expect_equal({
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005, max_time = time_change)
    qtimecov(luck = a$luck,a_fun = rate_exp2,dist = "exp", dt = 0.005, start_time=a$tte)$tte
    
  },{# Manually reproduce what qtimecov does from a$tte
    t <- 0
    luck <- init_luck
    dt <- 0.005
    repeat {
      t <- t + dt
      s_prev <- 1 - pexp(t - dt, rate = rate_exp(t - dt))
      s_curr <- 1 - pexp(t,       rate = rate_exp(t - dt))
      
      luck <- luck_adj(prevsurv = s_prev, cursurv = s_curr, luck = luck, condq = TRUE)
      
      res_tte <- qcond_exp(luck, rate = rate_exp(t))
      total_tte <- t + res_tte
      
      if (res_tte <= dt || total_tte <= t || t >= time_change) {
        break
      }
    }
    
    if (total_tte <= time_change) {
      return(total_tte)
    }
    
    # Phase 2: after change
    repeat {
      t <- t + dt
      s_prev <- 1 - pexp(t - dt, rate = rate_exp2(t - dt))
      s_curr <- 1 - pexp(t,       rate = rate_exp2(t - dt))
      
      luck <- luck_adj(prevsurv = s_prev, cursurv = s_curr, luck = luck, condq = TRUE)
      
      res_tte <- qcond_exp(luck, rate = rate_exp2(t))
      total_tte <- t  + res_tte
      
      if (res_tte <= dt || total_tte <= t) {
        break
      }
    }
    
    total_tte
  }, tolerance = 0.01)
  
  rate_exp <- function(.time) 0.1 + floor(.time)*0.01
  init_luck <- 0.95
  
  #30x faster 
  expect_equal({
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 1)$tte
    a
  },{
    time_vec <- 0:100
    new_luck <- init_luck
    for (i in 2:length(time_vec)) {
      new_luck <- luck_adj(prevsurv = 1 - pexp(q=time_vec[i-1],rate_exp(time_vec[i-1])),
                           cursurv = 1 - pexp(q=time_vec[i],rate_exp(time_vec[i])),
                           luck = new_luck,
                           condq = FALSE) 
      qexp(new_luck,rate_exp(1))
    }
    a
  }, tolerance = 0.01)
  
  
  expect_equal({
    bs_rate <- 0.08
    time_cov <- 0.0002271
    bs_rate_f <- function(.time) bs_rate + time_cov * floor(.time)
    luck <- 0.403
    t <- 0
    
    # First rate and quantile (from time 0)
    new_rate <- bs_rate_f(0)
    ttdeath <- qexp(luck, rate = new_rate)
    
    # Unconditional quantile loop
    while (t < ttdeath) {
      old_rate <- new_rate
      new_rate <- bs_rate_f(t)
      
      # Both survival values are evaluated at time t (NOT t+1)
      prev_surv <- 1 - pexp(q = t, rate = old_rate)
      cur_surv  <- 1 - pexp(q = t, rate = new_rate)
      
      # Use condq = FALSE for unconditional approach
      luck <- luck_adj(prevsurv = prev_surv, cursurv = cur_surv, luck = luck, condq = FALSE)
      
      # Unconditional quantile from time 0
      ttdeath <- qexp(luck, rate = new_rate)
      
      t <- t + 1
    }
    
    # Final time-to-event
    ttdeath
    
  },{
    new_rate <-bs_rate <-  0.08
    bs_rate_f <- function(.time){bs_rate + (time_cov * floor(.time))}
    luck <- 0.403
    time_cov <- 0.0002271
    t <- 0
    qtimecov(luck = luck,a_fun = bs_rate_f,dist = "exp", dt = 1)$tte
  }, tolerance = 0.01)
  
})



test_that("qtimecov throws error for unsupported distribution", {
  dummy <- function(.time) 1
  expect_error(qtimecov(runif(1), a_fun = dummy, dist = "beta"), "Unsupported distribution")
})

test_that("qtimecov respects max_time bound", {
  slow_fun <- function(.time) 0.00001
  tte <- qtimecov(
    luck = 0.999,
    a_fun = slow_fun,
    dist = "exp",
    max_time = 10
  )$tte
  expect_lte(tte, 10)
})


test_that("adj_val works as intended",{
  
  bs_age <- 1
  vec <- 1 - seq(from =0, to = 0.2, by = 0.02)
  expect_equal(adj_val(0,0,by=1, vec[floor(.time + bs_age)], discount = Inf), 0)
  expect_equal(adj_val(0,5,by=1, vec[floor(.time + bs_age)], discount = Inf), 0)
  expect_error(adj_val(0,20,by=1, vec[floor(.time + bs_age)]))
  expect_error(adj_val(0,20,by=1, vec[floor(.time + bs_age)]))
  expect_error(adj_val(0,-5,by=1, vec[floor(.time + bs_age)]))
  expect_equal(adj_val(0,0,by=1, vec[floor(.time + bs_age)]), 0)
  expect_equal(adj_val(0,0.1,by=1, vec[floor(.time + bs_age)], discount = 0),1)
  expect_equal(adj_val(0,1.1,by=1, vec[floor(.time + bs_age)] * .time, discount = 0),(1*0+0.98*0.1)/1.1)
  expect_equal(adj_val(0,0.1,by=1, vec[floor(.time + bs_age)] * .time, discount = 0),0)
  expect_equal(adj_val(0,0.1,by=1, vec[floor(.time + bs_age)] * .time, discount = 0.03),0)
  expect_equal(adj_val(8,9,by=0.2, vec[floor(.time + bs_age)], discount = 0),0.84)
  expect_equal(adj_val(8,9,by=0.5, vec[floor(.time + bs_age)], discount = Inf),0)
  expect_equal(adj_val(8,9,by=0.5, 1),1)
  expect_equal(adj_val(0,4,by=1, .time),1.5)
  expect_equal(adj_val(0,4,by=1, .time),adj_val(0,4,by=1, .time,vectorized_f = TRUE))
  
  expect_equal({
    val <- 0.8 * adj_val(0,5.2,by=1, vec[floor(.time + bs_age)], discount = 0.03)
    disc_ongoing_v(0.03,
                   0,
                   5.2,
                   val)
  },{
    a <- 0
    time_vec <- c(0,1,2,3,4,5,5.2)
    for (i in 1:(length(time_vec)-1)) {
      a <-  a + disc_ongoing_v(0.03,
                               time_vec[i],
                               time_vec[i+1],
                               0.8 * vec[floor(time_vec[i] + bs_age)])
    }
    a
  })
  
  #17x faster
  expect_equal({
    val <- 0.8 * adj_val(0,5.2,by=1, vec[floor(.time + bs_age)], discount = 0)
    disc_ongoing_v(0,
                   0,
                   5.2,
                   val)
  },{
    a <- 0
    time_vec <- c(0,1,2,3,4,5,5.2)
    for (i in 1:(length(time_vec)-1)) {
      a <-  a + disc_ongoing_v(0,
                               time_vec[i],
                               time_vec[i+1],
                               0.8 * vec[floor(time_vec[i] + bs_age)])
    }
    a
  })
  
})

test_that("qcondweibull and qcondweibullPH can give same results",
          {
            # Set base R parameters
            shape <- 2
            scale_base <- 10
            # Convert to PH parameter
            scale_PH <- (1 / scale_base)^shape  # = 0.01
            
            # Set common inputs
            p <- 0.5
            t0 <- 5
            
            # Evaluate both
            expect_equal(qcond_weibull(p, shape, scale_base, lower_bound = t0),
            qcond_weibullPH(p, shape, scale_PH, lower_bound = t0))
            }
          )




# shared inputs -----------------------------------------------------------

test_that("immutable mode: modify returns new object; originals unaffected", {
  a <- shared_input(5, constrained = FALSE)
  expect_s3_class(a, c("shared_input_val", "shared_input"))
  
  a2 <- a$modify(a$value() + 7)
  expect_equal(a$value(), 5)
  expect_equal(a2$value(), 12)
  
  # clone is a deep copy of the wrapper's current value
  a3 <- a2$clone()
  a2b <- a2$modify(99)
  expect_equal(a2$value(), 12)     # old handle unchanged
  expect_equal(a2b$value(), 99)
  expect_equal(a3$value(), 12)     # clone still has 12
})

test_that("immutable mode: reset and fork work", {
  a <- shared_input(10, constrained = FALSE)
  
  # reset goes back to initial and returns a fresh handle
  a2 <- a$modify(11)
  a3 <- a2$reset()
  expect_equal(a2$value(), 11) # old handle unchanged
  expect_equal(a3$value(), 10) # new handle reset to init
  
  # fork returns n independent deep clones (of the wrapper)
  forks <- a$fork(3)
  expect_type(forks, "list")
  expect_length(forks, 3)
  expect_true(all(vapply(forks, function(x) x$value(), numeric(1)) == 10))
  
  forks[[1]] <- forks[[1]]$modify(77)
  expect_equal(forks[[1]]$value(), 77)
  expect_equal(forks[[2]]$value(), 10) # independence across forks
})

test_that("shared mode: aliasing across handles and clone/reset break sharing", {
  make_shared <- function(val) {
    constrained <- TRUE  # parent flag consumed by shared_input()
    shared_input(val)
  }
  
  b1 <- make_shared(10)
  expect_s3_class(b1, c("shared_input_env", "shared_input"))
  
  # alias: same underlying state
  b2 <- b1
  expect_equal(b1$value(), 10)
  expect_equal(b2$value(), 10)
  
  # modify returns a fresh wrapper that *shares* the same state
  b1 <- b1$modify(b1$value() + 1)
  expect_equal(b1$value(), 11)
  expect_equal(b2$value(), 11)  # reflects shared state
  
  # clone: new independent state
  b3 <- b1$clone()
  b1 <- b1$modify(99)
  expect_equal(b1$value(), 99)
  expect_equal(b2$value(), 99)  # still sharing with b1
  expect_equal(b3$value(), 11)  # independent
  
  # reset: new independent state set to init
  b4 <- b1$reset()
  expect_equal(b4$value(), 10)  # initial value
  b1 <- b1$modify(123)
  expect_equal(b4$value(), 10)  # independent from b1/b2
  
  b5 <- b3$modify(b3$value() + 1)
  
  expect_equal(b5$value(), b3$value())  # same value
})

test_that("explicit constrained arg overrides parent lookup", {
  constrained <- TRUE
  x <- shared_input(1, constrained = FALSE)
  expect_s3_class(x, "shared_input_val")
  y <- shared_input(1, constrained = TRUE)
  expect_s3_class(y, "shared_input_env")
})

test_that("parent lookup defaults to immutable unless TRUE", {
  # no 'constrained' in parent -> immutable
  x <- (function() shared_input(1))()
  expect_s3_class(x, "shared_input_val")
  
  # constrained present but FALSE -> immutable
  x <- (function() { constrained <- FALSE; shared_input(1) })()
  expect_s3_class(x, "shared_input_val")
  
  # constrained present and TRUE -> shared
  x <- (function() { constrained <- TRUE; shared_input(1) })()
  expect_s3_class(x, "shared_input_env")
})

test_that("works with non-scalar values; note on reference types", {
  # list value copies fine (but elements can be reference types)
  a <- shared_input(list(x = 1L), constrained = FALSE)
  a2 <- a$modify(list(x = 2L))
  expect_equal(a$value(), list(x = 1L))
  expect_equal(a2$value(), list(x = 2L))
  
  # environment as value: wrapper deep-copies the handle, but *not* the env itself
  e <- new.env(parent = emptyenv()); e$k <- 1L
  a <- shared_input(e, constrained = FALSE)
  a2 <- a$clone()
  e$k <- 99L
  # both see 99 because the underlying env is shared by R semantics
  expect_equal(a$value()$k, 99L)
  expect_equal(a2$value()$k, 99L)
})



# add_item ---------------------------------------------------------------

test_that("add_item: standalone with input= returns correct block", {
  result <- add_item(input = { a <- 1; b <- 2 })
  expect_equal(result, quote({ a <- 1; b <- 2 }))
})

test_that("add_item: named ... args become assignments", {
  result <- add_item(x = 5, y = 10)
  expect_equal(result, quote({ x <- 5; y <- 10 }))
})

test_that("add_item: native pipe chains accumulate expressions", {
  result <- add_item(input = { a <- 1 }) |>
    add_item(input = { b <- 2 }) |>
    add_item(c = 3)
  expect_equal(result, quote({ a <- 1; b <- 2; c <- 3 }))
})


test_that("add_item: .data arg directly accepts prior block", {
  prior <- add_item(input = { x <- 10 })
  result <- add_item(.data = prior, y = 20)
  expect_equal(result, quote({ x <- 10; y <- 20 }))
})

# add_item: unnamed expression as first positional arg ----------------------------

test_that("add_item: unnamed function call is captured unevaluated", {
  # sensitivity_bool does not exist in this scope — must not be evaluated
  result <- add_item(some_fn(x = sensitivity_bool))
  expect_equal(result, quote({ some_fn(x = sensitivity_bool) }))
})

test_that("add_item: unnamed if() expression is captured unevaluated", {
  result <- add_item(if(x_bool) create_x() else rep(1, 5))
  expect_equal(result, quote({ if(x_bool) create_x() else rep(1, 5) }))
})

test_that("add_item: unnamed call followed by named args", {
  result <- add_item(some_fn(a = b), x = 5)
  expect_equal(result, quote({ some_fn(a = b); x <- 5 }))
})

test_that("add_item: native pipe after unnamed call accumulates correctly", {
  result <- add_item(some_fn(a = b)) |> add_item(x = 5)
  expect_equal(result, quote({ some_fn(a = b); x <- 5 }))
})

test_that("add_item: pick_val_v-style call with absent engine vars not evaluated", {
  # psa_bool, sensitivity_bool, indicators are not in scope — must not error
  result <- add_item(
    pick_val_v(base = list(1), psa = list(2), sens = list(3),
               psa_ind = psa_bool, sens_ind = sensitivity_bool,
               indicator = indicators, names_out = list("a"))
  ) |> add_item(x = 5)
  stmts <- as.list(result)[-1L]
  expect_length(stmts, 2L)
  expect_identical(stmts[[1L]][[1L]], as.name("pick_val_v"))
  expect_equal(stmts[[2L]], quote(x <- 5))
})

test_that("add_item: inline brace block as first positional becomes starting block", {
  result <- add_item({ a <- 1; b <- 2 }, c = 3)
  expect_equal(result, quote({ a <- 1; b <- 2; c <- 3 }))
})

test_that("add_item: variable holding prior block used positionally as start", {
  prior <- add_item(x = 10)
  result <- add_item(prior, y = 20)
  expect_equal(result, quote({ x <- 10; y <- 20 }))
})

test_that("add_item: unnamed call among named args is captured unevaluated", {
  # R matches the first unnamed positional to .data, so it leads the block;
  # named args follow in their call order
  result <- add_item(x = 1, some_fn(a = absent_var), y = 2)
  expect_equal(result, quote({ some_fn(a = absent_var); x <- 1; y <- 2 }))
})

test_that("add_item: three-item native pipe chain with unnamed first call", {
  result <- add_item(some_fn(a = b)) |> add_item(x = 5) |> add_item(y = 6)
  expect_equal(result, quote({ some_fn(a = b); x <- 5; y <- 6 }))
})


# pick_val_v grouped mode PSA bug fix --------------------------------------

test_that("pick_val_v grouped mode: indicator_psa=NULL uses PSA for all when psa_ind=TRUE", {
  set.seed(1)
  psa_vals <- list(rnorm(1, 2, 0.1), rnorm(1, 3, 0.1), rnorm(1, 4, 0.1))
  result <- pick_val_v(
    base = list(2, 3, 4),
    psa  = psa_vals,
    sens = list(10, 11, 12),
    psa_ind  = TRUE,
    sens_ind = TRUE,
    indicator = list(1, 2, 3),
    indicator_psa = NULL,
    names_out = c("a", "b", "c"),
    indicator_sens_binary = FALSE,
    sens_iterator = 99,
    distributions = list("rnorm", "rnorm", "rnorm"),
    covariances = list(0.1, 0.1, 0.1),
    deploy_env = FALSE
  )
  expect_equal(result, setNames(psa_vals, c("a", "b", "c")))
})

test_that("pick_val_v grouped mode: indicator_psa respects 0 entries", {
  set.seed(1)
  psa_vals <- list(rnorm(1, 2, 0.1), rnorm(1, 3, 0.1), rnorm(1, 4, 0.1))
  result <- pick_val_v(
    base = list(2, 3, 4),
    psa  = psa_vals,
    sens = list(10, 11, 12),
    psa_ind  = TRUE,
    sens_ind = TRUE,
    indicator = list(1, 2, 3),
    indicator_psa = list(1, 0, 1),
    names_out = c("a", "b", "c"),
    indicator_sens_binary = FALSE,
    sens_iterator = 99,
    distributions = list("rnorm", "rnorm", "rnorm"),
    covariances = list(0.1, 0.1, 0.1),
    deploy_env = FALSE
  )
  expect_equal(result$a, psa_vals[[1]])
  expect_equal(result$b, 3)
  expect_equal(result$c, psa_vals[[3]])
})

test_that("pick_val_v grouped mode: all indicator_psa=0 uses base when psa_ind=TRUE", {
  set.seed(1)
  psa_vals <- list(rnorm(1, 2, 0.1), rnorm(1, 3, 0.1))
  result <- pick_val_v(
    base = list(2, 3),
    psa  = psa_vals,
    sens = list(10, 11),
    psa_ind  = TRUE,
    sens_ind = TRUE,
    indicator = list(1, 2),
    indicator_psa = list(0, 0),
    names_out = c("a", "b"),
    indicator_sens_binary = FALSE,
    sens_iterator = 99,
    distributions = list("rnorm", "rnorm"),
    covariances = list(0.1, 0.1),
    deploy_env = FALSE
  )
  expect_equal(result$a, 2)
  expect_equal(result$b, 3)
})

test_that("pick_val_v grouped mode: vector indicator_psa applies element-wise", {
  set.seed(1)
  psa_v <- c(rnorm(1, 0.4, 0.05), rnorm(1, 0.3, 0.05), rnorm(1, 0.3, 0.05))
  result <- pick_val_v(
    base = list(2, c(0.4, 0.3, 0.3)),
    psa  = list(rnorm(1, 2, 0.1), psa_v),
    sens = list(10, c(0.6, 0.5, 0.5)),
    psa_ind  = TRUE,
    sens_ind = TRUE,
    indicator = list(1, c(2, 3, 4)),
    indicator_psa = list(0, c(1, 0, 1)),
    names_out = c("scalar", "vector"),
    indicator_sens_binary = FALSE,
    sens_iterator = 99,
    distributions = list("rnorm", "rdirichlet"),
    covariances = list(0.1, NULL),
    deploy_env = FALSE
  )
  expect_equal(result$scalar, 2)
  expect_equal(result$vector, c(psa_v[1], 0.3, psa_v[3]))
})

test_that("pick_val_v grouped mode: sens_iterator match updates value correctly", {
  result <- pick_val_v(
    base = list(2, 3, 4),
    psa  = list(2, 3, 4),
    sens = list(10, 11, 12),
    psa_ind  = FALSE,
    sens_ind = TRUE,
    indicator = list(1, 2, 3),
    indicator_psa = list(1, 1, 1),
    names_out = c("a", "b", "c"),
    indicator_sens_binary = FALSE,
    sens_iterator = 2,
    distributions = list("rnorm", "rnorm", "rnorm"),
    covariances = list(0.1, 0.1, 0.1),
    deploy_env = FALSE
  )
  expect_equal(result$a, 2)
  expect_equal(result$b, 11)
  expect_equal(result$c, 4)
})

test_that("pick_val_v binary mode: indicator_psa behavior unchanged by grouped fix", {
  set.seed(1)
  psa1 <- rnorm(1, 0, 0.1)
  psa2 <- rnorm(1, 0, 0.1)
  result <- pick_val_v(
    base = c(0, 0),
    psa  = c(psa1, psa2),
    sens = c(2, 3),
    psa_ind  = TRUE,
    sens_ind = FALSE,
    indicator = c(0, 0),
    indicator_psa = c(1, 0),
    deploy_env = FALSE
  )
  expect_equal(result[[1]], psa1)
  expect_equal(result[[2]], 0)
})


# input_block() ------------------------------------------------------------

test_that("input_block: returns a { } call", {
  blk <- input_block(
    base      = list(1, 2),
    psa       = list(rnorm(1), rnorm(1)),
    sens      = my_sens,
    names_out = c("a", "b")
  )
  expect_true(is.call(blk))
  expect_identical(blk[[1L]], as.name("{"))
})

test_that("input_block: has warden_block_meta attribute with correct structure", {
  blk <- input_block(
    base      = list(1, 2, 3),
    psa       = list(1, 2, 3),
    sens      = my_sens,
    names_out = c("a", "b", "c")
  )
  meta <- attr(blk, "warden_block_meta")
  expect_type(meta, "list")
  expect_equal(meta$mode, "binary")
  expect_equal(meta$n_params, 3L)
  expect_equal(meta$n_groups, 3L)
})

test_that("input_block: grouped mode sets correct metadata", {
  blk <- input_block(
    base            = list(1, 2, 3, 4),
    psa             = list(1, 2, 3, 4),
    sens            = my_sens,
    names_out       = c("a", "b", "c", "d"),
    sens_indicators = list(1L, 1L, 2L, 2L),
    distributions   = list("rnorm", "rnorm", "rnorm", "rnorm"),
    covariances     = list(0.1, 0.1, 0.1, 0.1)
  )
  meta <- attr(blk, "warden_block_meta")
  expect_equal(meta$mode, "grouped")
  expect_equal(meta$n_params, 4L)
  expect_equal(meta$n_groups, 2L)
})

test_that("input_block: generated expression contains pick_val_v call", {
  blk <- input_block(
    base      = list(1, 2),
    psa       = pick_psa(list("rnorm"), list(1), list(0), list(0.1)),
    sens      = my_sens,
    names_out = c("a", "b")
  )
  stmts <- as.list(blk)[-1L]
  expect_length(stmts, 1L)
  expect_identical(stmts[[1L]][[1L]], as.name("pick_val_v"))
})

test_that("input_block: binary DSA mode uses create_indicators with 0 offset", {
  blk <- input_block(
    base      = list(1, 2),
    psa       = list(1, 2),
    sens      = my_sens,
    names_out = c("a", "b"),
    indicator_sens_binary = TRUE,
    dsa_names = c("dsa_min")
  )
  # Structure: { if(is_dsa_cond, pvv_dsa, pvv_scenario) }
  pvv_call <- as.list(blk)[[2L]]
  expect_identical(pvv_call[[1L]], as.name("if"))
  pvv_dsa  <- pvv_call[[3L]]
  expect_identical(pvv_dsa[[1L]], as.name("pick_val_v"))
  ci_call  <- pvv_dsa$indicator
  expect_identical(ci_call[[1L]], as.name("create_indicators"))
  expect_equal(ci_call[[5L]], 0L)
})

test_that("input_block: sens is indexed by sens_name_used in the expression", {
  blk <- input_block(
    base      = list(1),
    psa       = list(1),
    sens      = l_inputs,
    names_out = "a"
  )
  pvv   <- as.list(blk)[[2L]]
  s_arg <- pvv$sens
  expect_identical(s_arg[[1L]], as.name("[["))
  expect_identical(s_arg[[2L]], as.name("l_inputs"))
  expect_identical(s_arg[[3L]], as.name("sens_name_used"))
})

test_that("input_block: pipe chaining with add_item preserves block content", {
  blk <- input_block(
    base      = list(1, 2),
    psa       = list(1, 2),
    sens      = my_sens,
    names_out = c("a", "b")
  ) |> add_item(extra = 5)
  stmts <- as.list(blk)[-1L]
  expect_length(stmts, 2L)
  expect_identical(stmts[[1L]][[1L]], as.name("pick_val_v"))
  expect_equal(stmts[[2L]], quote(extra <- 5))
})

test_that("input_block: pipe chaining preserves warden_block_meta", {
  blk <- input_block(
    base      = list(1, 2),
    psa       = list(1, 2),
    sens      = my_sens,
    names_out = c("a", "b")
  ) |> add_item(extra = 5)
  expect_false(is.null(attr(blk, "warden_block_meta")))
  expect_equal(attr(blk, "warden_block_meta")$n_params, 2L)
})

test_that("input_block: errors when names_out is empty", {
  expect_error(
    input_block(base = list(), psa = list(), sens = my_sens, names_out = character(0)),
    "names_out must have at least one element"
  )
})

test_that("input_block: errors in grouped mode when distributions is NULL", {
  expect_error(
    input_block(
      base            = list(1, 2),
      psa             = list(1, 2),
      sens            = my_sens,
      names_out       = c("a", "b"),
      sens_indicators = list(1L, 2L)
    ),
    "distributions cannot be NULL in grouped mode"
  )
})

test_that("input_block: .data prepended correctly", {
  prior <- add_item(x = 10)
  blk   <- input_block(
    .data     = prior,
    base      = list(1),
    psa       = list(1),
    sens      = my_sens,
    names_out = "a"
  )
  stmts <- as.list(blk)[-1L]
  expect_length(stmts, 2L)
  expect_equal(stmts[[1L]], quote(x <- 10))
  expect_identical(stmts[[2L]][[1L]], as.name("pick_val_v"))
})


# input_block() integration with run_sim() ---------------------------------

test_that("input_block: run_sim produces valid output", {
  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0
      end   <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  i_block <- input_block(
    base      = list(0.8, 0.5),
    psa       = list(0.8, 0.5),
    sens      = list(),
    names_out = c("util.sick", "util.sicker")
  )

  results <- run_sim(
    npats             = 5,
    n_sim             = 1,
    psa_bool          = FALSE,
    sensitivity_bool  = FALSE,
    arm_list          = c("arm1"),
    common_all_inputs = i_block,
    init_event_list   = make_evt_list(),
    evt_react_list    = make_react_list(),
    util_ongoing_list = "util.sick",
    seed              = 42
  )

  expect_type(results, "list")
  expect_length(results, 1L)
  expect_true(results[[1]][[1]]$total_qalys > 0)
})

test_that("input_block: run_sim auto-computes n_sensitivity from metadata", {
  assign("l_inputs_test", list(
    parameter_name = list("util.sick", "util.sicker"),
    base_value     = list(0.8, 0.5),
    PSA_dist       = list("rnorm", "rbeta_mse"),
    a              = list(0.8, 0.5),
    b              = list(0.04, 0.025),
    n              = list(1, 1),
    DSA_min        = list(0.6, 0.3),
    DSA_max        = list(0.9, 0.7),
    psa_indicators = list(1, 1)
  ), envir = globalenv())
  on.exit(rm("l_inputs_test", envir = globalenv()), add = TRUE)

  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0
      end   <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  i_block <- input_block(
    base                  = l_inputs_test[["base_value"]],
    psa                   = pick_psa(l_inputs_test[["PSA_dist"]], l_inputs_test[["n"]],
                                     l_inputs_test[["a"]],        l_inputs_test[["b"]]),
    sens                  = l_inputs_test,
    names_out             = l_inputs_test[["parameter_name"]],
    psa_indicators        = l_inputs_test[["psa_indicators"]],
    indicator_sens_binary = TRUE,
    dsa_names             = c("DSA_min", "DSA_max")
  )

  # n_sensitivity NOT passed — engine should auto-set it to 2 from metadata
  results <- run_sim(
    npats             = 5,
    n_sim             = 1,
    psa_bool          = FALSE,
    sensitivity_bool  = TRUE,
    arm_list          = c("arm1"),
    common_all_inputs = i_block,
    init_event_list   = make_evt_list(),
    evt_react_list    = make_react_list(),
    util_ongoing_list = "util.sick",
    sensitivity_names = c("DSA_min", "DSA_max"),
    seed              = 42
  )

  # 2 params x 2 DSA directions = 4 sensitivity iterations
  expect_length(results, 4L)
})

test_that("input_block: base case equivalent to add_item", {
  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0
      end   <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  manual_inputs <- add_item(util.sick = 0.8, util.sicker = 0.5)
  i_block <- input_block(
    base                  = list(0.8, 0.5),
    psa                   = list(0.8, 0.5),
    sens                  = list(),
    names_out             = c("util.sick", "util.sicker"),
    indicator_sens_binary = TRUE
  )

  r_manual <- run_sim(
    npats = 5, n_sim = 1, psa_bool = FALSE, sensitivity_bool = FALSE,
    arm_list = c("arm1"), common_all_inputs = manual_inputs,
    init_event_list = make_evt_list(), evt_react_list = make_react_list(),
    util_ongoing_list = "util.sick", seed = 42
  )
  r_block <- run_sim(
    npats = 5, n_sim = 1, psa_bool = FALSE, sensitivity_bool = FALSE,
    arm_list = c("arm1"), common_all_inputs = i_block,
    init_event_list = make_evt_list(), evt_react_list = make_react_list(),
    util_ongoing_list = "util.sick", seed = 42
  )

  expect_equal(r_manual[[1]][[1]]$total_qalys, r_block[[1]][[1]]$total_qalys)
})

test_that("input_block: DSA each iteration selects correct values", {
  # Verify that input_block correctly routes values for each DSA iteration
  # by comparing against reference runs with the exact expected input values.
  # 2 params x 2 DSA directions = 4 iterations:
  #   iter 1: DSA_min, param 1 -> util.sick=0.6, util.sicker=0.5 (base)
  #   iter 2: DSA_min, param 2 -> util.sick=0.8 (base), util.sicker=0.3
  #   iter 3: DSA_max, param 1 -> util.sick=0.9, util.sicker=0.5 (base)
  #   iter 4: DSA_max, param 2 -> util.sick=0.8 (base), util.sicker=0.7
  assign("l_dsa_test", list(
    DSA_min = list(0.6, 0.3),
    DSA_max = list(0.9, 0.7)
  ), envir = globalenv())
  on.exit(rm("l_dsa_test", envir = globalenv()), add = TRUE)

  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0
      end   <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  i_block <- input_block(
    base                  = list(0.8, 0.5),
    psa                   = list(0.8, 0.5),
    sens                  = l_dsa_test,
    names_out             = c("util.sick", "util.sicker"),
    indicator_sens_binary = TRUE,
    dsa_names             = c("DSA_min", "DSA_max")
  )

  r_block <- run_sim(
    npats = 5, n_sim = 1, psa_bool = FALSE, sensitivity_bool = TRUE,
    arm_list = c("arm1"),
    common_all_inputs = i_block,
    init_event_list = make_evt_list(), evt_react_list = make_react_list(),
    util_ongoing_list = "util.sick",
    sensitivity_names = c("DSA_min", "DSA_max"),
    seed = 42
  )

  expect_length(r_block, 4L)

  # Each iteration's reference: simple add_item with the expected exact values
  expected <- list(
    list(util.sick = 0.6, util.sicker = 0.5),
    list(util.sick = 0.8, util.sicker = 0.3),
    list(util.sick = 0.9, util.sicker = 0.5),
    list(util.sick = 0.8, util.sicker = 0.7)
  )

  for (idx in seq_along(expected)) {
    v <- expected[[idx]]
    r_ref <- run_sim(
      npats = 5, n_sim = 1, psa_bool = FALSE, sensitivity_bool = FALSE,
      arm_list = c("arm1"),
      common_all_inputs = eval(bquote(add_item(
        util.sick = .(v$util.sick), util.sicker = .(v$util.sicker)
      ))),
      init_event_list = make_evt_list(), evt_react_list = make_react_list(),
      util_ongoing_list = "util.sick",
      seed = 42
    )
    expect_equal(r_block[[idx]][[1]]$total_qalys, r_ref[[1]][[1]]$total_qalys,
                 label = paste("iteration", idx))
  }
})

test_that("input_block: equivalent to pick_val_v across psa/sensitivity/n_sim combinations", {
  assign("l_pvv", list(
    parameter_name = list("util.sick", "util.sicker"),
    base_value     = list(0.8, 0.5),
    DSA_min        = list(0.6, 0.3),
    DSA_max        = list(0.9, 0.7),
    PSA_dist       = list("rnorm", "rnorm"),
    psa_a          = list(0.8, 0.5),
    psa_b          = list(0.05, 0.05),
    n              = list(1, 1)
  ), envir = globalenv())
  on.exit(rm("l_pvv", envir = globalenv()), add = TRUE)

  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0; end <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  i_sensitivity_old <- add_item(
    iterator_sensitivity = sens_iterator(sens, n_sensitivity)
  ) |> add_item(
    indicators = if (sensitivity_bool) {
      create_indicators(iterator_sensitivity,
                        n_sensitivity * length(sensitivity_names),
                        rep(1, length(l_pvv[["base_value"]])))
    } else {
      rep(1, length(l_pvv[["base_value"]]))
    }
  )

  i_old <- add_item(
    pick_val_v(
      base      = l_pvv[["base_value"]],
      psa       = pick_psa(l_pvv[["PSA_dist"]], l_pvv[["n"]],
                           l_pvv[["psa_a"]], l_pvv[["psa_b"]]),
      sens      = l_pvv[[sens_name_used]],
      psa_ind   = psa_bool,
      sens_ind  = sensitivity_bool,
      indicator = indicators,
      names_out = l_pvv[["parameter_name"]]
    )
  )

  i_block <- input_block(
    base                  = list(0.8, 0.5),
    psa                   = pick_psa(l_pvv[["PSA_dist"]], l_pvv[["n"]],
                                     l_pvv[["psa_a"]], l_pvv[["psa_b"]]),
    sens                  = l_pvv,
    names_out             = c("util.sick", "util.sicker"),
    indicator_sens_binary = TRUE,
    dsa_names             = c("DSA_min", "DSA_max")
  )

  run_both <- function(psa_b, sens_b, n_s) {
    sn  <- if (sens_b) c("DSA_min", "DSA_max") else NULL
    nse <- if (sens_b) length(l_pvv[["base_value"]]) else 1L

    r_old <- run_sim(
      npats = 5, n_sim = n_s, psa_bool = psa_b, sensitivity_bool = sens_b,
      arm_list           = c("arm1"),
      sensitivity_inputs = i_sensitivity_old,
      common_all_inputs  = i_old,
      init_event_list    = make_evt_list(),
      evt_react_list     = make_react_list(),
      util_ongoing_list  = "util.sick",
      sensitivity_names  = sn,
      n_sensitivity      = nse,
      seed = 42
    )
    r_block <- run_sim(
      npats = 5, n_sim = n_s, psa_bool = psa_b, sensitivity_bool = sens_b,
      arm_list           = c("arm1"),
      common_all_inputs  = i_block,
      init_event_list    = make_evt_list(),
      evt_react_list     = make_react_list(),
      util_ongoing_list  = "util.sick",
      sensitivity_names  = sn,
      seed = 42
    )
    expect_length(r_old, length(r_block))
    for (s in seq_along(r_old)) {
      for (sim in seq_along(r_old[[s]])) {
        expect_equal(
          r_old[[s]][[sim]]$total_qalys,
          r_block[[s]][[sim]]$total_qalys,
          label = sprintf("psa=%s sens=%s n_sim=%d s=%d sim=%d",
                          psa_b, sens_b, n_s, s, sim)
        )
      }
    }
  }

  run_both(psa_b = FALSE, sens_b = FALSE, n_s = 1)  # base case
  run_both(psa_b = TRUE,  sens_b = FALSE, n_s = 2)  # PSA only, 2 sims
  run_both(psa_b = FALSE, sens_b = TRUE,  n_s = 2)  # DSA only, 2 sims
  run_both(psa_b = TRUE,  sens_b = TRUE,  n_s = 3)  # probabilistic DSA, 3 sims
})

test_that("input_block: non-binary - scalar params varied independently", {
  assign("l_nb_indep", list(
    DSA_min = list(0.6, 0.3),
    DSA_max = list(0.9, 0.7)
  ), envir = globalenv())
  on.exit(rm("l_nb_indep", envir = globalenv()), add = TRUE)

  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0; end <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  i_block <- input_block(
    base            = list(0.8, 0.5),
    psa             = list(0.8, 0.5),
    sens            = l_nb_indep,
    names_out       = c("util.sick", "util.sicker"),
    sens_indicators = list(1L, 2L),
    distributions   = list("rnorm", "rnorm"),
    dsa_names       = c("DSA_min", "DSA_max")
  )

  i_old <- add_item(
    pick_val_v(
      base      = list(0.8, 0.5),
      psa       = list(0.8, 0.5),
      sens      = l_nb_indep[[sens_name_used]],
      psa_ind   = psa_bool,
      sens_ind  = sensitivity_bool,
      indicator = list(1L, 2L),
      indicator_sens_binary = FALSE,
      sens_iterator = sens_iterator(sens, n_sensitivity),
      distributions = list("rnorm", "rnorm"),
      names_out = c("util.sick", "util.sicker")
    )
  )

  run_and_compare <- function(psa_b, n_s) {
    r_block <- run_sim(
      npats = 5, n_sim = n_s, psa_bool = psa_b, sensitivity_bool = TRUE,
      arm_list = c("arm1"), common_all_inputs = i_block,
      init_event_list = make_evt_list(), evt_react_list = make_react_list(),
      util_ongoing_list = "util.sick",
      sensitivity_names = c("DSA_min", "DSA_max"), seed = 42
    )
    r_old <- run_sim(
      npats = 5, n_sim = n_s, psa_bool = psa_b, sensitivity_bool = TRUE,
      arm_list = c("arm1"), common_all_inputs = i_old,
      init_event_list = make_evt_list(), evt_react_list = make_react_list(),
      util_ongoing_list = "util.sick",
      sensitivity_names = c("DSA_min", "DSA_max"), n_sensitivity = 2L, seed = 42
    )
    expect_length(r_block, length(r_old))
    for (s in seq_along(r_old)) {
      for (sim in seq_along(r_old[[s]])) {
        expect_equal(
          r_block[[s]][[sim]]$total_qalys,
          r_old[[s]][[sim]]$total_qalys,
          label = sprintf("psa=%s n_sim=%d s=%d sim=%d", psa_b, n_s, s, sim)
        )
      }
    }
  }

  run_and_compare(FALSE, 1)
  run_and_compare(TRUE, 2)
})

test_that("input_block: non-binary - scalar params varied as a group", {
  assign("l_nb_grp", list(
    DSA_min = list(0.6, 0.3),
    DSA_max = list(0.9, 0.7)
  ), envir = globalenv())
  on.exit(rm("l_nb_grp", envir = globalenv()), add = TRUE)

  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0; end <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  i_block <- input_block(
    base            = list(0.8, 0.5),
    psa             = list(0.8, 0.5),
    sens            = l_nb_grp,
    names_out       = c("util.sick", "util.sicker"),
    sens_indicators = list(1L, 1L),
    distributions   = list("rnorm", "rnorm"),
    dsa_names       = c("DSA_min", "DSA_max")
  )

  i_old <- add_item(
    pick_val_v(
      base      = list(0.8, 0.5),
      psa       = list(0.8, 0.5),
      sens      = l_nb_grp[[sens_name_used]],
      psa_ind   = psa_bool,
      sens_ind  = sensitivity_bool,
      indicator = list(1L, 1L),
      indicator_sens_binary = FALSE,
      sens_iterator = sens_iterator(sens, n_sensitivity),
      distributions = list("rnorm", "rnorm"),
      names_out = c("util.sick", "util.sicker")
    )
  )

  run_and_compare <- function(psa_b, n_s) {
    r_block <- run_sim(
      npats = 5, n_sim = n_s, psa_bool = psa_b, sensitivity_bool = TRUE,
      arm_list = c("arm1"), common_all_inputs = i_block,
      init_event_list = make_evt_list(), evt_react_list = make_react_list(),
      util_ongoing_list = "util.sick",
      sensitivity_names = c("DSA_min", "DSA_max"), n_sensitivity = 1L, seed = 42
    )
    r_old <- run_sim(
      npats = 5, n_sim = n_s, psa_bool = psa_b, sensitivity_bool = TRUE,
      arm_list = c("arm1"), common_all_inputs = i_old,
      init_event_list = make_evt_list(), evt_react_list = make_react_list(),
      util_ongoing_list = "util.sick",
      sensitivity_names = c("DSA_min", "DSA_max"), n_sensitivity = 1L, seed = 42
    )
    expect_length(r_block, length(r_old))
    for (s in seq_along(r_old)) {
      for (sim in seq_along(r_old[[s]])) {
        expect_equal(
          r_block[[s]][[sim]]$total_qalys,
          r_old[[s]][[sim]]$total_qalys,
          label = sprintf("psa=%s n_sim=%d s=%d sim=%d", psa_b, n_s, s, sim)
        )
      }
    }
  }

  run_and_compare(FALSE, 1)
  run_and_compare(TRUE, 2)
})

test_that("input_block: non-binary - vector param (mvrnorm) grouped c(3,3) vs independent c(3,4)", {
  # 3-parameter model: 2 scalars (groups 1, 2) + 1 two-element vector (groups
  # c(3,3) or c(3,4)). Identity covariance => partial DSA uses exact replacement.
  assign("l_mvn_test", list(
    DSA_min = list(0.6, 0.3, c(6.0, 3.0)),
    DSA_max = list(0.9, 0.7, c(9.0, 7.0))
  ), envir = globalenv())
  on.exit(rm("l_mvn_test", envir = globalenv()), add = TRUE)

  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0; end <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  run_and_compare <- function(i_block, i_old, n_sens, label) {
    r_block <- run_sim(
      npats = 5, n_sim = 1, psa_bool = FALSE, sensitivity_bool = TRUE,
      arm_list = c("arm1"), common_all_inputs = i_block,
      init_event_list = make_evt_list(), evt_react_list = make_react_list(),
      util_ongoing_list = "util.sick",
      sensitivity_names = c("DSA_min", "DSA_max"), n_sensitivity = n_sens, seed = 42
    )
    r_old <- run_sim(
      npats = 5, n_sim = 1, psa_bool = FALSE, sensitivity_bool = TRUE,
      arm_list = c("arm1"), common_all_inputs = i_old,
      init_event_list = make_evt_list(), evt_react_list = make_react_list(),
      util_ongoing_list = "util.sick",
      sensitivity_names = c("DSA_min", "DSA_max"), n_sensitivity = n_sens, seed = 42
    )
    expect_length(r_block, length(r_old))
    for (s in seq_along(r_old)) {
      expect_equal(r_block[[s]][[1]]$total_qalys, r_old[[s]][[1]]$total_qalys,
                   label = paste(label, "s", s))
    }
  }

  # Grouped: both elements of v_state change together (c(3,3)), n_sensitivity=3
  i_block_grp <- input_block(
    base            = list(0.8, 0.5, c(8.0, 5.0)),
    psa             = list(0.8, 0.5, c(8.0, 5.0)),
    sens            = l_mvn_test,
    names_out       = c("util.sick", "util.sicker", "v_state"),
    sens_indicators = list(1L, 2L, c(3L, 3L)),
    distributions   = list("rnorm", "rnorm", "mvrnorm"),
    covariances     = list(NULL, NULL, matrix(c(1, 0, 0, 1), 2, 2)),
    dsa_names       = c("DSA_min", "DSA_max")
  )
  i_old_grp <- add_item(
    pick_val_v(
      base      = list(0.8, 0.5, c(8.0, 5.0)),
      psa       = list(0.8, 0.5, c(8.0, 5.0)),
      sens      = l_mvn_test[[sens_name_used]],
      psa_ind   = psa_bool,
      sens_ind  = sensitivity_bool,
      indicator = list(1L, 2L, c(3L, 3L)),
      indicator_sens_binary = FALSE,
      sens_iterator = sens_iterator(sens, n_sensitivity),
      distributions = list("rnorm", "rnorm", "mvrnorm"),
      covariances   = list(NULL, NULL, matrix(c(1, 0, 0, 1), 2, 2)),
      names_out = c("util.sick", "util.sicker", "v_state")
    )
  )

  # Independent: each element of v_state has its own DSA step (c(3,4)), n_sensitivity=4
  i_block_ind <- input_block(
    base            = list(0.8, 0.5, c(8.0, 5.0)),
    psa             = list(0.8, 0.5, c(8.0, 5.0)),
    sens            = l_mvn_test,
    names_out       = c("util.sick", "util.sicker", "v_state"),
    sens_indicators = list(1L, 2L, c(3L, 4L)),
    distributions   = list("rnorm", "rnorm", "mvrnorm"),
    covariances     = list(NULL, NULL, matrix(c(1, 0, 0, 1), 2, 2)),
    dsa_names       = c("DSA_min", "DSA_max")
  )
  i_old_ind <- add_item(
    pick_val_v(
      base      = list(0.8, 0.5, c(8.0, 5.0)),
      psa       = list(0.8, 0.5, c(8.0, 5.0)),
      sens      = l_mvn_test[[sens_name_used]],
      psa_ind   = psa_bool,
      sens_ind  = sensitivity_bool,
      indicator = list(1L, 2L, c(3L, 4L)),
      indicator_sens_binary = FALSE,
      sens_iterator = sens_iterator(sens, n_sensitivity),
      distributions = list("rnorm", "rnorm", "mvrnorm"),
      covariances   = list(NULL, NULL, matrix(c(1, 0, 0, 1), 2, 2)),
      names_out = c("util.sick", "util.sicker", "v_state")
    )
  )

  run_and_compare(i_block_grp, i_old_grp, n_sens = 3L, label = "grouped")
  run_and_compare(i_block_ind, i_old_ind, n_sens = 4L, label = "independent")
})

test_that("input_block: dsa_names - binary mode, DSA iterations cycle params, scenario uses all", {
  # sensitivity_names = c("DSA_min","DSA_max","scenario_1")
  # dsa_names = c("DSA_min","DSA_max"): only these cycle through params 1-at-a-time.
  # scenario_1: all params should take the scenario values simultaneously.
  # Expected: 2 params x 2 DSA names + 1 scenario = 5 total iterations.
  assign("l_mixed", list(
    DSA_min    = list(0.6, 0.3),
    DSA_max    = list(0.9, 0.7),
    scenario_1 = list(0.55, 0.45)
  ), envir = globalenv())
  on.exit(rm("l_mixed", envir = globalenv()), add = TRUE)

  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0; end <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  i_block <- input_block(
    base                  = list(0.8, 0.5),
    psa                   = list(0.8, 0.5),
    sens                  = l_mixed,
    names_out             = c("util.sick", "util.sicker"),
    dsa_names             = c("DSA_min", "DSA_max"),
    indicator_sens_binary = TRUE
  )

  sn <- c("DSA_min", "DSA_max", "scenario_1")
  r_block <- run_sim(
    npats = 5, n_sim = 1, psa_bool = FALSE, sensitivity_bool = TRUE,
    arm_list = c("arm1"), common_all_inputs = i_block,
    init_event_list = make_evt_list(), evt_react_list = make_react_list(),
    util_ongoing_list = "util.sick",
    sensitivity_names = sn, seed = 42
  )

  # 2 params x 2 DSA names + 1 scenario = 5 total iterations
  expect_length(r_block, 5L)
  qalys <- vapply(r_block, \(s) s[[1]]$total_qalys, numeric(1))
  # DSA_min iter1 (util.sick=0.6) < DSA_min iter2 (util.sick=base=0.8)
  expect_lt(qalys[1], qalys[2])
  # DSA_max iter1 (util.sick=0.9) > DSA_max iter2 (util.sick=base=0.8)
  expect_gt(qalys[3], qalys[4])
  # Both DSA iter2s use base util.sick → equal QALYs
  expect_equal(qalys[2], qalys[4])
  # scenario_1 uses util.sick=0.55 < base=0.8
  expect_lt(qalys[5], qalys[2])
})

test_that("input_block: dsa_names - grouped mode, DSA iterates groups, scenario uses all", {

  # sens_indicators=list(1L,2L): each param is its own group (independent).
  # Expected: 2 groups x 2 DSA names + 1 scenario = 5 total iterations.
  assign("l_mixed_grp", list(
    DSA_min    = list(0.6, 0.3),
    DSA_max    = list(0.9, 0.7),
    scenario_1 = list(0.55, 0.45)
  ), envir = globalenv())
  on.exit(rm("l_mixed_grp", envir = globalenv()), add = TRUE)

  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0; end <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  i_block <- input_block(
    base            = list(0.8, 0.5),
    psa             = list(0.8, 0.5),
    sens            = l_mixed_grp,
    names_out       = c("util.sick", "util.sicker"),
    sens_indicators = list(1L, 2L),
    distributions   = list("rnorm", "rnorm"),
    dsa_names       = c("DSA_min", "DSA_max")
  )

  sn <- c("DSA_min", "DSA_max", "scenario_1")
  r_block <- run_sim(
    npats = 5, n_sim = 1, psa_bool = FALSE, sensitivity_bool = TRUE,
    arm_list = c("arm1"), common_all_inputs = i_block,
    init_event_list = make_evt_list(), evt_react_list = make_react_list(),
    util_ongoing_list = "util.sick",
    sensitivity_names = sn, seed = 42
  )

  # 2 groups x 2 DSA names + 1 scenario = 5 total iterations
  expect_length(r_block, 5L)
  qalys <- vapply(r_block, \(s) s[[1]]$total_qalys, numeric(1))
  # DSA_min iter1 (util.sick=0.6) < DSA_min iter2 (util.sick=base=0.8)
  expect_lt(qalys[1], qalys[2])
  # DSA_max iter1 (util.sick=0.9) > DSA_max iter2 (util.sick=base=0.8)
  expect_gt(qalys[3], qalys[4])
  # Both DSA iter2s use base util.sick → equal QALYs
  expect_equal(qalys[2], qalys[4])
  # scenario_1 uses util.sick=0.55 < base=0.8
  expect_lt(qalys[5], qalys[2])
})

test_that("input_block: dsa_names=NULL treats all sensitivity_names as scenarios (1 iter each)", {
  assign("l_scen", list(
    scen_a = list(0.6, 0.3),
    scen_b = list(0.9, 0.7)
  ), envir = globalenv())
  on.exit(rm("l_scen", envir = globalenv()), add = TRUE)

  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0; end <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  # dsa_names = NULL (default) → all names are scenarios
  i_block <- input_block(
    base      = list(0.8, 0.5),
    psa       = list(0.8, 0.5),
    sens      = l_scen,
    names_out = c("util.sick", "util.sicker")
  )

  r_block <- run_sim(
    npats = 5, n_sim = 1, psa_bool = FALSE, sensitivity_bool = TRUE,
    arm_list = c("arm1"), common_all_inputs = i_block,
    init_event_list = make_evt_list(), evt_react_list = make_react_list(),
    util_ongoing_list = "util.sick",
    sensitivity_names = c("scen_a", "scen_b"), seed = 42
  )

  # 2 scenario names → 2 iterations (not 2 params x 2 = 4)
  expect_length(r_block, 2L)
  qalys <- vapply(r_block, \(s) s[[1]]$total_qalys, numeric(1))
  # scen_a: util.sick=0.6 < scen_b: util.sick=0.9
  expect_lt(qalys[1], qalys[2])
})

test_that("input_block: sens_indicators with 0 skips that parameter in DSA", {
  assign("l_skip", list(
    DSA_min = list(0.6, 0.3),
    DSA_max = list(0.9, 0.7)
  ), envir = globalenv())
  on.exit(rm("l_skip", envir = globalenv()), add = TRUE)

  make_evt_list <- function() {
    add_tte(arm = c("arm1"), evts = c("start", "end"), input = {
      start <- 0; end <- 5
    })
  }
  make_react_list <- function() {
    add_reactevt(name_evt = "start", input = {}) |>
      add_reactevt(name_evt = "end", input = { curtime <- Inf })
  }

  # sens_indicators = list(1L, 0L): only param 1 varies, param 2 is fixed
  i_block <- input_block(
    base                  = list(0.8, 0.5),
    psa                   = list(0.8, 0.5),
    sens                  = l_skip,
    names_out             = c("util.sick", "util.sicker"),
    sens_indicators       = list(1L, 0L),
    indicator_sens_binary = TRUE,
    dsa_names             = c("DSA_min", "DSA_max")
  )

  # n_groups=1 (only param 1 active) → 1 iter per DSA name = 2 total
  r_block <- run_sim(
    npats = 5, n_sim = 1, psa_bool = FALSE, sensitivity_bool = TRUE,
    arm_list = c("arm1"), common_all_inputs = i_block,
    init_event_list = make_evt_list(), evt_react_list = make_react_list(),
    util_ongoing_list = "util.sick",
    sensitivity_names = c("DSA_min", "DSA_max"), seed = 42
  )

  expect_length(r_block, 2L)
  qalys <- vapply(r_block, \(s) s[[1]]$total_qalys, numeric(1))
  # DSA_min: util.sick=0.6 (lower), DSA_max: util.sick=0.9 (higher)
  expect_lt(qalys[1], qalys[2])
})
