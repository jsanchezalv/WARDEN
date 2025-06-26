test_that("Vector discounting equal to non-vector for single-length elements", {
  #Ongoing
  expect_equal(disc_ongoing_v(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500),
               disc_ongoing(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500))
  
  #Instant
  expect_equal(disc_instant_v(lcldr=0.035, lclcurtime=3, lclval=2500),
               disc_instant(lcldr=0.035, lclcurtime=3, lclval=2500))

 })


test_that("Discounting works with no discounting", {
  #Ongoing
  expect_equal(disc_ongoing_v(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500),
               2500)
  expect_equal(disc_ongoing(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500),
               2500)
  
  #Instant
  expect_equal(disc_instant_v(lcldr=0, lclcurtime=2, lclval=2500),
               2500)
  expect_equal(disc_instant(lcldr=0, lclcurtime=2, lclval=2500),
               2500)
  
  #Cycle
  expect_equal(disc_cycle_v(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500, cyclelength=1/12, starttime=0),
               12*2500)
  expect_equal(disc_cycle(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500, cyclelength=1/12, starttime=0),
               12*2500)
})



test_that("Discounting works with odd numbers", {
  #Ongoing
  expect_equal(disc_ongoing_v(lcldr=0,lclprvtime=0, lclcurtime=Inf, lclval=2500),
               Inf)
  expect_equal(disc_ongoing(lcldr=0,lclprvtime=0, lclcurtime=Inf, lclval=2500),
               Inf)
  expect_equal(disc_ongoing_v(lcldr=0.035,lclprvtime=0, lclcurtime=0, lclval=2500),
               0)
  expect_equal(disc_ongoing(lcldr=0.035,lclprvtime=0, lclcurtime=0, lclval=2500),
               0)
  expect_equal(disc_ongoing_v(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=2500),
               0)
  expect_equal(disc_ongoing(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=2500),
               0)
  #Inf*0 gives NaN, while the element-wise function just check wehterh prevtime and curtime are equal
  expect_equal(disc_ongoing_v(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=Inf),
               NaN)
  expect_equal(disc_ongoing(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=Inf),
               NaN) 
  
  #Instant
  expect_equal(disc_instant_v(lcldr=0, lclcurtime=Inf, lclval=2500),
               2500)
  expect_equal(disc_instant(lcldr=0, lclcurtime=Inf, lclval=2500),
               2500)
  expect_equal(disc_instant_v(lcldr=0.035, lclcurtime=0, lclval=2500),
               2500)
  expect_equal(disc_instant(lcldr=0.035, lclcurtime=0, lclval=2500),
               2500)
  expect_equal(disc_instant_v(lcldr=0.035, lclcurtime=5, lclval=2500),
               2104.93292)
  expect_equal(disc_instant(lcldr=0.035, lclcurtime=5, lclval=2500),
               2104.93292)
  #Inf*0 gives NaN
  expect_equal(disc_instant_v(lcldr=0.035, lclcurtime=5, lclval=Inf),
               Inf)
  expect_equal(disc_instant(lcldr=0.035, lclcurtime=5, lclval=Inf),
               Inf) 
  expect_equal(disc_instant_v(lcldr=5, lclcurtime=5, lclval=2500),
               0.32150206)
  expect_equal(disc_instant(lcldr=5, lclcurtime=5, lclval=2500),
               0.32150206) 
  
  #Cycle
  expect_equal(disc_cycle_v(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500, cyclelength=1/12, starttime=1),
               12*2500)
  expect_equal(disc_cycle_v(lcldr=0,lclprvtime=0, lclcurtime=0, lclval=2500, cyclelength=1/12, starttime=0),
               2500)
  expect_equal(disc_cycle_v(lcldr=0.035,lclprvtime=0, lclcurtime=0, lclval=2500, cyclelength=2, starttime=0),
               2500)
  #Inf*0 gives NaN
  expect_equal(disc_cycle_v(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=Inf, cyclelength=1/12, starttime=0),
               NaN)
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
               c(0,28215.4394,Inf))  
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
               covariances = list(0.1,0.1,NULL) ),
    list(util = 2, util2 = 3, dirichlet_vector = c(0.36, 0.54, 0.1))
    
  )
  
  expect_equal(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = FALSE,
      sens_ind = FALSE,
      indicator=c(1,0)
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
      indicator=c(1,0)
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
      indicator=c(0,1)
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
      indicator=c(1,1)
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
      indicator=c(1,1)
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
      indicator=c(1,5)
    )
  )
  
  expect_error(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = 5,
      sens_ind = TRUE,
      indicator=c(1,1)
    )
  )
  
  expect_error(
    pick_val_v(
      base = c(0,0),
      psa =c(rnorm(1,0,0.1),rnorm(1,0,0.1)),
      sens = c(2,3),
      psa_ind = TRUE,
      sens_ind = 3,
      indicator=c(1,1)
    )
  )
  
  expect_equal(
    pick_val_v(
      base = c(0,0),
      psa ={c(draw_tte(1,'norm',0,0.1,seed=1),draw_tte(1,'norm',0,0.1,seed=2))},
      sens = c(2,3),
      psa_ind = TRUE,
      sens_ind = TRUE,
      indicator=c(0,0)
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
      indicator=c(1,0)
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
    
    modify_item(list(afsa=ifelse(TRUE,"asda",NULL)))
    
    modify_item_seq(list(
      
      o_exn = o_exn + 1,
      
      a = NULL,
      
      b = if(a){"CZ"}else{"AW"},
      
      rnd_prob_exn_sev = runif(1),
      
      exn_sev = rnd_prob_exn_sev <= p_sev,
      
      o_exn_mod = o_exn_mod + if(exn_sev) { 0 } else { 1 },
      
      o_exn_sev = o_exn_sev + if(exn_sev) { 1 } else { 0 },
      
      o_rec_time_without_exn = (o_exn == 0) * 1,
      
      o_rec_time_without_exn_sev = (o_exn_sev == 0) * 1,
      
      o_c_exn = if(exn_sev) { c_sev } else { c_mod },
      
      o_other_c_exn_mod = if(exn_sev) { 0 } else { c_mod },
      
      o_other_c_exn_sev = if(exn_sev) { c_sev } else { 0 },
      
      o_qloss_exn = -if(exn_sev) { q_sev } else { q_mod },
      
      o_other_qloss_exn_mod = -if(exn_sev) { 0 } else { q_mod },
      
      o_other_qloss_exn_sev = -if(exn_sev) { q_sev } else { 0 },
      
      o_qloss_cg_exn = -if(exn_sev) { q_cg_sev } else { q_cg_mod },
      
      o_other_qloss_cg_exn_mod = -if(exn_sev) { 0 } else { q_cg_mod },
      
      o_other_qloss_cg_exn_sev = -if(exn_sev) { q_cg_sev } else { 0 },
      
      o_q = utility,
      
      o_other_q_gold1 = if(gold == 1) { utility } else { 0 },
      
      o_other_q_gold2 = if(gold == 2) { utility } else { 0 },
      
      o_other_q_gold3 = if(gold == 3) { utility } else { 0 },
      
      o_other_q_gold4 = if(gold == 4) { utility } else { 0 },
      
      o_other_q_on_dup = if(on_dup) { utility } else { 0 },
      
      n_exn = n_exn + 1,
      
      n_exn_mod = n_exn_mod + (1 - exn_sev),
      
      n_exn_sev = n_exn_sev + exn_sev,
      
      u_adj_exn_lt = u_adj_exn_lt + if(exn_sev) { u_adj_sev_lt } else { u_adj_mod_lt },
      
      utility = u_gold - u_adj_exn_lt - u_mace_lt,
      
      o_rec_utility = utility,
      
      rnd_exn = runif(1)
      
    ))
    
    if(a==1){
      modify_item(list(a=list(6+b)))
      
      modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
    } else{
      modify_event(list(e_exn = curtime + 14 / days_in_year + qexp(rnd_exn, r_exn)))
      if(a>6){
        modify_item(list(a=8))
      }
      
    }
    
    
    if (sel_resp_incl == 1 & on_dup == 1) {
      
      modify_event(list(e_response = curtime, z = 6))
      
    }
    
  })
  
  expect_length(ast_as_list(expr),13) 
  
  expect_type(ast_as_list(expr),"list")
  
  expect_equal(class(extract_elements_from_list(ast_as_list(expr))),"data.frame")
  
  expect_length(extract_elements_from_list(ast_as_list(expr)),4) #4 columns
  
  expect_equal(nrow(extract_elements_from_list(ast_as_list(expr))),43) #43 items/events changed, including assignments
  
  a <- add_reactevt(name_evt="example",
                    input={
                      a <- 5
                      modify_item(list(w=5))
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

test_that("modify_item modifies input items correctly", {
  input_list_arm <- environment()
  input_list_arm$qaly_default_instant = 100
  input_list_arm$debug = FALSE
  input_list_arm$accum_backwards = FALSE

  modify_item(list("qaly_default_instant" = 200))
  expect_equal(input_list_arm$qaly_default_instant, 200)
  
  modify_item(list(new_cost = 300))
  expect_equal(input_list_arm$new_cost, 300)
})

test_that("modify_event modifies events correctly", {
  input_list_arm <- environment()
  input_list_arm$qaly_default_instant = 100
  input_list_arm$debug = FALSE
  input_list_arm$cur_evtlist = FALSE
  input_list_arm$accum_backwards = FALSE
  input_list_arm$cur_evtlist = c(ae = 5, nat.death = 100)


  # Modify an existing event
  modify_event(list(ae = 10))
  expect_equal(input_list_arm$cur_evtlist[["ae"]], 10)
  
  # Create new event if not exists
  modify_event(list(new_event = 50), create_if_null = TRUE)
  expect_equal(input_list_arm$cur_evtlist[["new_event"]], 50)
  
  # Ignore non-existent event
  expect_warning(modify_event(list(nonexistent = 20), create_if_null = FALSE))
  expect_error(input_list_arm$cur_evtlist[["nonexistent"]])
})

test_that("new_event adds new events correctly", {
  input_list_arm <- environment()
  input_list_arm$qaly_default_instant = 100
  input_list_arm$debug = FALSE
  input_list_arm$cur_evtlist = FALSE
  input_list_arm$accum_backwards = FALSE
  input_list_arm$cur_evtlist = c()
  
  new_event(list("ae" = 5))
  expect_equal(input_list_arm$cur_evtlist[["ae"]], 5)
  
  expect_error(new_event(list("not_numeric" = "five")), 
               "New event times are not all numeric, please review")
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


test_that("modify_item_seq works sequentially", {
  
  input_list_arm <- environment()
  input_list_arm$a <- 1
  input_list_arm$b <- 2
  input_list_arm$curtime <- 1
  input_list_arm$accum_backwards <- FALSE
  input_list_arm$debug <- FALSE
  
  # Test sequential modification
  modify_item_seq(list(a = 3, b = a + 2))
  expect_equal(input_list_arm$a, 3)
  expect_equal(input_list_arm$b, 5)
  
  # Test debug mode
  input_list_arm$debug <- TRUE
  input_list_arm$log_list <- list()
  modify_item_seq(list(a = 4, c = b * 2))
  expect_equal(input_list_arm$a, 4)
  expect_equal(input_list_arm$c, 10)
  expect_true(length(input_list_arm$log_list) > 0)
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
  expect_type(tte, "double")
  expect_length(tte, 1)
  expect_gt(tte, 0)
  expect_lt(tte, 100)
})

test_that("qtimecov works for all supported distributions", {
  set.seed(1)
  param_fun_factory <- function(p0, p1, p2, p3) {
    function(t) p0 + p1 * t + p2 * t^2 + p3 * (floor(t) + 1)
  }
  
  # 1. Exponential
  rate_exp <- param_fun_factory(0.1, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = rate_exp, dist = "exp"))
  expect_equal(qtimecov(luck = 0.5,a_fun = rate_exp,dist = "exp", dt = 0.001
  ),qexp(0.5,0.1), tolerance = 0.01)
  
  # 2. Gamma
  shape <- param_fun_factory(2, 0, 0, 0)
  rate <- param_fun_factory(0.2, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = shape, b_fun = rate, dist = "gamma"))
  expect_equal(qtimecov(luck = 0.5,a_fun = shape, b_fun = rate, dist = "gamma", dt = 0.001
  ),qgamma(0.5,2,0.2), tolerance = 0.01)
  
  # 3. Lognormal
  meanlog <- param_fun_factory(log(10) - 0.5^2 / 2, 0, 0, 0)
  sdlog <- param_fun_factory(0.5, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = meanlog, b_fun = sdlog, dist = "lnorm"))
  expect_equal(qtimecov(0.5, a_fun = meanlog, b_fun = sdlog, dist = "lnorm",dt=0.01),
               qlnorm(0.5,log(10) - 0.5^2 / 2,0.5), tolerance = 0.01)
  
  # 4. Normal
  mean <- param_fun_factory(10, 0, 0, 0)
  sd <- param_fun_factory(2, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = mean, b_fun = sd, dist = "norm"))
  expect_equal(qtimecov(0.5, a_fun = mean, b_fun = sd, dist = "norm",dt=0.01),
               qnorm(0.5,10,2), tolerance = 0.01)
  
  # 5. Weibull
  shape <- param_fun_factory(2, 0, 0, 0)
  scale <- param_fun_factory(10, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = shape, b_fun = scale, dist = "weibull"))
  expect_equal(qtimecov(0.5, a_fun = shape, b_fun = scale, dist = "weibull",dt=0.01),
               qweibull(0.5,2,10), tolerance = 0.01)
  
  # 6. Loglogistic
  shape <- param_fun_factory(2.5, 0, 0, 0)
  scale <- param_fun_factory(7.6, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = shape, b_fun = scale, dist = "llogis"))
  expect_equal(qtimecov(0.5, a_fun = shape, b_fun = scale, dist = "llogis",dt=0.01),
               flexsurv::qllogis(0.5,2.5,7.6), tolerance = 0.01)
  
  # 7. Gompertz
  shape <- param_fun_factory(0.01, 0, 0, 0)
  rate <- param_fun_factory(0.091, 0, 0, 0)
  expect_silent(qtimecov(runif(1), a_fun = shape, b_fun = rate, dist = "gompertz"))
  expect_equal(qtimecov(0.5, a_fun = shape, b_fun = rate, dist = "gompertz",dt=0.01),
               qgompertz(0.5,0.01,0.091), tolerance = 0.01)
  
  rate_exp <- function(t) 0.1
  init_luck <- 0.95
  expect_equal(qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001),{
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001, return_luck = TRUE,max_time = 10)
    
    qtimecov(luck = a$luck,a_fun = rate_exp,dist = "exp", dt = 0.001, start_time=a$tte)},
    tolerance = 0.01)
  
  init_luck <- 0.99
  expect_equal(qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001),{
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001, return_luck = TRUE,max_time = 5)
    
    qtimecov(luck = a$luck,a_fun = rate_exp,dist = "exp", dt = 0.001, start_time=a$tte)},
    tolerance = 0.01)
  
  
  init_luck <- 0.3
  expect_equal(qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001),{
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.001, return_luck = TRUE,max_time = 1)
    
    qtimecov(luck = a$luck,a_fun = rate_exp,dist = "exp", dt = 0.001, start_time=a$tte)},
    tolerance = 0.01)
  
  
  rate_exp <- function(t) 0.1
  rate_exp2 <- function(t) 0.2
  time_change <- 10
  init_luck <- 0.95
  
  expect_equal({
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005, max_time = time_change, return_luck = TRUE)
    qtimecov(luck = a$luck,a_fun = rate_exp2,dist = "exp", dt = 0.005, start_time=a$tte)
    
  },{
    new_luck <- luck_adj(prevsurv = 1 - pexp(q=time_change,rate_exp(1)),
                         cursurv = 1 - pexp(q=time_change,rate_exp2(1)),
                         luck = init_luck,
                         condq = FALSE) #time 10 change
    qexp(new_luck,rate_exp2(1))
  }, tolerance = 0.01)
  
  

# time varying and an event -----------------------------------------------

  rate_exp <- function(t) 0.1 + 0.01*t * 0.00001*t^2
  rate_exp2 <- function(t) 0.2 + 0.02*t
  time_change <- 8
  init_luck <- 0.95
  
  expect_equal({
    a <- qtimecov(luck = init_luck,a_fun = rate_exp,dist = "exp", dt = 0.005, max_time = time_change, return_luck = TRUE)
    qtimecov(luck = a$luck,a_fun = rate_exp2,dist = "exp", dt = 0.005, start_time=a$tte)
    
  },{# Manually reproduce what qtimecov does from a$tte
    t <- 0
    luck <- init_luck
    dt <- 0.005
    repeat {
      t <- t + dt
      s_prev <- 1 - pexp(t - dt, rate = rate_exp(t - dt))
      s_curr <- 1 - pexp(t,       rate = rate_exp(t))
      
      luck <- luck_adj(prevsurv = s_prev, cursurv = s_curr, luck = luck, condq = TRUE)
      
      res_tte <- qcond_exp(luck, rate = rate_exp(t))
      total_tte <- t - dt + res_tte
      
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
      s_curr <- 1 - pexp(t,       rate = rate_exp2(t))
      
      luck <- luck_adj(prevsurv = s_prev, cursurv = s_curr, luck = luck, condq = TRUE)
      
      res_tte <- qcond_exp(luck, rate = rate_exp2(t))
      total_tte <- t - dt + res_tte
      
      if (res_tte <= dt || total_tte <= t) {
        break
      }
    }
    
    total_tte
  }, tolerance = 0.01)  
  
})



test_that("qtimecov throws error for unsupported distribution", {
  dummy <- function(t) 1
  expect_error(qtimecov(runif(1), a_fun = dummy, dist = "beta"), "Unsupported distribution")
})

test_that("qtimecov respects max_time bound", {
  slow_fun <- function(t) 0.00001
  tte <- qtimecov(
    luck = 0.999,
    a_fun = slow_fun,
    dist = "exp",
    max_time = 10
  )
  expect_lte(tte, 10)
})
