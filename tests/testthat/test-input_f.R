test_that("Vector discounting equal to non-vector for single-length elements", {
  #Ongoing
  expect_equal(disc_ongoing_v(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500),
               disc_ongoing(lcldr=0.035,lclprvtime=0.5, lclcurtime=3, lclval=2500))
  
  #Instant
  expect_equal(disc_instant_v(lcldr=0.035, lclcurtime=3, lclval=2500),
               disc_instant(lcldr=0.035, lclcurtime=3, lclval=2500))

  #Cycle
  expect_equal(disc_cycle_v(lcldr=0.035, lclcurtime=3, lclval=2500,lclprvtime=0, cyclelength=1/12,starttime=0),
               disc_cycle(lcldr=0.035, lclcurtime=3, lclval=2500,lclprvtime=0, cyclelength=1/12,starttime=0))  
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
  expect_equal(disc_cycle(lcldr=0,lclprvtime=1, lclcurtime=2, lclval=2500, cyclelength=1/12, starttime=0),
               12*2500)
  expect_equal(disc_cycle_v(lcldr=0,lclprvtime=0, lclcurtime=0, lclval=2500, cyclelength=1/12, starttime=0),
               2500)
  expect_equal(disc_cycle(lcldr=0,lclprvtime=0, lclcurtime=0, lclval=2500, cyclelength=1/12, starttime=0),
               2500)
  expect_equal(disc_cycle_v(lcldr=0.035,lclprvtime=0, lclcurtime=0, lclval=2500, cyclelength=2, starttime=0),
               2500)
  expect_equal(disc_cycle(lcldr=0.035,lclprvtime=0, lclcurtime=0, lclval=2500, cyclelength=2, starttime=0),
               2500)
  expect_equal(disc_cycle_v(lcldr=0.035,lclprvtime=4, lclcurtime=5, lclval=2500, cyclelength=1/12, starttime=4.5),
               disc_cycle(lcldr=0.035,lclprvtime=4, lclcurtime=5, lclval=2500, cyclelength=1/12, starttime=4.5))
  #Inf*0 gives NaN
  expect_equal(disc_cycle_v(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=Inf, cyclelength=1/12, starttime=0),
               NaN)
  expect_equal(disc_cycle(lcldr=0.035,lclprvtime=5, lclcurtime=5, lclval=Inf, cyclelength=1/12, starttime=0),
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
      
      rnd_prob_exn_sev = runif_stream(1, substream_prob_exn_sev),
      
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
      
      rnd_exn = runif_stream(1, substream_exn)
      
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
  
  expect_equal(nrow(extract_elements_from_list(ast_as_list(expr))),39) #39 items/events changed
  
  a <- add_reactevt(name_evt="example",
                    input={
                      a <- 5
                      modify_item(list(w=5))
                    })
  
  
  expect_equal(nrow(extract_from_reactions(a)),1) #1 items/events changed
  
  expect_equal(extract_from_reactions(a),
               data.table(event = "example",
                          name = "w",
                          type = "item",
                          conditional_flag = FALSE,
                          definition = "5")
               ) 
  
  
  
})
