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

