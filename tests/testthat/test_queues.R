# Run all unit tests for EventQueue
  cat("Running EventQueue Unit Tests...\n\n")
  
  test_that("Queue creation works", {
    event_priority <- c("death", "dropout", "visit")
    q <- queue_create(event_priority)
    
    expect_true(!is.null(q))
    expect_true(queue_empty(q))
    expect_equal(queue_size(q), 0L)
  })
  
  test_that("Single event operations work", {
    q <- queue_create(c("death", "dropout", "visit"))
    new_event2(c(death = 11), ptr = q, patient_id = 10)
    
    expect_false(queue_empty(q))
    expect_equal(queue_size(q), 1L)
    expect_true(has_event("death", ptr = q, patient_id = 10))
    
    evt <- next_event(ptr = q)
    expect_equal(evt$patient_id, 10L)
    expect_equal(evt$event_name, "death")
    expect_equal(evt$time, 11)
    
    pop_event(ptr = q)
    expect_true(queue_empty(q))
  })
  
  test_that("Multiple events for one patient work", {
    q <- queue_create(c("death", "dropout", "visit"))
    new_event2(c(death = 15, dropout = 8, visit = 5), ptr = q, patient_id = 10)
    
    expect_equal(queue_size(q), 3L)
    expect_equal(next_event(ptr = q)$event_name, "visit")
    pop_event(ptr = q)
    expect_equal(next_event(ptr = q)$event_name, "dropout")
    pop_event(ptr = q)
    expect_equal(next_event(ptr = q)$event_name, "death")
  })
  
  test_that("Priority ordering works for same time", {
    q <- queue_create(c("death", "dropout", "visit"))
    new_event2(c(death = 10, dropout = 10, visit = 10), ptr = q, patient_id = 10)
    
    expect_equal(next_event(ptr = q)$event_name, "death")
    pop_event(ptr = q)
    expect_equal(next_event(ptr = q)$event_name, "dropout")
    pop_event(ptr = q)
    expect_equal(next_event(ptr = q)$event_name, "visit")
  })
  
  test_that("Multiple patients work correctly", {
    q <- queue_create(c("death", "dropout", "visit"))
    new_event2(c(death = 20), ptr = q, patient_id = 1)
    new_event2(c(visit = 5), ptr = q, patient_id = 2)
    new_event2(c(dropout = 15), ptr = q, patient_id = 3)
    
    expect_equal(queue_size(q), 3L)
    
    evt <- next_event(ptr = q)
    expect_equal(evt$patient_id, 2L)
    expect_equal(evt$event_name, "visit")
  })
  
  test_that("Event modification works", {
    q <- queue_create(c("death", "dropout", "visit"))
    new_event2(c(death = 20, visit = 5), ptr = q, patient_id = 10)
    
    modify_event2(c(visit = 25), ptr = q, patient_id = 10)
    
    expect_equal(queue_size(q), 2L)
    evt <- next_event(ptr = q)
    expect_equal(evt$event_name, "death")
    pop_event(ptr = q)
    evt <- next_event(ptr = q)
    expect_equal(evt$event_name, "visit")
    expect_equal(evt$time, 25)
  })
  
  test_that("Create if missing works", {
    q <- queue_create(c("death", "dropout", "visit"))
    modify_event2(c(visit = 10), create_if_missing = TRUE, ptr = q, patient_id = 10)
    expect_equal(queue_size(q), 1L)
    expect_true(has_event("visit", ptr = q, patient_id = 10))
    
    modify_event2(c(death = 5), create_if_missing = FALSE, ptr = q, patient_id = 10)
    expect_false(has_event("death", ptr = q, patient_id = 10))
  })
  
  test_that("Event removal works", {
    q <- queue_create(c("death", "dropout", "visit"))
    new_event2(c(death = 15, dropout = 10, visit = 5), ptr = q, patient_id = 10)
    
    remove_event2(c("dropout", "visit"), ptr = q, patient_id = 10)
    
    expect_equal(queue_size(q), 1L)
    expect_true(has_event("death", ptr = q, patient_id = 10))
    expect_false(has_event("dropout", ptr = q, patient_id = 10))
    expect_false(has_event("visit", ptr = q, patient_id = 10))
  })
  
  test_that("Large scale operations work", {
    q <- queue_create(c("death", "dropout", "visit", "adverse_event"))
    for (i in 1:1000) {
      times <- runif(3, 0, 100)
      new_event2(c(death = times[1], dropout = times[2], visit = times[3]), ptr = q, patient_id = i)
    }
    
    expect_equal(queue_size(q), 3000L)
    
    for (i in 1:100) {
      modify_event2(c(visit = runif(1, 0, 50)), ptr = q, patient_id = i)
    }
    
    for (i in 1:50) {
      remove_event2("dropout", ptr = q, patient_id = i)
    }
    
    expect_equal(queue_size(q), 2950L)
  })
  
  test_that("Error handling works", {
    q <- queue_create(c("death", "dropout", "visit"))
    
    expect_error(new_event2(c(unknown_event = 5), ptr = q, patient_id = 10), "Unknown event type")
    expect_equal(next_event(ptr = q), list(patient_id = integer(0), event_name = character(0), time = numeric(0)))
    expect_error(pop_event(ptr = q), "Queue is empty")
    
    new_event2(c(death = 10), ptr = q, patient_id = 10)
    expect_true(has_event("death", ptr = q, patient_id = 10))
  })
  
  test_that("Memory management works under heavy modification", {
    q <- queue_create(c("death", "dropout", "visit"))
    
    for (i in 1:2000) {
      new_event2(c(visit = i), ptr = q, patient_id = 1L)
      for (j in 1:3) {
        modify_event2(c(visit = i + j * 0.1), ptr = q, patient_id = 1L)
      }
      if (i %% 10 == 0) {
        remove_event2("visit", ptr = q, patient_id = 1L)
        new_event2(c(visit = i + 100), ptr = q, patient_id = 1L)
      }
    }
    
    expect_equal(queue_size(q), 1L)
    evt <- next_event(ptr = q)
    expect_equal(evt$patient_id, 1L)
    expect_equal(evt$event_name, "visit")
  })
  
  test_that("Mixed operations work correctly", {
    q <- queue_create(c("death", "dropout", "visit", "lab_test"))
    new_event2(c(death = 100, visit = 10), ptr = q, patient_id = 1L)
    new_event2(c(dropout = 50, lab_test = 5), ptr = q, patient_id = 2L)
    
    expect_equal(queue_size(q), 4L)
    expect_equal(next_event(ptr = q)$event_name, "lab_test")
    pop_event(ptr = q)
    
    modify_event2(c(visit = 75), ptr = q, patient_id = 1L)
    expect_equal(next_event(ptr = q)$event_name, "dropout")
    pop_event(ptr = q)
    
    remove_event2("visit", ptr = q, patient_id = 1L)
    expect_equal(queue_size(q), 1L)
    expect_equal(next_event(ptr = q)$event_name, "death")
  })
  
  
  test_that("Use of cur_evtlist and i", {
    cur_evtlist <- queue_create(c("death", "dropout", "visit", "lab_test"))
    i <- 10
    new_event2(c(death = 100, visit = 10))
    new_event2(c(dropout = 50, lab_test = 5))
    
    expect_equal(queue_size(), 4L)
    expect_equal(length(next_event(10)$time), 4L)
    expect_equal(next_event()$event_name, "lab_test")
    pop_event()
    
    modify_event2(c(visit = 75))
    expect_equal(next_event()$event_name, "dropout")
    a <- next_event(10)
    expect_equal(a$time[a$event_name == "visit"], 75)
    pop_event()
    
    remove_event2("visit")
    expect_equal(queue_size(), 1L)
    expect_equal(next_event()$event_name, "death")
  })

