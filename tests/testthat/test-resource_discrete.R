library(testthat)

# Test Basic Functionality
test_that("basic functionality works correctly", {
  # Create resource
  beds <- resource_discrete(3)
  
  # Test initial state
  expect_equal(beds$size(), 3)
  expect_equal(beds$n_free(), 3)
  expect_equal(beds$queue_size(), 0)
  expect_equal(length(beds$patients_using()), 0)
  # Test Empty Resource (n=0)
  
})
  test_that("empty resource (n=0) works correctly", {
    empty_beds <- resource_discrete(0)
    
    # Test initial state
    expect_equal(empty_beds$size(), 0)
    expect_equal(empty_beds$n_free(), 0)
    expect_equal(empty_beds$queue_size(), 0)
    expect_equal(length(empty_beds$patients_using()), 0)
    
    # Any attempt to block should go to queue
    result <- empty_beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
    expect_false(result)
    expect_equal(empty_beds$queue_size(), 1)
    expect_equal(empty_beds$n_free(), 0)
    
    # Test status functions
    expect_false(empty_beds$is_patient_using(101))
    expect_true(empty_beds$is_patient_in_queue(101))
    
    # Test queue functions
    next_patient <- empty_beds$next_patient_in_line(1)
    expect_equal(next_patient[1], 101)
    
    queue_times <- empty_beds$queue_start_times()
    expect_equal(queue_times[1], 0.0)
    
    # Add resources to make it functional
    empty_beds$add_resource(2)
    expect_equal(empty_beds$size(), 2)
    expect_equal(empty_beds$n_free(), 2)
    
    # Now patient can get resource
    result2 <- empty_beds$attempt_block(patient_id = 101, priority = 1, start_time = 1.0)
    expect_true(result2)
    expect_equal(empty_beds$queue_size(), 0)
    expect_equal(empty_beds$n_free(), 1)
    expect_true(empty_beds$is_patient_using(101))
  })

# Test Blocking and Freeing
test_that("blocking and freeing work correctly", {
  beds <- resource_discrete(2)
  
  # Test blocking with explicit parameters
  result1 <- beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
  expect_true(result1)
  expect_equal(beds$n_free(), 1)
  expect_equal(beds$queue_size(), 0)
  expect_equal(beds$patients_using()[1], 101)
  
  # Block second resource
  result2 <- beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
  expect_true(result2)
  expect_equal(beds$n_free(), 0)
  expect_equal(length(beds$patients_using()), 2)
  
  # Try to block when full - should go to queue
  result3 <- beds$attempt_block(patient_id = 103, priority = 1, start_time = 2.0)
  expect_false(result3)
  expect_equal(beds$queue_size(), 1)
  
  # Free a resource
  beds$attempt_free(patient_id = 101)
  expect_equal(beds$n_free(), 1)
  expect_equal(length(beds$patients_using()), 1)
})

# Test Priority Queue
test_that("priority queue works correctly", {
  beds <- resource_discrete(1)
  
  # Block the only resource
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
  
  # Add patients to queue with different priorities
  beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)  # Low priority
  beds$attempt_block(patient_id = 103, priority = 3, start_time = 2.0)  # High priority
  beds$attempt_block(patient_id = 104, priority = 1, start_time = 3.0)  # Low priority
  
  expect_equal(beds$queue_size(), 3)
  
  # Check queue order - highest priority first
  next_patients <- beds$next_patient_in_line(3)
  expect_equal(next_patients[1], 103)  # Highest priority
  expect_equal(next_patients[2], 102)  # FIFO among same priority
  expect_equal(next_patients[3], 104)  # FIFO among same priority
  
  # Free resource - highest priority should get it
  beds$attempt_free(patient_id = 101)
  result <- beds$attempt_block(patient_id = 103, priority = 3, start_time = 4.0)
  expect_true(result)
  expect_equal(beds$queue_size(), 2)
})

# Test New Status Functions
test_that("status functions work correctly", {
  beds <- resource_discrete(2)
  
  # Block one resource
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
  expect_true(beds$is_patient_using(101))
  expect_false(beds$is_patient_in_queue(101))
  
  # Block second resource - fills capacity
  beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
  
  # Add to queue
  beds$attempt_block(patient_id = 103, priority = 1, start_time = 2.0)
  expect_true(beds$is_patient_in_queue(103))
  expect_false(beds$is_patient_using(103))
  
  # Test non-existent patient
  expect_false(beds$is_patient_using(999))
  expect_false(beds$is_patient_in_queue(999))
})

# Test Queue Start Times
test_that("queue start times work correctly", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)  # Uses resource
  
  # Add patients to queue at different times
  beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
  beds$attempt_block(patient_id = 103, priority = 2, start_time = 2.0)  # Higher priority
  
  queue_ids <- beds$next_patient_in_line(2)
  queue_times <- beds$queue_start_times()
  
  # Patient 103 should be first (higher priority) with start time 2.0
  expect_equal(queue_ids[1], 103)
  expect_equal(queue_times[1], 2.0)
  
  # Patient 102 should be second with start time 1.0
  expect_equal(queue_ids[2], 102)
  expect_equal(queue_times[2], 1.0)
})

# Test Priority Modification
test_that("priority modification works correctly", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
  
  # Add patients to queue
  beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
  beds$attempt_block(patient_id = 103, priority = 1, start_time = 2.0)
  
  # Modify priority of second patient
  beds$modify_priority(patient_id = 103, new_priority = 5)
  
  # Check that patient 103 is now first in line
  next_patient <- beds$next_patient_in_line(1)
  expect_equal(next_patient[1], 103)
  
  # Check that queue start time is preserved
  queue_times <- beds$queue_start_times()
  expect_equal(queue_times[1], 2.0)  # Original queue start time
})

# Test Resource Addition and Removal
test_that("resource addition and removal work correctly", {
  beds <- resource_discrete(2)
  
  # Fill resources
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
  beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
  
  # Add one more resource
  beds$add_resource(1)
  expect_equal(beds$size(), 3)
  expect_equal(beds$n_free(), 1)
  
  # Try to add patient 103 - should use the free resource directly (no queue needed)
  result <- beds$attempt_block(patient_id = 103, priority = 1, start_time = 2.0)
  expect_true(result)  # Should succeed and use the free resource
  expect_equal(beds$n_free(), 0)
  expect_equal(beds$queue_size(), 0)
  
  # Now add patient 104 - should go to queue since all resources are used
  result2 <- beds$attempt_block(patient_id = 104, priority = 1, start_time = 2.5)
  expect_false(result2)
  expect_equal(beds$queue_size(), 1)
  
  # Remove 2 resources - should move patient 102 and 103 to queue (keeping 101 using)
  beds$remove_resource(2, current_time = 5.0)
  expect_equal(beds$size(), 1)
  expect_equal(beds$n_free(), 0)  # Only patient 101 using
  expect_equal(beds$queue_size(), 3)  # 102, 103, and 104 in queue
  
  # Check that moved patients have new queue start time
  queue_ids <- beds$next_patient_in_line(3)
  queue_times <- beds$queue_start_times()
  
  # Find patients 102 and 103 in queue - they should have queue start time of 5.0 (when moved)
  patient_102_index <- which(queue_ids == 102)
  patient_103_index <- which(queue_ids == 103)
  
  expect_equal(queue_times[patient_102_index], 5.0)
  expect_equal(queue_times[patient_103_index], 5.0)
})

# Test Edge Cases
test_that("edge cases are handled correctly", {
  beds <- resource_discrete(2)
  
  # Test freeing non-existent patient (should not error)
  expect_silent(beds$attempt_free(patient_id = 999))
  
  # Test freeing from queue
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
  beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
  beds$attempt_block(patient_id = 103, priority = 1, start_time = 2.0)  # Goes to queue
  
  expect_equal(beds$queue_size(), 1)
  beds$attempt_free(patient_id = 103)  # Remove from queue
  expect_equal(beds$queue_size(), 0)
  
  # Test attempt_free_if_using
  beds$attempt_block(patient_id = 104, priority = 1, start_time = 3.0)  # Goes to queue
  beds$attempt_free_if_using(patient_id = 104)  # Should do nothing (not using)
  expect_equal(beds$queue_size(), 1)
  
  beds$attempt_free_if_using(patient_id = 101)  # Should free (is using)
  expect_equal(beds$n_free(), 1)
  
  # Test error cases - provide current_time explicitly
  expect_error(beds$remove_resource(10, current_time = 5.0), "Cannot remove more resources than available")
})

# Test Input Validation
test_that("input validation works correctly", {
  beds <- resource_discrete(2)
  
  # Test invalid resource creation
  expect_error(resource_discrete(-1), "n must be a single integer >= 0")
  expect_error(resource_discrete(c(1, 2)), "n must be a single integer >= 0")
  expect_error(resource_discrete("invalid"), "n must be a single integer >= 0")
  
  # Test valid edge case: n = 0
  empty_resource <- resource_discrete(0)
  expect_equal(empty_resource$size(), 0)
  expect_equal(empty_resource$n_free(), 0)
  expect_equal(empty_resource$queue_size(), 0)
  
  # Test invalid patient_id inputs
  expect_error(beds$is_patient_using(c(1, 2)), "patient_id must be a single number")
  expect_error(beds$is_patient_in_queue("invalid"), "patient_id must be a single number")
  
  # Test invalid priority inputs
  expect_error(beds$attempt_block(patient_id = 101, priority = c(1, 2), start_time = 0.0), 
               "priority must be a single number")
  expect_error(beds$attempt_block(patient_id = 101, priority = "invalid", start_time = 0.0), 
               "priority must be a single number")
  
  # Test invalid start_time inputs
  expect_error(beds$attempt_block(patient_id = 101, priority = 1, start_time = c(1, 2)), 
               "start_time must be a single number")
  expect_error(beds$attempt_block(patient_id = 101, priority = 1, start_time = "invalid"), 
               "start_time must be a single number")
  
  # Test invalid remove_all inputs
  expect_error(beds$attempt_free(patient_id = 101, remove_all = "invalid"), 
               "remove_all must be a single logical value")
  expect_error(beds$attempt_free(patient_id = 101, remove_all = c(TRUE, FALSE)), 
               "remove_all must be a single logical value")
  
  # Test invalid n inputs for next_patient_in_line
  expect_error(beds$next_patient_in_line(0), "n must be a single positive integer")
  expect_error(beds$next_patient_in_line(-1), "n must be a single positive integer")
  expect_error(beds$next_patient_in_line(c(1, 2)), "n must be a single positive integer")
  
  # Test invalid new_priority inputs
  expect_error(beds$modify_priority(101, c(1, 2)), "new_priority must be a single number")
  expect_error(beds$modify_priority(101, "invalid"), "new_priority must be a single number")
  
  # Test invalid n_to_add inputs
  expect_error(beds$add_resource(0), "n_to_add must be a single positive integer")
  expect_error(beds$add_resource(-1), "n_to_add must be a single positive integer")
  expect_error(beds$add_resource(c(1, 2)), "n_to_add must be a single positive integer")
  
  # Test invalid n_to_remove inputs
  expect_error(beds$remove_resource(0, current_time = 1.0), "n_to_remove must be a single positive integer")
  expect_error(beds$remove_resource(-1, current_time = 1.0), "n_to_remove must be a single positive integer")
  expect_error(beds$remove_resource(c(1, 2), current_time = 1.0), "n_to_remove must be a single positive integer")
  
  # Test invalid current_time inputs
  expect_error(beds$remove_resource(1, current_time = c(1, 2)), "current_time must be a single number")
  expect_error(beds$remove_resource(1, current_time = "invalid"), "current_time must be a single number")
})

# Test Default Parameter Behavior
test_that("default parameter behavior works correctly", {
  beds <- resource_discrete(2)
  
  # Test with variables in environment
  i <- 101
  curtime <- 5.0
  
  # Should use default i and curtime
  result <- beds$attempt_block()
  expect_true(result)
  expect_true(beds$is_patient_using(101))
  expect_equal(beds$patients_using_times()[1], 5.0)
  
  # Test attempt_free with default i
  i <- 101
  beds$attempt_free()
  expect_false(beds$is_patient_using(101))
  
  # Test remove_resource with default curtime
  beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
  beds$attempt_block(patient_id = 103, priority = 1, start_time = 2.0)
  
  curtime <- 10.0
  beds$remove_resource(1)  # Should use curtime = 10.0
  
  # Check that moved patient has curtime as queue start time
  queue_times <- beds$queue_start_times()
  expect_true(any(queue_times == 10.0))
})

# Test Missing Default Variables
test_that("missing default variables throw appropriate errors", {
  beds <- resource_discrete(2)
  
  # Create helper functions that run in clean environments
  test_missing_i <- function() {
    # This function has no 'i' or 'curtime' variables
    beds$attempt_block()
  }
  
  test_missing_curtime <- function() {
    # This function has 'i' but no 'curtime'
    i <- 101
    beds$attempt_block()
  }
  
  test_missing_curtime_remove <- function() {
    # This function has no 'curtime' for remove_resource
    beds$remove_resource(1)
  }
  
  # Test when i is missing
  expect_error(test_missing_i())
  
  # Test when curtime is missing (but i exists)
  expect_error(test_missing_curtime())
  expect_error(test_missing_curtime_remove())
  
  # Test other functions that need 'i'
  test_missing_i_free <- function() {
    if(exists("i")){rm(i)}
    beds$attempt_free()
  }
  
  test_missing_i_free_if_using <- function() {
    beds$attempt_free_if_using()
  }
  
  expect_error(test_missing_i_free())
  expect_error(test_missing_i_free_if_using())
})

# Test Lazy Deletion Cleanup
test_that("lazy deletion cleanup works correctly", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
  
  # Add many patients and modify priorities to trigger cleanup
  for (i in 1:10) {
    beds$attempt_block(patient_id = 100 + i, priority = 1, start_time = i)
  }
  
  # Modify priorities multiple times to create invalid entries
  for (i in 1:5) {
    beds$modify_priority(patient_id = 102, new_priority = i)
    beds$modify_priority(patient_id = 103, new_priority = i + 1)
  }
  
  # Check that queue still works correctly
  next_patients <- beds$next_patient_in_line(3)
  expect_equal(length(next_patients), 3)
  
  # Verify patients are still tracked correctly
  expect_true(beds$is_patient_in_queue(102))
  expect_true(beds$is_patient_in_queue(103))
})

# Test Remove All Functionality
test_that("remove_all functionality works correctly", {
  beds <- resource_discrete(3)
  
  # Block same patient multiple times (if allowed by your logic)
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 1.0)
  beds$attempt_block(patient_id = 102, priority = 1, start_time = 2.0)
  
  # Check initial state
  expect_equal(beds$n_free(), 0)
  expect_equal(length(beds$patients_using()), 3)
  
  # Free all instances of patient 101
  beds$attempt_free(patient_id = 101, remove_all = TRUE)
  
  # Should have freed 2 instances
  expect_equal(beds$n_free(), 2)
  expect_equal(length(beds$patients_using()), 1)
  expect_equal(beds$patients_using()[1], 102)
})

# Test Print Method
test_that("print method works correctly", {
  beds <- resource_discrete(3)
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
  beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
  
  # Fill all resources first, then add patient 103 to queue
  beds$attempt_block(patient_id = 103, priority = 1, start_time = 2.0)  # Uses last resource
  beds$attempt_block(patient_id = 104, priority = 2, start_time = 3.0)  # Goes to queue
  
  # Verify state before testing print
  expect_equal(beds$n_free(), 0)
  expect_equal(beds$queue_size(), 1)
  expect_equal(length(beds$patients_using()), 3)
  
  # Capture print output
  output <- capture.output(print(beds))
  
  expect_true(any(grepl("Discrete Resource:", output)))
  expect_true(any(grepl("Total capacity: 3", output)))
  expect_true(any(grepl("Free units: 0", output)))
  expect_true(any(grepl("Queue size: 1", output)))
  expect_true(any(grepl("Patients using: 3", output)))
})

# ── New features ─────────────────────────────────────────────────────────────

# F1.1 Constructor parameters
test_that("discipline and max_queue constructor parameters work", {
  expect_silent(resource_discrete(5, discipline = "FIFO"))
  expect_silent(resource_discrete(5, discipline = "LIFO"))
  expect_silent(resource_discrete(5, max_queue = 3))
  expect_silent(resource_discrete(5, discipline = "LIFO", max_queue = 2))
  expect_error(resource_discrete(5, discipline = "INVALID"))
  expect_error(resource_discrete(5, max_queue = -1))
})

# F1.2 LIFO discipline
test_that("LIFO serves last-queued patient first among same priority", {
  beds <- resource_discrete(1, discipline = "LIFO")
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 3)
  expect_equal(beds$next_patient_in_line(1L), 3L)
})

test_that("LIFO priority still takes precedence over insertion order", {
  beds <- resource_discrete(1, discipline = "LIFO")
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 3L, priority = 2L, start_time = 2)
  expect_equal(beds$next_patient_in_line(1L), 3L)
})

test_that("FIFO serves first-queued patient first", {
  beds <- resource_discrete(1, discipline = "FIFO")
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 3)
  expect_equal(beds$next_patient_in_line(1L), 2L)
})

# F1.3 Limited queue (max_queue)
test_that("max_queue rejects when queue is full", {
  beds <- resource_discrete(1, max_queue = 2)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  expect_false(beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1))
  expect_false(beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 2))
  result <- beds$attempt_block(patient_id = 4L, priority = 1L, start_time = 3)
  expect_true(is.na(result))
})

test_that("max_queue=0 rejects immediately when capacity full", {
  beds <- resource_discrete(1, max_queue = 0)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  result <- beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  expect_true(is.na(result))
  expect_equal(beds$queue_size(), 0L)
})

test_that("max_queue allows queuing after a queue slot opens", {
  beds <- resource_discrete(1, max_queue = 2)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 2)
  beds$attempt_free(patient_id = 2L)
  result <- beds$attempt_block(patient_id = 4L, priority = 1L, start_time = 3)
  expect_false(result)
})

test_that("unlimited queue (default) never rejects", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  for (j in 2:101) {
    expect_false(beds$attempt_block(patient_id = j, priority = 1L, start_time = j))
  }
})

# F1.4 Batch amount (single patient, multiple units)
test_that("attempt_block with amount blocks multiple units", {
  beds <- resource_discrete(5)
  result <- beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  expect_true(result)
  expect_equal(beds$n_free(), 2L)
})

test_that("attempt_block with amount queues when insufficient capacity", {
  beds <- resource_discrete(2)
  result <- beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  expect_false(result)
})

test_that("n_free correct after patients with different amounts", {
  beds <- resource_discrete(10)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1, amount = 4L)
  expect_equal(beds$n_free(), 3L)
})

test_that("attempt_free frees correct amount of units", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  beds$attempt_free(patient_id = 1L, amount = 3L)
  expect_equal(beds$n_free(), 5L)
  expect_false(beds$is_patient_using(1L))
})

test_that("queued patient with amount dequeues when enough slots available", {
  beds <- resource_discrete(3)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1, amount = 2L)
  beds$attempt_free(patient_id = 1L, amount = 3L)
  result <- beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2, amount = 2L)
  expect_true(result)
  expect_equal(beds$n_free(), 1L)
})

# F1.5 Queue wait time
test_that("queue_wait_time returns NA for patient that never queued", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  expect_equal(beds$queue_wait_time(1L), NA_real_)
})

test_that("queue_wait_time_current returns elapsed wait while in queue", {
  beds <- resource_discrete(1)
  b <- resource_discrete(1)
  i <- 1L
  curtime <- 5.0
  seize_all(list(beds, b))
  i <- 2L
  seize_all(list(beds, b))   # P2 queued for beds at t=5
  curtime <- 10.0
  # P2 still in queue: elapsed = 10 - 5 = 5
  expect_equal(beds$queue_wait_time_current(2L, 10.0), 5.0)
  # queue_wait_time (final-only) still NA while in queue
  expect_equal(beds$queue_wait_time(2L), NA_real_)
})

test_that("queue_wait_time returns NA while patient is still in queue", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 5)
  expect_equal(beds$queue_wait_time(2L), NA_real_)
})

test_that("queue_wait_time returns most recent wait for repeated queuing", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  beds$attempt_free(patient_id = 1L)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 7)  # wait=5
  beds$attempt_free(patient_id = 2L)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 7)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 10)
  beds$attempt_free(patient_id = 1L)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 12)  # wait=2
  expect_equal(beds$queue_wait_time(2L), 2.0)
})

# F1.6 had_to_queue
test_that("had_to_queue returns 0L for direct acquisition", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  expect_equal(beds$had_to_queue(1L), 0L)
})

test_that("had_to_queue returns 1L for patient that queued", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  expect_equal(beds$had_to_queue(2L), 1L)
})

test_that("had_to_queue stays 1L even after patient acquires", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_free(patient_id = 1L)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  expect_equal(beds$had_to_queue(2L), 1L)
})

test_that("had_to_queue returns 0L for patient that never interacted", {
  beds <- resource_discrete(5)
  expect_equal(beds$had_to_queue(99L), 0L)
})

# F1.7 time_in_use
test_that("time_in_use returns correct duration", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 3.0)
  expect_equal(beds$time_in_use(patient_id = 1L, current_time = 8.0), 5.0)
})

test_that("time_in_use returns NA for patient not currently using", {
  beds <- resource_discrete(5)
  expect_equal(beds$time_in_use(patient_id = 1L, current_time = 5.0), NA_real_)
})

# F1.8 Aggregate statistics
test_that("total_patients_blocked counts unique patients not occurrences", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_free(patient_id = 1L)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 1)
  expect_equal(beds$total_patients_blocked(), 1L)
})

test_that("total_patients_blocked counts multiple distinct patients", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 2)
  expect_equal(beds$total_patients_blocked(), 3L)
})

test_that("total_patients_queued counts unique queued patients", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 2)
  expect_equal(beds$total_patients_queued(), 2L)
})

test_that("total_patients_queued is 0 when all acquired directly", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 2)
  expect_equal(beds$total_patients_queued(), 0L)
})

# F1.9 utilization and n_using
test_that("utilization is 0 on fresh resource", {
  beds <- resource_discrete(5)
  expect_equal(beds$utilization(), 0)
})

test_that("utilization and n_using correct for partial use", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 2)
  expect_equal(beds$utilization(), 0.6)
  expect_equal(beds$n_using(), 3L)
})

test_that("utilization is 1.0 when fully occupied", {
  beds <- resource_discrete(3)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 2)
  expect_equal(beds$utilization(), 1.0)
})

test_that("utilization accounts for batch amounts", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  expect_equal(beds$utilization(), 0.6)
})

# F1.10 batch_seize
test_that("batch_seize acquires first patients and queues rest", {
  beds <- resource_discrete(3)
  result <- beds$batch_seize(patient_ids = 1:5, priority = 1L, start_time = 0)
  expect_equal(result[1:3], c(1L, 1L, 1L))
  expect_equal(result[4:5], c(0L, 0L))
})

test_that("batch_seize rejects when max_queue reached", {
  beds <- resource_discrete(3, max_queue = 1)
  result <- beds$batch_seize(patient_ids = 1:5, priority = 1L, start_time = 0)
  expect_equal(result[1:3], c(1L, 1L, 1L))
  expect_equal(result[4], 0L)
  expect_equal(result[5], -1L)
})

test_that("batch_seize returns integer vector of correct length", {
  beds <- resource_discrete(5)
  result <- beds$batch_seize(patient_ids = 1:10, priority = 1L, start_time = 0)
  expect_type(result, "integer")
  expect_length(result, 10L)
})

test_that("batch_seize with amount_each uses correct units per patient", {
  beds <- resource_discrete(5)
  result <- beds$batch_seize(patient_ids = 1:3, priority = 1L, start_time = 0, amount_each = 2L)
  expect_equal(result[1:2], c(1L, 1L))
  expect_equal(result[3], 0L)
})

# F1.11 no-arg is_patient_using / is_patient_in_queue
test_that("is_patient_using no-arg uses i from parent frame", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  i <- 1L
  expect_true(beds$is_patient_using())
  i <- 99L
  expect_false(beds$is_patient_using())
})

test_that("is_patient_in_queue no-arg uses i from parent frame", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  i <- 2L
  expect_true(beds$is_patient_in_queue())
  i <- 1L
  expect_false(beds$is_patient_in_queue())
})

# ── End new features ──────────────────────────────────────────────────────────

# Test Complex Scenario
test_that("complex scenario works correctly", {
  beds <- resource_discrete(2)
  
  # Fill resources
  beds$attempt_block(patient_id = 101, priority = 1, start_time = 0.0)
  beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
  
  # Add to queue with different priorities
  beds$attempt_block(patient_id = 103, priority = 3, start_time = 2.0)  # High priority
  beds$attempt_block(patient_id = 104, priority = 1, start_time = 3.0)  # Low priority
  beds$attempt_block(patient_id = 105, priority = 2, start_time = 4.0)  # Medium priority
  
  # Verify queue order
  queue_order <- beds$next_patient_in_line(3)
  expect_equal(queue_order[1], 103)  # Highest priority
  expect_equal(queue_order[2], 105)  # Medium priority
  expect_equal(queue_order[3], 104)  # Lowest priority
  
  # Modify priority to change order
  beds$modify_priority(patient_id = 104, new_priority = 5)  # Make highest priority
  
  # Check new order
  new_queue_order <- beds$next_patient_in_line(3)
  expect_equal(new_queue_order[1], 104)  # Now highest priority
  
  # Free a resource - patient 104 should get it
  beds$attempt_free(patient_id = 101)
  result <- beds$attempt_block(patient_id = 104, priority = 5, start_time = 5.0)
  expect_true(result)
  expect_equal(beds$queue_size(), 2)  # 103 and 105 left in queue
  
  # Add more resources and verify they don't automatically assign
  beds$add_resource(2)
  expect_equal(beds$size(), 4)
  expect_equal(beds$n_free(), 2)  # 2 new free resources
  expect_equal(beds$queue_size(), 2)  # Queue unchanged
  
  # Remove more resources than capacity allows - should error
  expect_error(beds$remove_resource(5, current_time = 6.0),
               "Cannot remove more resources than available")
})

# ── Amount validation ─────────────────────────────────────────────────────────

test_that("attempt_block amount validation", {
  beds <- resource_discrete(5)
  expect_error(beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 0L))
  expect_error(beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = -1L))
  expect_error(beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 2.5))
})

test_that("indivisible release: amount mismatch errors", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  expect_error(beds$attempt_free(patient_id = 1L, amount = 1L))
  beds2 <- resource_discrete(5)
  beds2$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  expect_silent(beds2$attempt_free(patient_id = 1L, amount = 3L))
  expect_equal(beds2$n_free(), 5L)
})

test_that("release(amount = NULL) means release(amount = 1)", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  expect_error(beds$attempt_free(patient_id = 1L, amount = NULL))
  beds2 <- resource_discrete(5)
  beds2$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 1L)
  expect_silent(beds2$attempt_free(patient_id = 1L, amount = NULL))
  expect_equal(beds2$n_free(), 5L)
})

test_that("remove_all releases all entries regardless of amount", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 2L)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 1, amount = 3L)
  expect_equal(beds$n_free(), 0L)
  beds$attempt_free(patient_id = 1L, remove_all = TRUE)
  expect_equal(beds$n_free(), 5L)
  expect_false(beds$is_patient_using(1L))
})

test_that("multiple seize entries are independent", {
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 2L)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 1, amount = 3L)
  expect_equal(beds$n_free(), 0L)
  beds$attempt_free(patient_id = 1L, amount = 2L)
  expect_equal(beds$n_free(), 2L)
  beds$attempt_free(patient_id = 1L, amount = 3L)
  expect_equal(beds$n_free(), 5L)
})

# ── allow_multiple_queue ──────────────────────────────────────────────────────

test_that("allow_multiple_queue = FALSE rejects repeated queue", {
  beds <- resource_discrete(1, allow_multiple_queue = FALSE)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  result <- beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  expect_true(is.na(result))
  expect_equal(beds$queue_size(), 1L)
})

test_that("allow_multiple_queue = FALSE allows acquisition if first in line", {
  beds <- resource_discrete(1, allow_multiple_queue = FALSE)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_free(patient_id = 1L)
  result <- beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  expect_true(result)
  expect_equal(beds$queue_size(), 0L)
})

test_that("allow_multiple_queue = TRUE allows repeated queuing", {
  beds <- resource_discrete(1, allow_multiple_queue = TRUE)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  expect_equal(beds$queue_size(), 2L)
})

# ── seize_all all_or_none ─────────────────────────────────────────────────────

test_that("seize_all all_or_none queues ALL bottlenecks", {
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))
  i <- 2L; curtime <- 1.0
  result <- seize_all(list(beds, chairs))
  expect_false(result)
  expect_equal(beds$queue_size(), 1L)
  expect_equal(chairs$queue_size(), 1L)
})

test_that("seize_all rejects duplicate resources", {
  beds <- resource_discrete(5)
  i <- 1L; curtime <- 0.0
  expect_error(seize_all(list(beds, beds)), "duplicate")
})

test_that("seize_all force_unblock = FALSE returns FALSE on deadlock", {
  # True deadlock setup: P2 pre-queued on both, capacity now free, P3 still first on chairs
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))    # P1 acquires both
  i <- 3L; curtime <- 0.5
  seize_all(list(chairs))          # P3 queues for chairs (P3 is first in chairs queue)
  i <- 2L; curtime <- 1.0
  seize_all(list(beds, chairs))    # P2 queues on both: first on beds, second on chairs
  i <- 1L; curtime <- 2.0
  beds$attempt_free(patient_id = 1L)     # P1 releases beds (P2 now first, capacity free)
  chairs$attempt_free(patient_id = 1L)   # P1 releases chairs (P3 still first, P2 second)
  i <- 2L; curtime <- 3.0
  result <- seize_all(list(beds, chairs), force_unblock = FALSE)
  expect_false(result)             # deadlock: P3 ahead on chairs, but force_unblock=FALSE
  expect_false(beds$is_patient_using(2L))    # P2 did NOT acquire beds
  expect_false(chairs$is_patient_using(2L))  # P2 did NOT acquire chairs
})

test_that("seize_all force_unblock = TRUE acquires on deadlock", {
  # Same true deadlock setup: P2 pre-queued on both, P3 first on chairs
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))    # P1 acquires both
  i <- 3L; curtime <- 0.5
  seize_all(list(chairs))          # P3 queues for chairs (P3 is first)
  i <- 2L; curtime <- 1.0
  seize_all(list(beds, chairs))    # P2 queues on both: first on beds, second on chairs
  i <- 1L; curtime <- 2.0
  beds$attempt_free(patient_id = 1L)     # P1 releases beds
  chairs$attempt_free(patient_id = 1L)   # P1 releases chairs (P3 still first, P2 second)
  i <- 2L; curtime <- 3.0
  result <- seize_all(list(beds, chairs), force_unblock = TRUE)
  expect_true(result)              # P2 jumps to front of chairs, acquires both
  expect_true(beds$is_patient_using(2L))
  expect_true(chairs$is_patient_using(2L))
})

test_that("seize_all force_unblock has no effect when capacity insufficient", {
  # Both resources full — no deadlock, just insufficient capacity
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))    # P1 acquires both (both full, no free capacity)
  i <- 3L; curtime <- 0.5
  seize_all(list(beds))            # P3 queues for beds
  i <- 2L; curtime <- 1.0
  result <- seize_all(list(beds, chairs), force_unblock = TRUE)
  expect_false(result)  # capacity insufficient: force_unblock has no effect
})

# ── release_all all_or_none ───────────────────────────────────────────────────

test_that("release_all only acts if using ALL resources", {
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))
  i <- 2L; curtime <- 1.0
  seize_all(list(beds, chairs))   # P2 queues on both (both full)
  i <- 1L; curtime <- 2.0
  beds$attempt_free(patient_id = 1L)   # P1 releases beds; P2 in queue but not auto-acquired
  expect_equal(beds$n_free(), 1L)      # capacity freed, P2 still in queue
  expect_equal(chairs$n_free(), 0L)    # chairs still held by P1
  # P1 is no longer using beds → all_or_none policy → message + no action
  expect_message(release_all(list(beds, chairs)))
  expect_equal(chairs$n_free(), 0L)    # chairs NOT released (policy blocked it)
})

test_that("release_all with NULL amounts purges queue entries", {
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))
  i <- 2L; curtime <- 1.0
  seize_all(list(beds, chairs))
  i <- 1L
  release_all(list(beds, chairs), amounts = NULL)
  expect_equal(beds$n_free(), 1L)
  expect_equal(chairs$n_free(), 1L)
  expect_false(beds$is_patient_in_queue(1L))
  expect_false(chairs$is_patient_in_queue(1L))
})

test_that("release_all with amounts does NOT purge queue entries", {
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))
  i <- 2L; curtime <- 1.0
  seize_all(list(beds, chairs))
  i <- 1L
  release_all(list(beds, chairs), amounts = c(1L, 1L))
  expect_equal(beds$n_free(), 1L)
  expect_equal(chairs$n_free(), 1L)
  expect_true(beds$is_patient_in_queue(2L))
  expect_true(chairs$is_patient_in_queue(2L))
})

test_that("release_all errors when amounts do not match seized amounts", {
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))           # P1 seizes both with amount=1
  expect_error(
    release_all(list(beds, chairs), amounts = c(3L, 1L)),  # amount=3 ≠ seized=1 for beds
    "amount mismatch"
  )
  expect_equal(beds$n_free(), 0L)         # capacity NOT freed (error occurred first)
})

test_that("release_all purges queue entries when patient not using resources", {
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))            # P1 acquires both

  i <- 2L; curtime <- 1.0
  acquired <- seize_all(list(beds, chairs)) # P2 queues on both
  expect_false(acquired)
  expect_equal(beds$queue_size(),   1L)
  expect_equal(chairs$queue_size(), 1L)

  # P2 "dies" while queued: release_all(amounts = NULL) should purge queues
  release_all(list(beds, chairs))

  expect_equal(beds$queue_size(),   0L)
  expect_equal(chairs$queue_size(), 0L)

  # P1 releases: queue is empty so no dead patient rescheduled
  i <- 1L
  release_all(list(beds, chairs))
  expect_equal(beds$n_free(),   1L)
  expect_equal(chairs$n_free(), 1L)
})

# ── release_if_using ──────────────────────────────────────────────────────────

test_that("release_if_using does nothing when not using", {
  beds <- resource_discrete(3)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  i <- 2L
  expect_silent(release_if_using(beds))
  expect_equal(beds$n_free(), 2L)
})

test_that("release_if_using releases when using", {
  beds <- resource_discrete(3)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  i <- 1L
  release_if_using(beds)
  expect_equal(beds$n_free(), 3L)
  expect_false(beds$is_patient_using(1L))
})

test_that("release_if_using never touches queue entries", {
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  i <- 2L
  release_if_using(beds)
  expect_true(beds$is_patient_in_queue(2L))
  expect_equal(beds$queue_size(), 1L)
})

# ── Queued immutable amount ───────────────────────────────────────────────────

test_that("queued amount is immutable: retry with different amount does not dequeue", {
  # P1 holds 3 of 4 capacity; n_free=1; P2 queues for 3. Retry with amount=1:
  # n_free(1) < queued_amount(3) so P2 cannot dequeue — queues again
  beds <- resource_discrete(4)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1, amount = 3L)
  result <- beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2, amount = 1L)
  expect_false(result)
  expect_equal(beds$queue_size(), 2L)            # original entry(3) + new entry(1)
})

test_that("queued amount is immutable: allow_multiple_queue=TRUE adds second entry on retry", {
  beds <- resource_discrete(4, allow_multiple_queue = TRUE)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1, amount = 3L)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2, amount = 1L)
  expect_equal(beds$queue_size(), 2L)
  expect_true(beds$is_patient_in_queue(2L))
})

test_that("queued amount is immutable: allow_multiple_queue=FALSE rejects retry", {
  beds <- resource_discrete(4, allow_multiple_queue = FALSE)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1, amount = 3L)
  result <- beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2, amount = 1L)
  expect_true(is.na(result))
  expect_equal(beds$queue_size(), 1L)
})

test_that("queued amount is honored on dequeue: patient acquires with original amount", {
  beds <- resource_discrete(3)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 3L)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1, amount = 3L)
  beds$attempt_free(patient_id = 1L, amount = 3L)  # frees 3 units
  # Now beds n_free=3, P2 queued for 3, P2 is first in queue
  result <- beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2, amount = 1L)
  expect_true(result)          # dequeue-and-acquire using original queued amount (3)
  expect_equal(beds$n_free(), 0L)  # all 3 units consumed by P2
})

# ── Multiple queue entries ────────────────────────────────────────────────────

test_that("allow_multiple_queue=TRUE: queue_size reflects multiple entries per patient", {
  beds <- resource_discrete(1, allow_multiple_queue = TRUE)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  expect_equal(beds$queue_size(), 2L)
})

test_that("next_patient_in_line returns same patient twice when queued twice", {
  beds <- resource_discrete(1, allow_multiple_queue = TRUE)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  result <- beds$next_patient_in_line(2L)
  expect_equal(result[1], 2L)
  expect_equal(result[2], 2L)
})

test_that("attempt_free(remove_all=FALSE) removes one queue entry per call", {
  beds <- resource_discrete(1, allow_multiple_queue = TRUE)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  expect_equal(beds$queue_size(), 2L)
  beds$attempt_free(patient_id = 2L)
  expect_equal(beds$queue_size(), 1L)
  beds$attempt_free(patient_id = 2L)
  expect_equal(beds$queue_size(), 0L)
})

test_that("attempt_free(remove_all=TRUE) removes all queue entries at once", {
  beds <- resource_discrete(1, allow_multiple_queue = TRUE)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 2)
  expect_equal(beds$queue_size(), 2L)
  beds$attempt_free(patient_id = 2L, remove_all = TRUE)
  expect_equal(beds$queue_size(), 0L)
  expect_false(beds$is_patient_in_queue(2L))
})

test_that("modify_priority updates all active queue entries, start times preserved", {
  beds <- resource_discrete(1, allow_multiple_queue = TRUE)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 2)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 3)  # P2 has 2 entries
  beds$modify_priority(patient_id = 2L, new_priority = 5L)             # elevate P2
  result <- beds$next_patient_in_line(2L)
  expect_equal(result[1], 2L)   # P2 is now first (highest priority)
  expect_equal(result[2], 2L)   # both of P2's entries are first
  times <- beds$queue_start_times()
  expect_true(all(c(1.0, 3.0) %in% times[1:2]))  # original start times preserved
})

# ── Fix B: accum_queue parameter ─────────────────────────────────────────────

test_that("seize_all accum_queue=FALSE prevents duplicate queue entries", {
  beds   <- resource_discrete(1, allow_multiple_queue = TRUE)
  chairs <- resource_discrete(1, allow_multiple_queue = TRUE)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))
  i <- 2L; curtime <- 1.0
  seize_all(list(beds, chairs), accum_queue = FALSE)
  expect_equal(beds$queue_size(),   1L)
  expect_equal(chairs$queue_size(), 1L)
  curtime <- 2.0
  result <- seize_all(list(beds, chairs), accum_queue = FALSE)
  expect_false(result)
  expect_equal(beds$queue_size(),   1L)
  expect_equal(chairs$queue_size(), 1L)
})

test_that("seize_all accum_queue=TRUE (default) adds entries on retry", {
  beds   <- resource_discrete(1, allow_multiple_queue = TRUE)
  chairs <- resource_discrete(1, allow_multiple_queue = TRUE)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))
  i <- 2L; curtime <- 1.0
  seize_all(list(beds, chairs))
  curtime <- 2.0
  seize_all(list(beds, chairs))
  expect_equal(beds$queue_size(),   2L)
  expect_equal(chairs$queue_size(), 2L)
})

test_that("seize_all accum_queue=FALSE still queues on first attempt", {
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))
  i <- 2L; curtime <- 1.0
  result <- seize_all(list(beds, chairs), accum_queue = FALSE)
  expect_false(result)
  expect_equal(beds$queue_size(),   1L)
  expect_equal(chairs$queue_size(), 1L)
})

# ── Fix C: allow_multiple_queue=FALSE + seize_all ────────────────────────────

test_that("seize_all does not reject when allow_multiple_queue=FALSE and already queued", {
  beds   <- resource_discrete(1, allow_multiple_queue = FALSE)
  chairs <- resource_discrete(1, allow_multiple_queue = FALSE)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))
  i <- 2L; curtime <- 1.0
  result1 <- seize_all(list(beds, chairs))
  expect_false(result1)
  expect_equal(beds$queue_size(),   1L)
  expect_equal(chairs$queue_size(), 1L)
  curtime <- 2.0
  result2 <- seize_all(list(beds, chairs))
  expect_false(result2)
  expect_equal(beds$queue_size(),   1L)
  expect_equal(chairs$queue_size(), 1L)
})

test_that("seize_all allow_multiple_queue=FALSE acquires when resources free up", {
  beds   <- resource_discrete(1, allow_multiple_queue = FALSE)
  chairs <- resource_discrete(1, allow_multiple_queue = FALSE)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))
  i <- 2L; curtime <- 1.0
  seize_all(list(beds, chairs))
  i <- 1L
  beds$attempt_free(patient_id = 1L)
  chairs$attempt_free(patient_id = 1L)
  i <- 2L; curtime <- 2.0
  result <- seize_all(list(beds, chairs))
  expect_true(result)
  expect_true(beds$is_patient_using(2L))
  expect_true(chairs$is_patient_using(2L))
})

test_that("seize_all rejects when not queued and queue full", {
  beds   <- resource_discrete(1, allow_multiple_queue = FALSE, max_queue = 1)
  chairs <- resource_discrete(1, allow_multiple_queue = FALSE, max_queue = 1)
  i <- 1L; curtime <- 0.0
  seize_all(list(beds, chairs))
  i <- 2L; curtime <- 1.0
  seize_all(list(beds, chairs))
  i <- 3L; curtime <- 2.0
  result <- seize_all(list(beds, chairs))
  expect_true(is.na(result))
})

# ── Fix E: priority inflation headroom ───────────────────────────────────────

test_that("force_unblock works correctly after many priority bumps", {
  beds   <- resource_discrete(1)
  chairs <- resource_discrete(1)
  for (iter in 1:100) {
    i <- 1L; curtime <- 0.0
    beds2   <- resource_discrete(1)
    chairs2 <- resource_discrete(1)
    beds2$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
    chairs2$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
    beds2$attempt_block(patient_id = 3L, priority = 1L, start_time = 0.5)
    i <- 2L; curtime <- 1.0
    seize_all(list(beds2, chairs2))
    i <- 1L
    beds2$attempt_free(patient_id = 1L)
    chairs2$attempt_free(patient_id = 1L)
    i <- 2L; curtime <- 2.0
    result <- seize_all(list(beds2, chairs2), force_unblock = TRUE)
    expect_true(result)
  }
})

# Item 1: seize() errors on seize_all()-specific arguments
test_that("seize errors when accum_queue is passed", {
  res <- resource_discrete(2)
  i <- 1L; curtime <- 0.0
  expect_error(seize(res, accum_queue = FALSE), "seize_all")
})

test_that("seize errors when force_unblock is passed", {
  res <- resource_discrete(2)
  i <- 1L; curtime <- 0.0
  expect_error(seize(res, force_unblock = TRUE), "seize_all")
})

test_that("seize errors on unknown extra arguments", {
  res <- resource_discrete(2)
  i <- 1L; curtime <- 0.0
  expect_error(seize(res, foo = 1), "unexpected argument")
})

test_that("seize works with priority argument", {
  res <- resource_discrete(1)
  i <- 1L; curtime <- 0.0
  result <- seize(res, priority = 2L)
  expect_true(result)
})

test_that("seize still works with default arguments", {
  res <- resource_discrete(2)
  i <- 1L; curtime <- 0.0
  expect_true(seize(res))
  expect_true(seize(res, amount = 1L))
})