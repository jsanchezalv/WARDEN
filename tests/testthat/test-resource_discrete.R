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
})
  # Test Empty Resource (n=0)
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
  
  # Remove any existing i or curtime variables
  if (exists("i")) rm(i)
  if (exists("curtime")) rm(curtime)
  
  # Should throw error when i is missing
  expect_error(beds$attempt_block(), "'i' not found in parent frame")
  expect_error(beds$attempt_free(), "'i' not found in parent frame")
  expect_error(beds$attempt_free_if_using(), "'i' not found in parent frame")
  
  # Should throw error when curtime is missing
  expect_error(beds$attempt_block(patient_id = 101), "'curtime' not found in parent frame")
  expect_error(beds$remove_resource(1), "'curtime' not found in parent frame")
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