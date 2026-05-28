library(testthat)

# Helper tests need i, curtime, cur_evtlist in the calling scope.
# Setting these directly in test_that bodies works because parent.frame()
# inside seize/release resolves to the test_that evaluation environment.

# F2.1 seize() ----------------------------------------------------------------

test_that("seize returns TRUE when capacity available", {
  i <- 1L
  curtime <- 0.0
  beds <- resource_discrete(3)
  expect_true(seize(beds))
  expect_true(beds$is_patient_using(1L))
})

test_that("seize returns FALSE when capacity full", {
  i <- 2L
  curtime <- 1.0
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  expect_false(seize(beds))
  expect_true(beds$is_patient_in_queue(2L))
})

test_that("seize returns NA when capacity full and max_queue reached", {
  i <- 3L
  curtime <- 2.0
  beds <- resource_discrete(1, max_queue = 1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  result <- seize(beds)
  expect_true(is.na(result))
})

test_that("seize with amount=2 blocks 2 units", {
  i <- 1L
  curtime <- 0.0
  beds <- resource_discrete(5)
  seize(beds, amount = 2L)
  expect_equal(beds$n_free(), 3L)
})

# F2.2 release() --------------------------------------------------------------

test_that("release frees resource and returns TRUE when patient was using", {
  i <- 1L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("sick", "sicker", "death"))
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  result <- release(beds)
  expect_true(result)
  expect_equal(beds$n_free(), 1L)
})

test_that("release removes from queue and returns FALSE when patient was queued", {
  i <- 2L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("sick", "sicker", "death"))
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  result <- release(beds)
  expect_false(result)
  expect_equal(beds$queue_size(), 0L)
})

test_that("release when patient not in resource returns FALSE with no error", {
  i <- 99L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("sick", "sicker", "death"))
  beds <- resource_discrete(3)
  result <- release(beds)
  expect_false(result)
})

test_that("release with resume_event triggers next queued patient", {
  i <- 1L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("sicker", "death"))
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  release(beds, resume_event = "sicker")
  expect_equal(queue_size(cur_evtlist), 1L)
  evt <- next_event(ptr = cur_evtlist)
  expect_equal(evt$patient_id, 2L)
  expect_equal(evt$event_name, "sicker")
  expect_equal(evt$time, 5.0)
})

test_that("release with resume_event does not trigger if no one queued", {
  i <- 1L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("sicker", "death"))
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  release(beds, resume_event = "sicker")
  expect_equal(queue_size(cur_evtlist), 0L)
})

test_that("release does not trigger event when patient was in queue (not using)", {
  i <- 2L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("sicker", "death"))
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  beds$attempt_block(patient_id = 3L, priority = 1L, start_time = 2)
  release(beds, resume_event = "sicker")  # patient 2 was in queue, not using
  expect_equal(queue_size(cur_evtlist), 0L)
})

test_that("release with amount=2 frees 2 units", {
  i <- 1L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("sicker", "death"))
  beds <- resource_discrete(5)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0, amount = 2L)
  release(beds, amount = 2L)
  expect_equal(beds$n_free(), 5L)
})

# F2.3 seize_all() ------------------------------------------------------------

test_that("seize_all all_or_none acquires all resources when available", {
  i <- 1L
  curtime <- 0.0
  beds <- resource_discrete(2)
  surgeons <- resource_discrete(2)
  result <- seize_all(list(beds, surgeons))
  expect_true(result)
  expect_true(beds$is_patient_using(1L))
  expect_true(surgeons$is_patient_using(1L))
})

test_that("seize_all all_or_none queues only for first unavailable resource", {
  i <- 1L
  curtime <- 0.0
  beds <- resource_discrete(1)
  surgeons <- resource_discrete(1)
  beds$attempt_block(patient_id = 99L, priority = 1L, start_time = 0)
  result <- seize_all(list(beds, surgeons))
  expect_false(result)
  expect_true(beds$is_patient_in_queue(1L))
  expect_false(surgeons$is_patient_using(1L))
  expect_false(surgeons$is_patient_in_queue(1L))
})

test_that("seize_all sequential holds partial acquisitions", {
  i <- 1L
  curtime <- 0.0
  beds <- resource_discrete(1)
  surgeons <- resource_discrete(1)
  surgeons$attempt_block(patient_id = 99L, priority = 1L, start_time = 0)
  result <- seize_all(list(beds, surgeons), policy = "sequential")
  expect_false(result)
  expect_true(beds$is_patient_using(1L))   # holds beds
  expect_true(surgeons$is_patient_in_queue(1L))  # queued for surgeons
})

test_that("seize_all with custom amounts works", {
  i <- 1L
  curtime <- 0.0
  beds <- resource_discrete(5)
  nurses <- resource_discrete(10)
  result <- seize_all(list(beds, nurses), amounts = c(2L, 3L))
  expect_true(result)
  expect_equal(beds$n_free(), 3L)
  expect_equal(nurses$n_free(), 7L)
})

test_that("seize_all with custom priorities works", {
  i <- 1L
  curtime <- 0.0
  beds <- resource_discrete(1)
  beds$attempt_block(patient_id = 99L, priority = 1L, start_time = 0)
  result <- seize_all(list(beds), priorities = c(5L))
  expect_false(result)
  expect_equal(beds$next_patient_in_line(1L), 1L)
})

# F2.4 release_all() ----------------------------------------------------------

test_that("release_all frees all resources", {
  i <- 1L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("surgery", "death"))
  beds <- resource_discrete(2)
  surgeons <- resource_discrete(2)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  surgeons$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  release_all(list(beds, surgeons))
  expect_equal(beds$n_free(), 2L)
  expect_equal(surgeons$n_free(), 2L)
})

test_that("release_all with single resume_event triggers queued patients for each resource", {
  i <- 1L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("surgery", "death"))
  beds <- resource_discrete(1)
  surgeons <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  surgeons$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  surgeons$attempt_block(patient_id = 3L, priority = 1L, start_time = 1)
  release_all(list(beds, surgeons), resume_event = "surgery")
  expect_equal(queue_size(cur_evtlist), 2L)
})

test_that("release_all with per-resource resume_events uses correct event per resource", {
  i <- 1L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("surgery", "consult", "death"))
  beds <- resource_discrete(1)
  surgeons <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  surgeons$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  beds$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  surgeons$attempt_block(patient_id = 3L, priority = 1L, start_time = 1)
  release_all(list(beds, surgeons), resume_event = c("surgery", "consult"))
  expect_equal(queue_size(cur_evtlist), 2L)
  e1 <- next_event(ptr = cur_evtlist); pop_event(ptr = cur_evtlist)
  e2 <- next_event(ptr = cur_evtlist); pop_event(ptr = cur_evtlist)
  all_names <- c(e1$event_name, e2$event_name)
  expect_true("surgery" %in% all_names)
  expect_true("consult" %in% all_names)
})

test_that("release_all only triggers events for resources with queued patients", {
  i <- 1L
  curtime <- 5.0
  cur_evtlist <- queue_create(c("surgery", "death"))
  beds <- resource_discrete(2)
  surgeons <- resource_discrete(1)
  beds$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  surgeons$attempt_block(patient_id = 1L, priority = 1L, start_time = 0)
  surgeons$attempt_block(patient_id = 2L, priority = 1L, start_time = 1)
  release_all(list(beds, surgeons), resume_event = "surgery")
  expect_equal(queue_size(cur_evtlist), 1L)  # only surgeons had a queue
})

# F2.x seize_all + queue_wait_time end-to-end --------------------------------

test_that("seize_all all_or_none: queue_wait_time correct after dequeue via release_all", {
  # P1 seizes both at t=0. P2 tries at t=1 (both full) -> queued for beds (first unavailable).
  # P1 releases at t=5 -> triggers P2. P2 acquires both.
  # Expected: specialists$had_to_queue(2)==0, beds$had_to_queue(2)==1,
  #           beds$queue_wait_time(2)==4.0
  i <- 1L; curtime <- 0.0
  cur_evtlist <- queue_create(c("sicker", "death"))
  beds        <- resource_discrete(1)
  specialists <- resource_discrete(1)

  seize_all(list(beds, specialists))  # P1 acquires both

  i <- 2L; curtime <- 1.0
  r2 <- seize_all(list(beds, specialists))  # P2 queued for beds (first full)
  expect_false(r2)
  expect_equal(beds$had_to_queue(2L), 1L)
  expect_equal(specialists$had_to_queue(2L), 0L)
  expect_equal(beds$queue_wait_time(2L), NA_real_)  # still waiting

  i <- 1L; curtime <- 5.0
  release_all(list(beds, specialists), resume_event = "sicker")  # P1 releases

  # P2's sicker re-fires at t=5
  i <- 2L; curtime <- 5.0
  r2b <- seize_all(list(beds, specialists))
  expect_true(r2b)
  expect_equal(beds$had_to_queue(2L), 1L)         # permanent flag
  expect_equal(beds$queue_wait_time(2L), 4.0)     # 5 - 1 = 4
  expect_equal(specialists$queue_wait_time(2L), NA_real_)  # never queued for specialists
})

test_that("seize_all all_or_none: queue_wait_time for specialist bottleneck", {
  # Beds have plenty of capacity (10), specialists are the bottleneck (1).
  # P1 seizes both. P2 tries -> gets beds directly (Phase 1) -> specialists full -> queued
  # for specialists.
  i <- 1L; curtime <- 0.0
  cur_evtlist <- queue_create(c("sicker", "death"))
  beds        <- resource_discrete(10)
  specialists <- resource_discrete(1)

  seize_all(list(beds, specialists))  # P1 acquires both

  i <- 2L; curtime <- 2.0
  r2 <- seize_all(list(beds, specialists))  # P2: beds available, specialists full
  expect_false(r2)
  expect_equal(beds$had_to_queue(2L), 0L)         # beds NOT the bottleneck
  expect_equal(specialists$had_to_queue(2L), 1L)  # specialists IS the bottleneck

  i <- 1L; curtime <- 7.0
  release_all(list(beds, specialists), resume_event = "sicker")

  i <- 2L; curtime <- 7.0
  r2b <- seize_all(list(beds, specialists))
  expect_true(r2b)
  expect_equal(specialists$queue_wait_time(2L), 5.0)  # 7 - 2 = 5
  expect_equal(beds$queue_wait_time(2L), NA_real_)    # P2 never queued for beds
})

# F2.5 shared_incr / shared_decr ----------------------------------------------

test_that("shared_incr increments and returns new value", {
  constrained <- TRUE
  obj <- shared_input(10)
  expect_equal(shared_incr(obj, 3), 13)
  expect_equal(obj$value(), 13)
})

test_that("shared_decr decrements and returns new value", {
  constrained <- TRUE
  obj <- shared_input(10)
  shared_incr(obj, 5)
  expect_equal(shared_decr(obj, 2), 13)
  expect_equal(obj$value(), 13)
})

test_that("shared_incr default delta is 1", {
  constrained <- TRUE
  obj <- shared_input(10)
  expect_equal(shared_incr(obj), 11)
})

test_that("shared_incr works with non-integer values", {
  constrained <- TRUE
  obj <- shared_input(1.5)
  expect_equal(shared_incr(obj, 0.3), 1.8, tolerance = 1e-10)
})

# Regression: seize_all partial acquisition bug ---------------------------------

test_that("seize_all all_or_none does not partially acquire when queue exists on one resource", {
  beds <- resource_discrete(2)
  specialists <- resource_discrete(1)
  # P1 acquires both
  i <- 1L; curtime <- 1.0
  seize_all(list(beds, specialists))
  # P2: beds free (n_free=1), specialists full → queued for specialists
  i <- 2L; curtime <- 2.0
  r2 <- seize_all(list(beds, specialists))
  expect_false(r2)
  expect_true(specialists$is_patient_in_queue(2L))
  expect_false(beds$is_patient_using(2L))
  # Free specialist (P1) — now specialists has 1 free but P2 in queue
  i <- 1L
  specialists$attempt_free()
  # State: beds: 1 free, no queue. specialists: 1 free, P2 in queue.
  # P3 tries seize_all — should NOT partially acquire beds
  i <- 3L; curtime <- 5.0
  r3 <- seize_all(list(beds, specialists))
  expect_false(r3)
  expect_equal(beds$n_free(), 1L)
  expect_false(beds$is_patient_using(3L))
  expect_true(specialists$is_patient_in_queue(3L))
  i <- 2L; curtime <- 6.0
  r2 <- seize_all(list(beds, specialists))
  expect_true(r2)
  expect_false(specialists$is_patient_in_queue(2L))
  expect_true(beds$is_patient_using(2L))
  expect_true(specialists$is_patient_using(2L))
})
