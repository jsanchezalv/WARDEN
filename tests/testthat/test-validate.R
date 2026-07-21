# Tests for R/validate.R

test_that("validate_model passes with valid inputs", {
  evt_react <- add_reactevt(name_evt = "sick", input = {})
  evt_react <- add_reactevt(evt_react, name_evt = "death", input = {})

  init_evts <- add_tte(arm = c("int", "noint"), evts = c("sick", "death"), input = {
    sick <- 1
    death <- 5
  })

  expect_true(
    validate_model(
      arm_list = c("int", "noint"),
      init_event_list = init_evts,
      evt_react_list = evt_react,
      npats = 10,
      n_sim = 1,
      seed = 1
    )
  )
})

test_that("validate_model errors on duplicate arm_list", {
  expect_error(
    validate_model(arm_list = c("int", "int")),
    "duplicates"
  )
})

test_that("validate_model errors on missing arm in init_event_list", {
  init_evts <- add_tte(arm = "int", evts = c("sick", "death"), input = {
    sick <- 1
    death <- 5
  })
  evt_react <- add_reactevt(name_evt = "sick", input = {})
  evt_react <- add_reactevt(evt_react, name_evt = "death", input = {})

  expect_error(
    validate_model(
      arm_list = c("int", "noint"),
      init_event_list = init_evts,
      evt_react_list = evt_react
    ),
    "missing arm"
  )
})

test_that("validate_model errors on missing reaction", {
  init_evts <- add_tte(arm = c("int", "noint"), evts = c("sick", "death"), input = {
    sick <- 1
    death <- 5
  })
  evt_react <- add_reactevt(name_evt = "sick", input = {})

  expect_error(
    validate_model(
      arm_list = c("int", "noint"),
      init_event_list = init_evts,
      evt_react_list = evt_react
    ),
    "missing from"
  )
})

test_that("validate_model errors on seed overflow", {
  expect_error(
    validate_model(seed = 99999999, n_sim = 500, npats = 1000),
    "overflow"
  )
})

test_that("validate_model warns when n_sim > 1 and psa_bool = FALSE", {
  expect_output(
    validate_model(n_sim = 10, psa_bool = FALSE),
    "WARN.*psa_bool"
  )
})
