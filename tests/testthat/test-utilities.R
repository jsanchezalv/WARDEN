# Tests for R/utilities.R

test_that("waning_hr returns base_hr before waning starts", {
  expect_equal(waning_hr(0.7, 0, 2, 3), 0.7)
  expect_equal(waning_hr(0.7, 1.99, 2, 3), 0.7)
})

test_that("waning_hr returns 1 after waning ends", {
  expect_equal(waning_hr(0.7, 5, 2, 3), 1.0)
  expect_equal(waning_hr(0.7, 100, 2, 3), 1.0)
})

test_that("waning_hr linear interpolates correctly", {
  expect_equal(waning_hr(0.7, 3.5, 2, 3), 0.7 + 0.5 * 0.3)
  expect_equal(waning_hr(0.7, 2, 2, 3), 0.7)
})

test_that("waning_hr exponential waning works", {
  result <- waning_hr(0.5, 3.0, 2, 2, type = "exponential")
  frac <- (3.0 - 2) / 2
  expected <- 0.5 * exp(frac * log(1 / 0.5))
  expect_equal(result, expected)
})

test_that("waning_hr instant waning jumps to 1", {
  expect_equal(waning_hr(0.7, 2.01, 2, 3, type = "instant"), 1.0)
  expect_equal(waning_hr(0.7, 1.99, 2, 3, type = "instant"), 0.7)
})

test_that("waning_hr errors on unknown type", {
  expect_error(waning_hr(0.7, 3, 2, 3, type = "cubic"), "unknown type")
})

test_that("draw_categorical picks correct bin", {
  expect_equal(draw_categorical(0, c(0.3, 0.5, 0.2)), 1L)
  expect_equal(draw_categorical(0.29, c(0.3, 0.5, 0.2)), 1L)
  expect_equal(draw_categorical(0.31, c(0.3, 0.5, 0.2)), 2L)
  expect_equal(draw_categorical(0.79, c(0.3, 0.5, 0.2)), 2L)
  expect_equal(draw_categorical(0.81, c(0.3, 0.5, 0.2)), 3L)
  expect_equal(draw_categorical(0.99, c(0.3, 0.5, 0.2)), 3L)
})

test_that("draw_categorical normalizes probs", {
  expect_equal(draw_categorical(0.5, c(1, 1, 1)), 2L)
  expect_equal(draw_categorical(0.01, c(2, 2, 2)), 1L)
})

test_that("draw_categorical errors on negative probs", {
  expect_error(draw_categorical(0.5, c(0.3, -0.1, 0.8)), "must be >= 0")
})

test_that("qcategorical is an alias for draw_categorical", {
  expect_equal(qcategorical(0.5, c(0.3, 0.5, 0.2)), draw_categorical(0.5, c(0.3, 0.5, 0.2)))
})
