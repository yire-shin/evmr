# ------------------------------------------------------------
# testthat: basic tests for evmr package
# These tests verify that the evmr() function runs correctly
# for all supported distributions using randomly generated data.
# ------------------------------------------------------------

library(testthat)
library(evmr)

# Fix random seed for reproducibility
set.seed(123)

# ------------------------------------------------------------
# Helper function: check structure of generated r-largest sample
# ------------------------------------------------------------
check_rmat <- function(x, n, r) {

  # Generated object should be a list
  expect_true(is.list(x))

  # The list should contain rmat (matrix of r-largest data)
  expect_true("rmat" %in% names(x))

  xr <- as.matrix(x$rmat)

  # Check dimension of the generated sample
  expect_equal(nrow(xr), n)
  expect_equal(ncol(xr), r)

  # All values should be finite
  expect_true(all(is.finite(xr)))
}

# ------------------------------------------------------------
# Helper function: check structure of evmr() output
# ------------------------------------------------------------
check_evmr_output <- function(fit, model_name, r) {

  # evmr() should return a data.frame
  expect_true(is.data.frame(fit))

  # Number of rows should match number of r orders
  expect_equal(nrow(fit), r)

  # Required columns that every model must contain
  expect_true(all(c(
    "model", "r", "nllh", "mu", "sigma",
    "rl20", "rl20.se",
    "rl50", "rl50.se",
    "rl100", "rl100.se",
    "rl200", "rl200.se"
  ) %in% names(fit)))

  # Model name should match requested model
  expect_true(all(fit$model == model_name))

  # r column should be 1:r
  expect_equal(fit$r, seq_len(r))

  # Parameter values should be finite
  expect_true(all(is.finite(fit$nllh)))
  expect_true(all(is.finite(fit$mu)))
  expect_true(all(is.finite(fit$sigma)))

  # Scale parameter must be positive
  expect_true(all(fit$sigma > 0))

  # Return level estimates should be finite
  expect_true(all(is.finite(fit$rl20)))
  expect_true(all(is.finite(fit$rl20.se)))
  expect_true(all(is.finite(fit$rl50)))
  expect_true(all(is.finite(fit$rl50.se)))
  expect_true(all(is.finite(fit$rl100)))
  expect_true(all(is.finite(fit$rl100.se)))
  expect_true(all(is.finite(fit$rl200)))
  expect_true(all(is.finite(fit$rl200.se)))
}

# ------------------------------------------------------------
# Test 1: rk4d distribution
# ------------------------------------------------------------
test_that("evmr works for rk4d random sample", {

  # Generate r-largest sample from rk4d distribution
  x <- rk4dr(n = 50, r = 2, loc = 10, scale = 3,
             shape1 = 0.1, shape2 = 0.1)

  # Check generated sample structure
  check_rmat(x, n = 50, r = 2)

  # Run evmr estimation
  fit <- evmr(x$rmat, models = "rk4d", num_inits = 10)

  # Validate output
  check_evmr_output(fit, model_name = "rk4d", r = 2)
})

# ------------------------------------------------------------
# Test 2: rglo distribution
# ------------------------------------------------------------
test_that("evmr works for rglo random sample", {

  # Generate r-largest sample from rglo distribution
  x <- rglor(n = 50, r = 2, loc = 10, scale = 3, shape = 0.1)

  check_rmat(x, n = 50, r = 2)

  # Run evmr estimation
  fit <- evmr(x$rmat, models = "rglo", num_inits = 10)

  check_evmr_output(fit, model_name = "rglo", r = 2)
})

# ------------------------------------------------------------
# Test 3: rggd distribution
# ------------------------------------------------------------
test_that("evmr works for rggd random sample", {

  # Generate r-largest sample from rggd distribution
  x <- rggdr(n = 50, r = 2, loc = 10, scale = 3, shape = 0.1)

  check_rmat(x, n = 50, r = 2)

  # Run evmr estimation
  fit <- evmr(x$rmat, models = "rggd", num_inits = 10)

  check_evmr_output(fit, model_name = "rggd", r = 2)
})

# ------------------------------------------------------------
# Test 4: rgd distribution
# ------------------------------------------------------------
test_that("evmr works for rgd random sample", {

  # Generate r-largest sample from rgd distribution
  x <- rgdr(n = 50, r = 2, loc = 10, scale = 2)

  check_rmat(x, n = 50, r = 2)

  # Run evmr estimation
  fit <- evmr(x$rmat, models = "rgd", num_inits = 10)

  check_evmr_output(fit, model_name = "rgd", r = 2)
})

# ------------------------------------------------------------
# Test 5: rld distribution
# ------------------------------------------------------------
test_that("evmr works for rld random sample", {

  # Generate r-largest sample from rld distribution
  x <- rldr(n = 50, r = 2, loc = 10, scale = 2)

  check_rmat(x, n = 50, r = 2)

  # Run evmr estimation
  fit <- evmr(x$rmat, models = "rld", num_inits = 10)

  check_evmr_output(fit, model_name = "rld", r = 2)
})
