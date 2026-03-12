test_that("evmr returns a combined data frame", {
  data(bangkok, package = "evmr")
  out <- evmr(bangkok, models = c("rk4d", "rglo", "rggd", "rld", "rgd"), num_inits = 100)

  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 1)
  expect_true("model" %in% names(out))
  expect_true(all(c("rk4d", "rglo", "rggd", "rld", "rgd") %in% unique(out$model)))
})
