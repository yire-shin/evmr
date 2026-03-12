test_that("rk4d.summary returns a data frame", {
  data(bangkok, package = "evmr")
  out <- rk4d.summary(bangkok, num_inits = 100)
  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 1)
  expect_true(all(c("r", "nllh", "mu", "sigma", "xi", "h") %in% names(out)))
})
