test_that("rggd.summary returns a data frame", {
  data(bangkok, package = "evmr")
  out <- rggd.summary(bangkok, num_inits = 100)
  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 1)
  expect_true(all(c("r", "nllh", "mu", "sigma", "xi") %in% names(out)))
})
