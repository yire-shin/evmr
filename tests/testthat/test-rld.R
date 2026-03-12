test_that("rld.summary returns a data frame", {
  data(bangkok, package = "evmr")
  out <- rld.summary(bangkok, num_inits = 5)
  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 1)
  expect_true(all(c("r", "nllh", "mu", "sigma") %in% names(out)))
})
