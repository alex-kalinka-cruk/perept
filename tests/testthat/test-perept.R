context("test-perept.R")

# Data on p.95 of Wong et al. (2007) A Comprehensive Analysis of Common Copy-Number Variations in the Human Genome [Wetal2007].
# This data was re-analysed by Jakobsdottir and Weeks (2007) [JE2007] with their updated results on p.1111:
# theta_0 = 0.006299; p10 = 0.002283; p01 = 0.4891.
# Wong et al.'s raw negative result counts per clone:
wong_etal.n_neg <- c(rep(6,23911),rep(5,340),rep(4,50),rep(3,46),rep(2,15),rep(1,15),rep(0,15))


### Tests
# Compare against Jakobsdottir and Weeks re-analysis of Wong et al.'s data.
testthat::test_that("we get the same result as JE2007 using Wetal2007", {
  res <- EM.perept(wong_etal.n_neg, 6, wong_etal.n_neg/6)
  # Tolerance impacted by limited number of decimal places reported in JE2007.
  expect_equal(0.006299, res$prevalence$theta_1, tol=1e-6)
  expect_equal(0.002283, res$performance$p10, tol=1e-6)
  expect_equal(0.4891, res$performance$p01, tol=1e-4)
})
