context("test-perept.R")

# Data on p.95 of Wong et al. (2007) A Comprehensive Analysis of Common Copy-Number Variations in the Human Genome [Wetal2007].
# This data was re-analysed by Jakobsdottir and Weeks (2007) [JW2007] with their updated results on p.1111:
# theta_0 = 0.006299; p10 = 0.002283; p01 = 0.4891.
# Wetal2007's raw positive result counts per clone (across 6 tests):
wong_etal.n_neg <- c(rep(0,23911),rep(1,340),rep(2,50),rep(3,46),rep(4,15),rep(5,15),rep(6,15))


### Tests
# Compare against JW2007 re-analysis of Wetal2007's data.
testthat::test_that("we get the same result as JW2007 using Wetal2007", {
  # Run EM algorithm on Wetal2007 data.
  res <- EM.perept(wong_etal.n_neg, 6, wong_etal.n_neg/6)
  # Tolerance set to number of decimal places reported in JW2007.
  expect_equal(0.006299, res$prevalence$theta_1, tol=1e-6)
  expect_equal(0.002283, res$performance$p10, tol=1e-6)
  expect_equal(0.4891, res$performance$p01, tol=1e-4)
})
