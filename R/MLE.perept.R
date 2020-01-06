#' MLE.perept
#'
#' Maximum Likelihood Estimation (MLE) of the performance of a test when there are repeated tests performed on the same subjects.
#'
#' @param pos_counts A data frame of positive counts across a set of repeated tests with `n` columns. The first column is expected to contain subject IDs, columns 2 to `n-1` should contain positive counts across subjects for the `n-2` tests, and the `n`'th column should contain the total positive count across tests for each subject (the row sums).
#' @param delta_pos A vector of real numbers (in the range [0,1]) giving initial estimates of the probability of each subject being positive - must be the same length as `n_pos`. If `NULL`, then `n_pos`/`N` is used. Defaults to `NULL`.
#'
#' @param maxiter A positive integer specifying the maximum allowable number of iterations without convergence occurring. Defaults to 1e5.
#' @param tol A real number specifying the maximum ODL difference between successive steps of the EM algorithm below which convergence occurs. Defaults to 1e-6.
#'
#' @return An object of class `MLE.perept`, which is a list with the following elements:
#' @references Jakobsdottir and Weeks (2007). Estimating prevalence, false-positive rate, and false-negative rate with use of repeated testing when true responses are unknown. Am J Hum Genet 81:1111-1113.
#' @examples
#' @export
MLE.perept <- function(pos_counts, delta_pos = NULL, maxiter = 1e5, tol = 1e-6){
  if(class(pos_counts) != "data.frame")
    stop("'pos_counts' must be a data frame of positive result counts from repeated tests")
  n <- ncol(pos_counts)
  N <- n-2
  if(N < 3) stop(paste("a minimum of 3 repeated tests per subject is required: found",N))
  pos_tests <- unlist(pos_counts[n])

  # Perform Expectation Maximisation algorithm for performance MLE.
  em <- perept::EM.perept(pos_tests, N, delta_pos, maxiter, tol)
  return(em)
}
