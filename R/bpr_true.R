#' bpr_true
#'
#' Calculates a Bayesian posterior probability that a given test subject (e.g. a gene) has a given expected result (either negative or positive).
#'
#' @param n An integer specifying the number of tests that produced a particular result (either negative or positive) for the focal subject.
#' @param N Integer value giving the total number of repeated tests for the given subject.
#' @param expected An integer value (either 0 or 1) indicating whether the expected result is negative (0) or positive (1).
#' @param theta_0 A real number (in range [0,1]) giving the prevalence of negative results for a set of test subjects.
#' @param theta_1 A real number (in range [0,1]) giving the prevalence of positive results for a set of test subjects.
#' @param p00 A real number (in range [0,1]) giving the probability that a test is negative given that the true result is negative (true-negative rate) for a set of test subjects.
#' @param p11 A real number (in range [0,1]) giving the probability that a test is positive given that the true result is positive (true-positive rate) for a set of test subjects.
#' @param p01 A real number (in range [0,1]) giving the probability that a test is negative given that the true result is positive (false-negative rate) for a set of test subjects.
#' @param p10 A real number (in range [0,1]) giving the probability that a test is positive given that the true result is negative (false-positive rate) for a set of test subjects.
#'
#' @result A real number giving the Bayesian posterior probability that the test subject has a given result.
#' @examples
#' @export
bpr_true <- function(n, N, expected, theta_0, theta_1, p00, p11, p01, p10){
  if(!expected %in% 0:1) stop("expected result can only be 0 (negative) or 1 (positive)")
  if(N < n) stop("N must be greater than or equal to n")
  n_opp <- N-n
  tryCatch({
    if(expected == 0){
      num <- theta_0 * (p00^n) * (p10^n_opp)
      denom <- num + (theta_1 * (p01^n) * (p11^n_opp))
    }else{
      num <- theta_1 * (p01^n_opp) * (p11^n)
      denom <- num + (theta_0 * (p00^n_opp) * (p10^n))
    }
    res <- num/denom
  },
  error = function(e) stop(patse("'bpr_true' error:",e))
  )
  # Sanity check result.
  if(res < 0) stop(paste("Negative posterior probability returned:",res))
  if(res > 1) stop(paste("Posterior probability exceeds 1:",res))
  return(res)
}
