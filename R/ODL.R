#' ODL
#'
#' Calculate the logarithm of the Observed Data Likelihood (ODL) when true responses are unknown and repeated tests are performed on a set of test subjects (e.g. genes).
#'
#' @param n_pos Vector of integers specifying the number of tests that produced a positive result for each subject.
#' @param N Integer value giving the total number of repeated tests per subject.
#' @param theta_0 A real number (in range [0,1]) giving the prevalence of negative results for a set of test subjects.
#' @param theta_1 A real number (in range [0,1]) giving the prevalence of positive results for a set of test subjects.
#' @param p00 A real number (in range [0,1]) giving the probability that a test is negative given that the true result is negative (true-negative rate) for a set of test subjects.
#' @param p11 A real number (in range [0,1]) giving the probability that a test is positive given that the true result is positive (true-positive rate) for a set of test subjects.
#' @param p01 A real number (in range [0,1]) giving the probability that a test is negative given that the true result is positive (false-negative rate) for a set of test subjects.
#' @param p10 A real number (in range [0,1]) giving the probability that a test is positive given that the true result is negative (false-positive rate) for a set of test subjects.
#'
#' @return A real value giving the `log(ODL)`.
#' @references Jakobsdottir and Weeks (2007). Estimating prevalence, false-positive rate, and false-negative rate with use of repeated testing when true responses are unknown. Am J Hum Genet 81:1111-1113.
#' @examples
#' @export
ODL <- function(n_pos, N, theta_0, theta_1, p00, p11, p01, p10){
  if(any(N < n_pos)) stop("N must be greater than or equal to n_pos")
  tryCatch({
    n_neg <- N-n_pos
    neg <- theta_0 * (p00^n_neg) * (p10^n_pos)
    pos <- theta_1 * (p01^n_neg) * (p11^n_pos)
    odl <- sum(log(neg + pos))
  },
  error = function(e) stop(paste("'ODL' error:",e))
  )
  return(odl)
}
