#' estPerf
#'
#' Estimate performance (false-positive rate, false-negative rate, true-positive rate, true-negative rate) when there are repeated tests for a set of test subjects (e.g. genes).
#'
#' @param n Vector of integers specifying the number of tests that produced a particular result (either negative or positive) for each subject.
#' @param N Integer value giving the total number of repeated tests per subject.
#' @param delta Vector of real numbers giving the probability that a given gene has the expected result (either negative or positive) - must be the same length as `n`.
#'
#' @return The estimated performance, a real value in the range [0,1].
#' @references Jakobsdottir and Weeks (2007). Estimating prevalence, false-positive rate, and false-negative rate with use of repeated testing when true responses are unknown. Am J Hum Genet 81:1111-1113.
#' @examples \dontrun{estPerf(n,N,delta)}
#' @export
estPerf <- function(n, N, delta){
  if(length(n) != length(delta)) stop("n and delta must be the same length")
  tryCatch({
    n_opp <- N-n
    num <- sum(n*delta)
    denom <- num + sum(n_opp*delta)
    res <- num/denom
  },
  error = function(e) stop(paste("'estPerf' error:",e))
  )
  # Sanity check result.
  if(res < 0) stop(paste("Negative rate returned:",res))
  if(res > 1) stop(paste("Rate exceeds 1:",res))
  return(res)
}
