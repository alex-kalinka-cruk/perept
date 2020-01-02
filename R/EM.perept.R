#' EM.perept
#'
#' Expectation Maximisation (EM) algorithm to estimate false-positive and false-negative rates of a test when the true response is unknown and repeated tests are available for a set of test subjects (e.g. genes).
#'
#' @param n_neg Vector of integers specifying the number of tests that produced a negative result for each subject.
#' @param N Integer value giving the total number of repeated tests per subject.
#' @param delta_neg A vector of real numbers (in the range [0,1]) giving initial estimates of the probability of each subject being negative - must be the same length as `n_neg`.
#' @param maxiter A positive integer specifying the maximum allowable number of iterations without convergence occurring.
#' @param tol A real number specifying the maximum ODL difference between successive steps of the EM algorithm below which convergence occurs. Defaults to $1e-6$.
#'
#' @result A list with the following elements:
#' @examples
#' @export
EM.perept <- function(n_neg, N, delta_neg, maxiter = 1e5, tol = 1e-6){
  if(maxiter <= 1) stop("please specify a maxiter value greater than 1")
  if(tol >= 100) stop("please specify a tol value less than 100")
  if(length(n_neg) != length(delta_neg)) stop("n_neg and delta_neg must be the same length")
  if(any(N < n_neg)) stop("N must be greater than n_neg")
  if(any(delta_neg < 0) || any(delta_neg > 1)) stop("delta_neg must be in the range [0,1]")
  i <- 1
  odl_diff <- 100
  odl_last <- -1e6
  odl_list <- NULL
  ret <- NULL
  tryCatch({
    n_pos <- N-n_neg
    delta_pos <- 1-delta_neg
    while(i < maxiter && odl_diff > tol){
      # 1. Performance estimates.
      p00 <- estPerf(n_neg, N, delta_neg)
      p11 <- estPerf(n_pos, N, delta_pos)
      p01 <- estPerf(n_neg, N, delta_pos)
      p10 <- estPerf(n_pos, N, delta_neg)
      # 2. Prevalence estimates.
      theta_0 <- estPrev(delta_neg)
      theta_1 <- estPrev(delta_pos)
      # 3. Update individual subject classification probabilities.
      delta_neg <- sapply(n_neg, bpr_true,
                          N, 0, theta_0, theta_1, p00, p11, p01, p10)
      delta_pos <- 1-delta_neg
      # 4. ODL estimate.
      odl <- ODL(n_neg, N, theta_0, theta_1, p00, p11, p01, p10)
      odl_list <- append(odl_list, odl)
      odl_diff <- odl - odl_last
      odl_last <- odl
      i <- i+1
    }
  },
  error = function(e) stop("EM.perept error:",e)
  )
  ret$performance <- list(p00 = p00, p11 = p11, p01 = p01, p10 = p10)
  ret$prevalence <- list(theta_0 = theta_0, theta_1 = theta_1)
  ret$delta <- list(delta_neg = delta_neg, delta_pos = delta_pos)
  ret$ODL <- odl_list
  return(ret)
}
