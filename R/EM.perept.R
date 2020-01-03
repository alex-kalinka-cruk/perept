#' EM.perept
#'
#' Expectation Maximisation (EM) algorithm to estimate false-positive and false-negative rates of a test when the true response is unknown and repeated tests are available for a set of test subjects (e.g. genes).
#'
#' @param n_pos Vector of integers specifying the number of tests that produced a positive result for each subject.
#' @param N Integer value giving the total number of repeated tests per subject.
#' @param delta_pos A vector of real numbers (in the range [0,1]) giving initial estimates of the probability of each subject being positive - must be the same length as `n_pos` (using `n_pos`/`N` works well).
#' @param maxiter A positive integer specifying the maximum allowable number of iterations without convergence occurring. Defaults to 1e5.
#' @param tol A real number specifying the maximum ODL difference between successive steps of the EM algorithm below which convergence occurs. Defaults to 1e-6.
#'
#' @return A list with the following elements:
#' @examples
#' @export
EM.perept <- function(n_pos, N, delta_pos, maxiter = 1e5, tol = 1e-6){
  if(maxiter <= 1) stop("please specify a 'maxiter' value greater than 1")
  if(tol >= 100) stop("please specify a 'tol' value less than 100")
  if(length(n_pos) != length(delta_pos)) stop("'n_pos' and 'delta_pos' must be the same length")
  if(any(N < n_pos)) stop("'N' must be greater than or equal to 'n_pos'")
  if(any(delta_pos < 0) || any(delta_pos > 1)) stop("'delta_pos' must be in the range [0,1]")
  i <- 1
  odl_diff <- 100
  odl_last <- -1e6
  odl_list <- NULL
  ret <- NULL
  tryCatch({
    n_neg <- N-n_pos
    delta_neg <- 1-delta_pos
    while(i < maxiter && odl_diff > tol){
      # 1. Performance estimates.
      p00 <- perept::estPerf(n_neg, N, delta_neg)
      p11 <- perept::estPerf(n_pos, N, delta_pos)
      p01 <- perept::estPerf(n_neg, N, delta_pos)
      p10 <- perept::estPerf(n_pos, N, delta_neg)
      # 2. Prevalence estimates.
      theta_0 <- perept::estPrev(delta_neg)
      theta_1 <- perept::estPrev(delta_pos)
      # 3. Update individual subject classification probabilities.
      delta_neg <- sapply(n_neg, perept::bpr_true,
                          N, 0, theta_0, theta_1, p00, p11, p01, p10)
      delta_pos <- 1-delta_neg
      # 4. ODL estimate.
      odl <- perept::ODL(n_neg, N, theta_0, theta_1, p00, p11, p01, p10)
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
  class(ret) <- "perept"
  return(ret)
}
