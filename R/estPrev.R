#' estPrev
#'
#' Estimate the prevalence of a particular result for a set of test subjects when there are repeated tests per subject (e.g. gene).
#'
#' @param delta Vector of real numbers giving the probability that a given gene has the expected result (either negative or positive).
#'
#' @return A real number giving the estimated prevalence of a given result (either negative or positive) for the set of subjects.
#' @references Jakobsdottir and Weeks (2007). Estimating prevalence, false-positive rate, and false-negative rate with use of repeated testing when true responses are unknown. Am J Hum Genet 81:1111-1113.
#' @examples \dontrun{estPrev(delta)}
#' @export
estPrev <- function(delta){
  tryCatch(estPrev <- sum(delta)/length(delta),
           error = function(e) stop(paste("'estPrev' error:",e)))
  # Sanity check result.
  if(estPrev < 0) stop(paste("Negative prevalence returned:",estPrev))
  if(estPrev > 1) stop(paste("Prevalence exceeds 1:",estPrev))
  return(estPrev)
}
