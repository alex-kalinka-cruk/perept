#' estPrev
#'
#' Estimate the prevalence of a particular result for a set of test subjects when there are repeated tests per subject (e.g. gene).
#'
#' @param delta Vector of real numbers giving the probability that a given gene has the expected result (either negative or positive).
#'
#' @result A real number giving the estimated prevalence of a given result (either negative or positive) for the set of subjects.
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
