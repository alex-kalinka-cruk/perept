#' print.EM.perept
#'
#' print S3 method for `EM.perept`` objects.
#'
#' @param x An object of class `EM.perept`.
#' @param ... Other arguments to be passed to print.
#'
#' @return Prints a summary to the console screen.
#' @export
print.EM.perept <- function(x, ...){
  if(class(x) != "EM.perept") stop("expecting an object of class 'EM.perept'")
  cat(paste("### EM.perept summary ###\nFalse positive rate: ",round(x$performance$p10,6),
            "\nFalse negative rate: ",round(x$performance$p01,6),
            "\nPrevalence of positives in population (%): ",round(100*x$prevalence$theta_1,6),
            "\nNumber of iterations before convergence: ",length(x$ODL),
            sep=""))
}
