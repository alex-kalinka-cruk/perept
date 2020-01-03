#' print.perept
#'
#' print S3 method for perept.
#'
#' @param x An object of class `perept`.
#' @param ... Other arguments to be passed to print.
#'
#' @return Prints a summary to the console screen.
#' @export
print.perept <- function(x, ...){
  if(class(x) != "perept") stop("expecting an object of class 'perept'")
  cat(paste("### perept summary ###\nFalse positive rate: ",round(x$performance$p10,6),
            "\nFalse negative rate: ",round(x$performance$p01,6),
            "\nPrevalence of positives in population (%): ",round(100*x$prevalence$theta_1,6),
            "\nNumber of iterations before convergence: ",length(x$ODL),
            sep=""))
}
