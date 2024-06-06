#' Updates the USAFACalc package, if update is available, from github. Does not
#' seek to update dependencies.
#' @export
updateUSAFACalc = function(){
  devtools::install_github("jeichholz/USAFACalc",upgrade="never")
}
