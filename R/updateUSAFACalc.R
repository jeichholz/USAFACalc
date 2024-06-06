#' Updates the USAFACalc package, if update is available, from github. Does not
#' seek to update dependencies.
#' @param repo the repo to install from. Defaults to proper value.
#' @inheritParams devtools::install_github
#' @param ... will be passed into devtools::install_github
#' @export
updateUSAFACalc = function(repo="jeichholz/USAFACalc",upgrade="never",...){
  devtools::install_github(repo,upgrade=upgrade,...)
}
