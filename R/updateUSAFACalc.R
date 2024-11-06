#' Updates the USAFACalc package, if update is available, from github. Does not
#' seek to update dependencies.
#' @param repo the repo to install from. Defaults to proper value.
#' @param upgrade -- upgrade option to install_github.  When set to never it will not ask about upgrading packages upon which USAFACalc depends.
#' @param ... will be passed into devtools::install_github
#' @export
updateUSAFACalc = function(repo="jeichholz/USAFACalc",upgrade="never",...){
  devtools::install_github(repo,upgrade=upgrade,...)
}
