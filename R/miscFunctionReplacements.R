#' integrate an expression over a finite interval.  A wrapper for stats::integrate
#' @param expression  the expression to integrate
#' @param xlim the interval, as a list, over which to integrate.
#' @returns The definite integral requested.
#' @examples
#' integrate(x^2~x, xlim=c(-4,4))
#' #Integrate returns just a plain number now, so you can do stuff like:
#' int=integrate(x^2~x,xlim=c(-4,4))
#' int
#' @examples
#' integrate(t^2~t,tlim=c(0,5))
#'
#' @export
integrate=function(expression,xlim=NA,...){
  dots=list(...)

  #Find all the variables in the expression.
  rhsVars=all.vars(mosaicCore::rhs(expression))

  lims <- mosaic::inferArgs(dots = dots, vars = rhsVars, defaults = list(xlim = xlim))

  xlim=lims$xlim


  if (any(is.na(xlim))){
    stop("You must provide an interval to integrate over using xlim=...")
  }

  if (any(is.infinite(xlim))){
    stop("You must provide finite limits of integration.")
  }

  func=mosaicCore::makeFun(expression);

  ans=stats::integrate(func,xlim[[1]],xlim[[2]])
  print(ans);
  return(ans$value)
}
