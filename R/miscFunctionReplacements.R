#' integrate an expression over a finite interval.  A wrapper for stats::integrate
#' @param expression  the expression to integrate
#' @param xlim the interval, as a list, over which to integrate.
#' @returns The definite integral requested.
#' @examples
#' integrate(x^2~x, xlim=c(-4,4))
#' #Integrate returns just a plain number now, so you can do stuff like:
#' int=integrate(x^2~x,xlim=c(-4,4))
#' int
#' @export
integrate=function(expression,xlim=NA){

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
