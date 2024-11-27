#' integrate an expression over a finite interval.  A wrapper for stats::integrate
#' @param expression  the expression to integrate
#' @param xlim the interval, as a list, over which to integrate.
#' @param ... instead of xlim, you may provide a limit that reflects the variable name.  For instance, if the variable is t, then tlim=c(-1,1)
#' You may also provide additional arguments that are passed to stats::integrate.
#' @param quietly if TRUE, then don't print the message about possible error.
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
integrate=function(expression,xlim=NA,quietly=FALSE,...){
  dots=list(...)

  #Find all the variables in the expression.
  rhsVars=all.vars(mosaicCore::rhs(expression))

  #Get the limits of integration.  inferArgs lets us accept limits that match the variable name, like tlim or whatever.
  lims <- mosaic::inferArgs(dots = dots, vars = rhsVars, defaults = list(xlim = xlim))
  xlim=lims$xlim

  if (any(is.na(xlim))){
    stop("You must provide an interval to integrate over.")
  }



  #integrate doesn't accept any arguments that end in lim.  So delete anything that ends in lim from dots, assuming
  #that it was a limit definition.
  if (length(dots)!=0){
    dots[endsWith(names(dots),"lim")]<-NULL
  }

  #Default the max number of subintervals way higher than the stats::integrate default.
  if (!("subdivisions" %in% names(dots))){
    dots$subdivisions=1e5
  }

  func=mosaicCore::makeFun(expression);
  ans=do.call(stats::integrate,c(f=func,lower=xlim[[1]],upper=xlim[[2]],dots))
  if (!quietly){
    print(ans);
  }

  return(ans$value)
}



