#' A replacement for mosaic::antiD that uses sympy.
#' WARNING: sympy is very powerful and might return the integral
#' in terms of functions that R doesn't know about, like it might
#' return things in terms of Ci, the cosine integral.  In that case,
#' antiD will succeed, but you won't be able to evaluate the resulting function.
#' @param expression the expression to integrate.
#' @param ... additional arguments defining parameters that might exist in expression.
#' @examples
#' antiD(sin(x)*x~x)
#' antiD(exp(-z^2)~z)

#' @export
antiD=function(expression,...){


    tryCatch({
    #Grab the left- and right- sides of the tilde expression.
    lh=mosaicCore::lhs(expression)
    rh=mosaicCore::rhs(expression)
    #Get the name of the variable. Parameters not supported.
    vars=all.vars(rh)

    if (length(vars)>1){
      stop("antiDv1 only supports single integrals.\n")
    }

    #Create a function out of the expression, for the purposes of substituting in the sympy symbol.
    f=mosaic::makeFun(expression,...)

    x=caracas::symbol(vars[[1]])

    #Integrate, so easy!
    F=caracas::int(f(x),x);

    #Turn into a string.
    Fstr=caracas::as_character(F);

    failed=FALSE
    #If sympy returned an unevaluted integral, tell the user to give up and use sympy.
    if (grepl("Integral",Fstr)){
      warning("Can't find an antiderivative.  Use integrate() instead.")
      failed=TRUE
    }

    #Start creating the R function.  Make a string to evaluate to define it.
    GDeclaration=paste("G=function(",vars[[1]],",C=0)")
    GDeclaration=paste(GDeclaration,"{")

    #One common function that sympy might return is erf.  That isn't defined in R, but
    #pnorm is.  If this happens, define erf in terms of pnorm locally inside your anti derivative.
    #One day maybe work on reliably simplifying this.
    if (grepl("erf",Fstr)){
      GDeclaration=paste(GDeclaration,"erf <- function(x){ 2*pnorm(x*sqrt(2))-1};")
    }
    if (grepl("erfi",Fstr)){
      GDeclaration=paste(GDeclaration,"erfi <- function(x) {qnorm((1 + x)/2)/sqrt(2)};")
    }

    #Add on the constant of integration.
    GDeclaration=paste(GDeclaration,Fstr,"+C")
    GDeclaration=paste(GDeclaration,"}")

    #Create the function, and return.
    eval(str2expression(GDeclaration))
    if (!failed){
      return(G)
    }
    },
    error = function(msg){
      warning("USAFACalc antiD ran into an error:")
      warning(paste(msg))
      warning("Falling back to mosaicCalc antiD")
      return(antiDFromMosaic(expression,...))
    })


}


#' Calculate an antiderivative, error if no symbolic antiderivative  A wrapper around
#' mosaicCalc::antiD.
#' @param tilde  the expression to integrate
#' @returns A function that is the symbolic antiderivative of tilde. Error if no
#' symbolic antiderivative found.
#' @examples
#' #Succeeds
#' #antiD(x^2~x)
#'
#' #Succeeds
#' #antiD(cos(2*x)~x)
#'
#' #This will fail:
#' #antiD(dnorm(x)~x)
#'
#' #If you really want to use numerical anti differentiation for an antiderivative
#' #do this instead:
#' #mosaic::makeFun(integrate(dnorm(s)~s,slim=c(0,x))~x)

antiDFromMosaic=function(tilde){

  #Call antiD.  This won't fail!
  F=mosaicCalc::antiD(tilde);

  #If it failed by resorting to a numerical anti derivative, then throw an error.
  if (grepl("evalFun",paste(deparse(F),collapse=""))){
    stop("Can't find an anti-derivative. Try numerical antiderivative using integrate().")
  }

  return(F)
}

