hello=function(){
  print("Hello there.")
}

#' Run Euler's method on an ODE.
#'
#' @param dydt an expression of t and y, in that order,  defining the right hand side of the ODE.
#' @param tlim a list containing the start and end of the timeframe to integrate over
#' @param y0 the initial condition
#' @param stepSize length of timestep to take.
#' @examples
#' soln=Euler(y*(y-1) ~ t&y, tlim=c(0,10),y0=0.99)
#' soln
#'
#' results=Euler(5*y*(y-1)~t&y,tlim=c(0,10),y0=0.99,stepSize=0.1)
#' results
#' mosaic::plotPoints(y~t,data=results)
#' @returns A data frame containing the approximations of the solutions and the corresponding times.
#' @export

Euler=function(dydt,tlim,y0,stepSize=(tlim[[2]]-tlim[[1]])/10){
  t0=tlim[[1]];
  t_final=tlim[[2]]
  dydtfunc=mosaicCore::makeFun(dydt);
  numSteps = ceiling((t_final-t0)/stepSize)
  return(mosaicCalc::Iterate(function(t,y) c(t+stepSize,y+stepSize*dydtfunc(t,y)),x0=c(t0,y0),n=numSteps))

}

