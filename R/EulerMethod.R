hello=function(){
  print("Hello there.")
}

#' Run Euler's method on an ODE.
#'
#' @param dydt an expression defining the right hand side of the ODE. If this is a system of ODE,
#' then the order of the equations is the order in which the variables are listed in the expression. t must be the time variable.
#' @param tlim a list containing the start and end of the timeframe to integrate over
#' @param ic the initial condition.  If this is a scalar ODE then you can do ic=4. If this is a system then you
#' must provide a list of values in the correct order as determined by the listing of varaibles in dydt. You may also
#' provide a list with names, in that case order doesn't matter.
#' @param stepSize length of timestep to take.
#' @examples
#' soln=Euler(y*(y-1) ~ t&y, tlim=c(0,10),ic=0.99)
#' soln
#'
#' results=Euler(5*y*(y-1)~t&y,tlim=c(0,10),ic=0.99,stepSize=0.1)
#' results
#' mosaic::plotPoints(y~t,data=results)
#'
#'
#' results=Euler(c(-y,x)~x&y,tlim=c(0,10),ic=c(1,0),stepSize=0.01)
#' head(results)
#' mosaic::plotPoints(y~x,data=results)
#'
#' results=Euler(c(-y,x)~x&y,tlim=c(0,10),ic=c(y=0,x=1),stepSize=0.01)
#' head(results)
#' mosaic::plotPoints(y~x,data=results)
#' @returns A data frame containing the approximations of the solutions and the corresponding times.
#' @export

Euler=function(dydt,tlim,ic,stepSize=(tlim[[2]]-tlim[[1]])/10){
  t0=tlim[[1]];
  t_final=tlim[[2]]

  #It is really important that we get the order of the state variables in the correct order.
  #Therefore, I'm going to parse this myself rather than relying on built in functions, because I
  #am unsure of the return order.

  dydtstr=as.character(dydt);

  #Get the variable declarations.
  RHS=dydtstr[[3]]

  #split at the & and delete whitespace.
  allvars=strsplit(gsub("[[:space:]]", "", RHS),split="&")[[1]]

  #Check to see if t was a listed variable.  If it is, remove it from the list
  #of state variables. If it is not, add it to the function declaration.


  stateVars=c();
  isAutonomous=TRUE
  for (i in 1:length(allvars)){
    if (allvars[[i]]!="t"){
      stateVars=c(stateVars,allvars[[i]])
    }
    else{
      isAutonomous=FALSE;
    }
  }

  #If the function was declared without t as a variable, add t to the list of variables
  #so that calling is simpler later.
  if (isAutonomous){
    dydt=as.formula(paste(dydtstr[[2]],"~t&",dydtstr[[3]],collapse=" "))
  }

  #now get the initial condition right.  If there were no names given in the IC,
  #then we assume that the IC is in the "correct" order, i.e., in the order the variables
  #were given.  If there are names, then create y0 in the correct order.

  if (is.null(names(ic))){
    y0=ic;
  }
  else{
    y0=list()
    for (i in 1:length(stateVars)){
      y0[[i]]=ic[[stateVars[[i]]]]
    }
    y0=unlist(y0);
  }

  dydtfunc=mosaicCore::makeFun(dydt);
  numSteps = ceiling((t_final-t0)/stepSize)
  solution=data.frame(matrix(ncol=length(stateVars)+1,nrow=numSteps+1))
  names(solution)=c("t",stateVars)


  solution[1,]=c(t0,y0)
  yi=y0;
  ti=t0;
  for (i in 2:(numSteps+1)){
    yi=yi+stepSize*do.call(dydtfunc,args=as.list(setNames(c(ti,yi),c("t",stateVars))))
    ti=ti+stepSize;
    solution[i,]=c(ti,yi);
  }

  #return(mosaicCalc::Iterate(function(t,y) c(t+stepSize,y+stepSize*dydtfunc(t,y)),x0=c(t0,y0),n=numSteps))
  return(solution)
}

