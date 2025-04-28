#' Run Euler's method on an ODE.
#'
#' @param dydt The definition of the ODE.  This may be just the right-hand side
#' of the ODE,
#' or you may specify an equation using _t to denote derivative. See examples.
#' @param tlim a list containing the start and end of the timeframe to integrate
#' over
#' @param ic the initial condition.  If this is a scalar ODE then you can do
#' ic=4. If this is a system then you
#' must provide a list of values in the correct order as determined by the
#' listing of varaibles in dydt. You may also provide a list with names, in that
#' case order doesn't matter.
#' @param stepSize length of timestep to take.
#' @examples
#' #Consider the ODE
#' #dy/dt=y*(y-1)
#'
#' #There are a lot of ways to specify this ode.
#'
#' #You can specify just the right-hand side of the ODE.
#' soln=Euler(y*(y-1) ~ t&y, tlim=c(0,10),ic=0.99)
#' soln
#'
#' #If the equation is autonomous, you don't need to specify t.
#' soln=Euler(y*(y-1) ~ y, tlim=c(0,10),ic=0.99)
#' soln
#'
#' #You can write the left-hand side as y_t.  Note that z_t or something else
#' #that doesn't match the state variable won't work.
#' results=Euler(y_t=y*(y-1)~t&y,tlim=c(0,10),ic=0.99)
#' results
#'
#' #This will fail
#' #results=Euler(z_t=y*(y-1)~t&y,tlim=c(0,10),ic=0.99)
#'
#' #It is fine to make a function out of the right hand side.
#' ode=mosaic::makeFun(y*(y-1)~y)
#' results=Euler(ode(y)~y, tlim=c(0,10),ic=0.99)
#' results
#'
#' #Or to name it again,
#' results=Euler(dydt=ode(y)~y, tlim=c(0,10),ic=0.99)
#' results
#'
#' #Adjust the stepsize with stepSize argument.
#' results=Euler(dydt=ode(y)~y, tlim=c(0,10),ic=0.99,stepSize = 0.1)
#' head(results)
#' mosaic::plotPoints(y~t,data=results)
#'
#' #Consider the system of differential equations
#' #dx/dt=-y
#' #dy/dt=x
#'
#' #Best is to combine the equations using c(), and name each equation using
#' #the naming convention _t.
#' #Name the initial conditions too.
#' results=Euler(c(x_t=-y,y_t=x)~x&y,tlim=c(0,10),ic=c(x=1,y=0))
#' results
#'
#' #If we name the equations, then the order of the equations and initial
#' #conditions does not matter.
#' results=Euler(c(y_t=x,x_t=-y)~x&y,tlim=c(0,10),ic=c(x=1,y=0))
#' results
#'
#' results=Euler(c(y_t=x,x_t=-y)~x&y,tlim=c(0,10),ic=c(y=0,x=1))
#' results
#'
#' #If either the equations or the initial condition is unlabels, then we infer
#' #the order from the order of the state variables.
#' #This corresponds to
#' #dx/dt=-y, x(0)=0
#' #dy/dt=x,  y(0)=1
#' results=Euler(c(-y,x)~x&y,tlim=c(0,10),ic=c(0,1))
#' results
#'
#' #But this corresponds to
#' #dy/dt=-y, y(0)=0
#' #dxdt=x,   x(0)=1
#' results=Euler(c(-y,x)~y&x,tlim=c(0,10),ic=c(0,1))
#' results
#' @returns A data frame containing the approximations of the solutions and the
#' corresponding times.
#' @export
Euler=function(dydt,tlim,ic,stepSize=(tlim[[2]]-tlim[[1]])/10,...){
  t0=tlim[[1]];
  t_final=tlim[[2]]
  #Ok, first, if ... is not empty, and the name of the variable ends in _t, then the user is specifying the DE in the format
  #y_t=t*y~y or something. We should just set dydt to an expression like c(y_t=t*y)~y
  dots=list(...)
  if (length(dots)>0){
    #check to see if this is specifying a derivative.
    derivstr=names(dots)[1]
    if (substr(derivstr,nchar(derivstr)-1,nchar(derivstr))=="_t"){
      #here is the formula for the derivative.
      dform=deparse(mosaic::lhs(dots[[derivstr]]))
      #here are the variables
      dformvars=deparse(mosaic::rhs(dots[[derivstr]]))
      #Put this in the form c(y_t=t*y)~y
      dydt=stats::as.formula(paste("c(",derivstr,"=",dform,")~",dformvars))
    }
  }

  #It is really important that we get the order of the state variables in the correct order. If the statefunction does not
  #return a vector with names, we are going to assume that it is returning derivatives in the order specified by the inputs.
  #Therefore, keeping this order is important. So I'll parse myself

  dydtstr=as.character(dydt);

  #Get the variable declarations.
  RHS=dydtstr[[3]]

  #split at the & and delete whitespace.
  allInputVars=strsplit(gsub("[[:space:]]", "", RHS),split="&")[[1]]

  #Make a list of the state variables.  The order of the state variables comes
  #from the order in which they were listed in input. Check to see if t was listed
  #as an input variable.
  stateVarsInputOrder=c();
  isAutonomous=TRUE
  for (i in 1:length(allInputVars)){
    if (allInputVars[[i]]!="t"){
      stateVarsInputOrder=c(stateVarsInputOrder,allInputVars[[i]])
    }
    else{
      isAutonomous=FALSE;
    }
  }

  #If the function was declared without t as a variable, add t to the list of variables
  #so that calling is simpler later.
  if (isAutonomous){
    dydt=stats::as.formula(paste(dydtstr[[2]],"~t&",dydtstr[[3]],collapse=" "))
  }

  #now get the initial condition right.


  #Were the correct number of initial conditions given?
  if (length(ic)!=length(stateVarsInputOrder)){
    cat(paste("Error. You need to supply ",length(stateVarsInputOrder)," initial conditions. You supplied ",length(ic),"."))
    return()
  }

  #If there were no names given in the IC, then we assume that the IC is in the "correct" order, i.e., in the order the state variables
  #were given as input.
  if (is.null(names(ic))){
    names(ic)=stateVarsInputOrder
  }

  #Now, check to see that the names given in the initial condition match the state variables.
  if (length(setdiff(names(ic),stateVarsInputOrder))>0 |
      length(setdiff(stateVarsInputOrder,names(ic)))>0){
    cat(paste("Error. You must give initial conditions with names "),paste(stateVarsInputOrder,collapse=",")," you gave initial conditions with names ",paste(names(ic),collapse=","))
    return()
  }

  #Now make sure the initial condition is given in an order that matches the state variable order.
  y0=ic[stateVarsInputOrder]

  #Make a function out of dydt.
  dydtfunc=mosaicCore::makeFun(dydt);

  #All the derivatives should be names as stateVar_t. Make a list of what all the derivatives should be.
  derivnames=paste0(stateVarsInputOrder,"_t")

  #Do some testing on dydt to make sure-ish that it is well-defined. To do this evaluate it at the initial condition.
  dydtvalue=do.call(dydtfunc,args=as.list(c(t=t0,y0)))

  #First, is it returning the correct number of outputs?
  if (length(dydtvalue)!=length(stateVarsInputOrder)){
    cat(paste("Error.  You supplied ",length(stateVarsInputOrder)," state variables, but gave ", length(dydtvalue), "differential equations."))
    return();
  }

  #dydtfunc might be returning vectors that are unnamed. To test, evaluate at the initial condition.
  #If it is, create a wrapper that gives the return values the correct names.
  basedydtfunc=dydtfunc
  if (length(names(dydtvalue))==0){
    dydtfunc=function(...){
      output=basedydtfunc(...)
      names(output)=derivnames
      return(output)
    }
  }

  dydtvalue=do.call(dydtfunc,args=as.list(c(t=t0,y0)))

  #Ok, now dydtfunc is for sure returning the right number of things, and they are named. Are they named correctly?
  if (length(setdiff(names(dydtvalue),derivnames))!=0 |
      length(setdiff(derivnames,names(dydtvalue)))!=0){
    cat(paste("Error, your ODEs should have ",paste(derivnames,collapse=","),"on the left-hand side. Instead they have ",paste(names(dydtvalue),collapse=",")))
    return()
  }

  #Ok, now dydtvalue has the right number of entries, named correctly, now are they coming out in the same order as the input goes in?
  #If not, create a wrapper that reorders the outputs into the input order.
  basedydtfunc2=dydtfunc
  if (!all(names(dydtvalue)==derivnames)){
    dydtfunc=function(...){
      output=basedydtfunc2(...)
      output=output[derivnames]
      return(output)
    }
  }

  #At this point, dydtfunc should be ready to go.  It returns outputs in the same order as inputs, has the right number, etc.
  numSteps = ceiling((t_final-t0)/stepSize)
  solution=matrix(ncol=length(stateVarsInputOrder)+1,nrow=numSteps+1)

  solution[1,]=c(t0,y0)
  yi=y0;
  ti=t0;
  for (i in 2:(numSteps+1)){
    yi=yi+stepSize*do.call(dydtfunc,args=as.list(c(t=ti,yi)))
    ti=ti+stepSize;
    solution[i,]=c(ti,yi);
  }
  solution = as.data.frame(solution)
  names(solution)=c("t",stateVarsInputOrder)
  return(solution)
}



