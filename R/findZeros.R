#
#' findZeros
#' @description
#' Find a zero of an expression, or system of expressions.
#'
#' findZeros is a drop-in replacement for mosaic::findZeros.
#'
#' The goals are as follows:
#'
#' 1) Correctly find roots of all polynomials, even ones that have a repeated root.
#' 2) Return accurate results at all times, that is, if xstar is returned as a root, it should indeed have the property that f(xstar) is near 0.
#' 3) Be able to handle nested functions correctly.  For instance, f(x)=x+1, g(x)=f(x)^2, h(x)=g(x)*f(x)-3.  If we try to find a root of h(x), that should still be treated as a polynomial.
#' 4) Be able to handle functions that are defined as the output of some other, built in, function.  For instance, model=fitModel(...), we should be able to find roots of model.
#' 5) Be able to solve systems of equations, in particular with regard to multivariable optimization and Lagrange multipliers.
#'
#' The strategy is:
#' If the expression is single-variable:
#' 1) Treat expr as a polynomial, and use sympy to solve.  If this succeeds, then we are done.
#' 2) If treating as a polynomial fails, then do a traditional multi-start numerical solve.  This is based on a Broyden (Newton-type) method.  Unlike a bracketing method like uniroot,
#' Newton-type methods can find repeated roots.  The disadvantage is that the make no guarantee about absolutely anything. However, we can state at an absolute ton of points, and it
#' seems to work reasonably well.
#' 3) If forced to do so, use the general sympy solver before the numerical one.  This has some advantages, for instance it can solve rational functions, etc., exactly.  However,
#' the general solve command is known to miss some roots and *not* generate an error.  For instance, only two roots are found for sin(x).  Further, it might generate things like asin(2)
#' as an answer, which is fine as a complex number, but by default R will generate NA when converting this to floating point.  This is not an error, per se, so right now we return
#' a table full of N/A when trying to solve sin(x)+2=0.
#' The solver is based on a Broyden method, not a bracketing technique, so it will find repeated roots.
#'
#' If the expression is multi-variable:
#' 1) Attempt to use a the sympy symbolic solver.  This might miss solutions in some cases. It performs well on multinomial and similar problems
#' 2) If symbolic fails, or if forced to do so, try a numerical solver based on a Broyden method.
#'
#' @inheritParams mosaic::findZeros
#' @param trySymbolicSingleVar if TRUE, try using a symbolic solver for solving the single-variable equation rather than skipping right to numeric.
#' @param forceMultivariableNumeric if TRUE, skip the symbolic solver for your multivariable system and go right to numeric.
#' @param verbose if TRUE, print out information about progress.
#' @param roundDigits if x1 and x2 are found numerically, if x1 and x2 round to the same
#' number at roundDigits place, then treat x1 and x2 as the same root.
#' @importFrom grDevices col2rgb rgb topo.colors
#' @importFrom stats aggregate dist na.exclude rnorm
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom utils str
#' @importFrom mosaicCore rhs
#' @importFrom mosaicCore lhs
#' @importFrom magrittr %>%
#' @examples
#' # solve sin(x)=cos(x)
#' findZeros(sin(x)-cos(x)~x)
#'
#' @export
findZeros=function(expr, ..., xlim = c(near - within, near + within),
                           near = 0, within = Inf, nearest = 10, npts = 1000, iterate = 1,
                           sortBy = c("byx", "byy", "radial"),trySymbolicSingleVar=FALSE,forceMultivariableNumeric=FALSE,
                          verbose=FALSE,roundDigits=5){

  #first things first, get the dots arguments.  This is for the purpose of calling the original findZeros,
  #should the sympy method fail.
  dots=list(...)

  #Find all the variables in the expression.
  varNames=all.vars(mosaicCore::rhs(expr))
  #Ok, make a function which takes an expression and returns a corresponding function.
  expression2function=function(expr,variables,...){
    dots=list(...);
    f=function(xvec){
      mydots=dots;
      for (i in 1:length(variables)){
        mydots[[variables[[i]]]]=xvec[[i]];
      }
      return(eval(expr,envir=mydots,enclos=parent.frame()));
    }
    return(f);
  }

  #Let's make a function for our expression(s):
  pfun=expression2function(mosaicCore::lhs(expr),varNames)

  #This will extract a numeric R value from a caracas number
  #There are two options.  First, convert the symbolic number to an exact string, something like sqrt(2)+3i, then have R evaluate that number.
  #That might fail sometimes.  For instance, if you are solving a complicated system then caracas might not know that the number is complex, and write something like
  #sqrt(sqrt(4)-5), for instance.  Evaluating this in R produces NaN, because R doesn't know that we want it as a complex number.
  #The second option is to let python evaluate it to a double, and then get that value back into R.  The only issue with that is that as.character() returns complex numbers
  #with I, not i. R likes complex numbers written as 3i, not 3*I.  So I replace I with 1.0i in the string and then evaluate, which seems to work.
  #Hopefully by using one strategy or antother we can always return a nice R number.
  extractNumberFromSymbolic=function (symbolicNumber){
    method1=function (){
      myz=eval(parse(text=caracas::as_character(symbolicNumber)));
      return(myz)
    }
    method2=function (){
      myz=as.character(symbolicNumber$pyobj$evalf());
      myz=gsub("I","1.0i",myz);
      return(eval(parse(text=myz)))
    }
    tryCatch(return(method1()),
             warning=function(e) {return(method2())},
             error=function (e) {return(method2())})

  }

  #Polynomial solver. Rather than seeing if we have a polynomial and then
  #running a polynomial solver, I suppose that the modern way to do things is just call the solver
  #and then hope for the best.
  #
  #If the function is polynomial then it is my expectation that this will work, possibly returning complex numbers.
  #There should be no case in which this will return, e.g., null.
  polynomialSolver=function(func,symbol,symbolName){
    if (verbose){
      print('Attempting polynomial solve')
    }

    sympy_roots=caracas::sympy_func(func(symbol),"nroots");
    df=data.frame(unlist(lapply(sympy_roots,extractNumberFromSymbolic)));
    names(df)=symbolName;
    return(df);
  }


  #Ok, attempt a generic solve from sympy. THIS IS DANGEROUS AND MIGHT MISS ROOTS!! For example, it misses solutions to sin(x)=0.

  #IT IS UNCLEAR THAT WE SHOULD USE THIS.

  #We might consider using the solveset function in sympy. It promises to return a set containing *all*, even infinitely many, solutions
  #should that be required.  However, at the moment it seems difficult to actually access those solutions and turn them into numbers.

  #This function may return no solutions!

  #This function may also return NaN for a solution -- I should probably just delete those.

  singleVariableSymbolicSolver=function(func,symbol,symbolName){
    if (verbose){
      print("Attempting single-variable symbolic solve")
    }
    sympy_roots=caracas::solve_sys(func(symbol),symbol)
    df=data.frame(unlist(lapply(sympy_roots, function(x){ extractNumberFromSymbolic(x[[symbolName]])})))
    names(df)=symbolName;
    return(df)
  }


  #This is a symbolic solver for systems of equations.  It relies on the symbolic solver from sympy.
  #As such, it might miss solutions, just like the symbolic general solver.  Indeed, there is almost zero difference between the two.
  #We might consider a numerical system solver as an alternative, but we'll give this a go for now.
  symbolicSystemSolver=function(func,symbols,symbolNames){
    if (verbose){
      print("Attempting symbolic system solve")
    }
    sympy_solns=caracas::solve_sys(func(symbols),symbolNames)
    df=data.frame(matrix(NA,nrow=length(sympy_solns),ncol=length(varNames)))
    names(df)=symbolNames;
    for (v in names(df)){
      df[[v]] =  unlist(lapply(sympy_solns,function(x){extractNumberFromSymbolic(x[[v]])}))
    }
    return(df);
  }


  #Ok, this is an attempt at a default solver.  We'll do a multi-start numerical method, the same as the original findZeros did.
  #We'll base the actual solver on sympy's nsolve, as opposed to R's uniroot, because nsolve will find repeated roots.  We do, of course, give up
  #any promise of finding roots even if they are simple, however, we'll retain the ability to to "zoom in" on a region and will be able to find
  #any roots so long as we zoom in close enough.

  #More concretely:
  #attempt to solve
  #expr=0 using a multi-start method.

  #Most of the original findZeros options are basically respected.  In particular, if you think there is a solution near x0 plus or minus 10, then use near=x0 and within=10.
  #Alternatively, if you think that there are solutions in the interval (20,25) then use within=c(20,25).
  #If you want more starting points for your search, which will be slower but have a better chance at finding your root, use npts=whatever you like.  Default is 1000 which seems like
  #a lot already.
  #orderBy, iterate,nearest are ignored.
  #
  #This function does not promise to return solutions within the specified interval, it promises to start the search within the specified interval.
  #
  #Solutions are currently rounded to 5 decimal places -- can think about the correct way to do that later.
  #

  singleVariableNumericSolver=function(func,varName){

    if (verbose){
      print("Attempting single-variable numeric solve")
    }

    #Make sure the provided limits are ok.
    tryCatch(xlim <- range(xlim), error = function(e) stop(paste("Improper limits value --",e)))
    if (xlim[1] >= xlim[2])
      stop(paste("Left limit (", xlim[1], ") must be less than right limit (",
                 xlim[2], ")."))

    #internal.near is the center of where we will be looking.
    internal.near <- near
    if (internal.near < xlim[1] || internal.near > xlim[2]) {
      internal.near <- mean(xlim[1], xlim[2])
    }

    #Mx is max distance from internal.near to the edge of the serach interval. mx and mn will be different only if internal.near is for some reason not set as the
    #center of the search interval.
    mx <- max(xlim - internal.near)
    #min distance to edge of search interval.
    mn <- min(xlim - internal.near)
    if (mx < 0)
      stop("Bug alert: near outside search interval.")
    if (mn > 0)
      stop("Bug alert: near outside search interval.")

    #Set up "big" numbers.
    verybig <- .Machine$double.xmax*0.95
    plainbig <- npts^(0.75)

    #mx might currently be Inf.  Set mx to be either mx (if it is finite) or verybig.  I'm not sure why the max(...,-verybig) is there. I think it is impossible that
    #mx or verybig could be negative, see check above, so comparing to -verybig seems silly.  Keep it for compatibility?
    mx <- max(min(verybig, mx), -verybig)
    mn <- min(max(-verybig, mn), verybig)

    rightseq <- NULL
    leftseq <- NULL
    middleseq <- NULL

    #Ok, divide up our search into three "sections".  [-mx,-plainbig], [-plainbig,plainbig], [plainbig,mx]
    #the middle part likely has the root, so use npts evenly spaced starting points.
    #The two edges likely don't, so use npts exponentially increasing starting points, so that initial points are more bunched up near the center interval.
    if (mx > plainbig) {
      rightseq <- exp(seq(log(max(plainbig, mn)), log(mx),
                          length = npts))
    }
    middleseq <- seq(max(-plainbig, mn), min(plainbig, mx),
                     length = npts)
    if (mn < -plainbig) {
      leftseq <- -exp(seq(log(-mn), log(-min(-plainbig, mx)),
                          length = npts))
    }

    #Finally, shift the search space to be centered at internal.near, make them unique, and sort.
    searchx <- sort(unique(internal.near + c(0, leftseq, middleseq,
                                             rightseq)))

    #Perform a search for roots from each starting point in searchx.
    #If nleqslv throws an error, return NA.  If nleqslv returns 2 (that means it stalled),
    #check the function value and make sure it is small. Note that nleqslv has an option to take multiple
    #starting values.  However, this is working ok.  The only reason I see to change would be performance.
    roots=unlist(lapply(searchx,function(z) {
      tryCatch({
        mysoln=nleqslv::nleqslv(z,func)
        if (mysoln$termcd==1){
          return(mysoln$x)
        }
        else if (mysoln$termcd==2 & abs(pfun(mysoln$x))<1e-10){
          return(mysoln$x)
        }
        else{
          return(NA)
        }
      },
      error=function (w) { return(NA)}
      )

    }))


    `%>%` <- magrittr::`%>%`

    #Make a data frame to hold the roots.
    roots=data.frame(roots)

    #Exclude NAs, and exlude solutions that are outside of the original window.
    roots=roots %>% na.exclude() %>% dplyr::filter(roots>=xlim[[1]]) %>% dplyr::filter(roots<=xlim[[2]]);

    #if the dataframe is not empty get rid of duplicate solutions.
    #Solutions are considered duplicate if they round to the same number at 5 decimal places.
    #The average of all duplicate solutions is considered the winner.
    if (dim(roots)[[1]]>0){
      roots = roots %>% mutate(rounded_value=round(roots,roundDigits))
      roots= roots %>% aggregate(by=list(roots$rounded_value),FUN=mean)
      roots=data.frame(roots$roots)
    }

    #Finally, stick the correct column name on there.
    names(roots)=varName;
    return(roots);
  }


  multiVariableNumericSolver=function(func,varnames,near=rep(0,length(varnames)),within=Inf,npts=1000){

    numVars=length(varnames);

      #if near happens to be a scalar, fix it for the user.
      if (length(near)==1){
        near=rep(near,numVars);
      }

      #Check to make sure that near is of the right dimension.
      if (length(near) != numVars){
        stop(paste("Length of initial guess (",length(near),") must match number of variables to solve for (",numVars,")."))
      }

      #if within is a scalar, take it to mean that within is the same for each dimension
      if (length(within)==1){
        within=rep(within[[1]],numVars);
      }

      #Ok, generate the limits for each dimension.
      xtop=near+within;
      xbottom=near-within;

      #Replace infinities with large numbers.
      verybig=1e6;

      xtop[xtop==Inf]=verybig;
      xbottom[xbottom==Inf]=verybig;
      xtop[xtop==-Inf]=-verybig;
      xbottom[xbottom==-Inf]=-verybig;

      #Within might have been Inf.  Let's just update that to the practical within.
      within=abs(xtop-near);

      #Create a matrix of initial guesses.  One row per guess.
      #Each variable is normally distributed with mean at "near" and standard deviation within[[i]]/3
      x0=matrix(0,npts,numVars);

      for (i in 1:numVars){
        x0[,i]=rnorm(npts,near[[i]],within[[i]]/2)
      }

      #Where are the initial guesses?
      #if (numVars==2){
      #  plot(data.frame(x0))
      #  abline(v=xtop[[1]])
      #  abline(v=xbottom[[1]])
      #  abline(h=xtop[[2]])
      #  abline(h=xbottom[[2]])
      #}

      #Ok, now use searchZeros to run a search from each starting point.
      nlsolv.solns=nleqslv::searchZeros(x0,func,digits=roundDigits);
      solns=data.frame(data=matrix(nrow=0,ncol=numVars));

      solns=rbind(solns,nlsolv.solns$x);

      `%>%` <- magrittr::`%>%`

      #Filter out solutions that we don't want.
      for (i in 1:numVars){
        solns = solns %>% filter(solns[,i]>=xbottom[[i]])
        solns = solns %>% filter(solns[,i]<=xtop[[i]])
      }

      names(solns)=varNames;
      return(solns)



  }

  #This is what to do if the symbolic general solver fails.  Report the error,
  #then call the numeric solver.
  singleVariableSolverErrorHandler= function(e){
    if (verbose){
      print("   Error in symbolic solver, sending to numeric solver.")
    }
    return(singleVariableNumericSolver(pfun,varNames[[1]]))
  }



  #This defines what to do if the symbolic polynomial solver fails.
  #if forced to do so, try the symbolic general solver.
  #If that fails, or if it didn't have to do so in the first place, then do the numeric solver.
  polynomialErrorHandler=function(e){
    if (trySymbolicSingleVar){
      if (verbose){
        print("   Error in polynomial solver. Sending to single-variable general symbolic solver")
      }
      tryCatch(
        return(singleVariableSymbolicSolver(pfun,varSymbols,varNames)),
        warning=singleVariableSolverErrorHandler,
        error=singleVariableSolverErrorHandler)
    }
    else{
      if (verbose){
        print("   Error in polynomial solver. Sending to numeric solver.")
      }
      return(singleVariableNumericSolver(pfun,varNames[[1]]))
    }
  }











  #Here is the actual solving.  Almost all of the above is defining how to do things, not doing them.

  #Make a sympy symbol for each variable.
  varSymbols=lapply(varNames,function(x) {return(caracas::symbol(x))})

  numberOfVariables=length(varNames);
  isMultiVariable=numberOfVariables>1;


  #If the system is multivariable, first try the symbolic solver, then numeric. go straight to numeric if forceMultivariableNumeric=TRUE
  if (isMultiVariable){
    if (forceMultivariableNumeric){
      print("Attempting multivariable numeric solve -- forced to do so")
      return(multiVariableNumericSolver(pfun,varNames,near,within,npts))
    }
    else{
      tryCatch({
        return(symbolicSystemSolver(pfun,varSymbols,varNames))
      },
      error=function(e){
        print("multivariable symbolic solve failed. Attempting numeric solve.")
        return(multiVariableNumericSolver(pfun,varNames,near,within,npts))
      })
    }
  }

  #If the system is single-variable, first try the polynomial solver. If it fails, go to the default numeric solver unless
  #trySymbolicSingleVariable=T.
  else{
    tryCatch(
      return(polynomialSolver(pfun,varSymbols,varNames)),
      warning=polynomialErrorHandler,
      error=polynomialErrorHandler)
  }
}



