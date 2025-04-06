#' Plot a vector field
#' Scales all vectors to have same length
#' @param expression A function of two variables that returns two ouputs to plot.
#' @param xlim the limits of the horizontal axis
#' @param ylim the limits of the vertical axis
#' @param N the number of horizontal and vertical rows of arrows to draw
#' @param normalize if true, then plot all vectors as unit length
#' @param lwd arrow line width
#' @param col arrow color
#' @param add T/F variable indicating whether this should be added to plot
#' @param plot the plot to which this should be added, default is the current plot.
#' @param ... additional arguments to pass to xyplot arrows.
#' @examples
#' # example code
#' f=mosaic::makeFun(c(-y,x)~x&y)
#' plotVectorField(f(x,y)~x&y, xlim=c(-3,3),ylim=c(-4,4))
#'
#' #You can set the color, line width, or line type as well.
#' plotVectorField(f(x,y)~x&y, xlim=c(-3,3),ylim=c(-4,4),col="green",lwd=3)
#'
#' #You can add it to other plots if you so choose.
#' plotFun(sin(x)~x,xlim=c(0,5))
#' plotVectorField(f(x,y)~x&y, xlim=c(-3,3),ylim=c(-4,4),col="green",lwd=3,add=TRUE)
#' @export
plotVectorField = function(expression,xlim=c(-5,5),ylim=c(-5,5),N=20,col="cornflowerblue",lwd=2, normalize=FALSE,
                           add=FALSE,plot = lattice::trellis.last.object(), ...){

  #expression should be an expression which takes as input two variables and returns a list of length 2 as output.

  allVars=all.vars(mosaic::rhs(expression))

  dots=list(...)

  if (length(allVars)!=2){
    stop(paste("Vector function must take two variables as input.  You supplied ",as.character(length(allVars)),
               "(",paste0(as.character(allVars),collapse=","),")"))
  }

  FUN=mosaic::makeFun(expression);

  if (add){
    xlim=plot$x.limits;
    ylim=plot$y.limits;
  }

  # grid points
  seqx = seq(xlim[[1]],xlim[[2]],length.out=N)
  seqy = seq(ylim[[1]],ylim[[2]],length.out=N)

  #radius=0.8*max(seqx[[2]]-seqx[[1]],seqy[[2]]-seqy[[1]])

  radius_scrunits=0.8*1/N;

  grid=expand.grid(seqx,seqy);
  grid$Fx=0;
  grid$Fy=0;
  grid$Color=col;

  for (i in 1:dim(grid)[[1]]){
    x=grid[i,]$Var1;
    y=grid[i,]$Var2;
    vec=FUN(x,y);
    grid$Fx[i]=vec[[1]];
    grid$Fy[i]=vec[[2]];
  }

  Lx=xlim[[2]]-xlim[[1]];
  Ly=ylim[[2]]-ylim[[1]];

  grid$nrm=sqrt(grid$Fx^2+grid$Fy^2)
  grid$nrm_scrunits=sqrt((grid$Fx/Lx)^2+(grid$Fy/Ly)^2)
  maxrootlen_scrunits=max(sqrt(grid$nrm_scrunits));

  if (normalize){
    grid$displayLen_scrunits=radius_scrunits;
  }
  else{
    grid$displayLen_scrunits=radius_scrunits*sqrt(grid$nrm_scrunits)/maxrootlen_scrunits
  }

  grid$toVar1=grid$Var1+ifelse(grid$nrm==0,0,grid$displayLen_scrunits/grid$nrm_scrunits)*grid$Fx
  grid$toVar2=grid$Var2+ifelse(grid$nrm==0,0,grid$displayLen_scrunits/grid$nrm_scrunits)*grid$Fy

  #browser()
  if (!add){
    return(lattice::xyplot(Var2~Var1,
                           panel=function(...){
                              lattice::panel.arrows(x0=grid$Var1,
                                                    y0=grid$Var2,
                                                    x1=grid$toVar1,
                                                    y1=grid$toVar2,
                                                    length=grid::unit(0.3*grid$displayLen,"npc"),
                                                    col=grid$Color,lwd=lwd,
                                                    xlab=xlab,ylab=ylab,...)},
                          data=grid,
                          xlab=allVars[[1]],ylab=allVars[[2]],
                          xlim=xlim,ylim=ylim,...))
  }

  if (add){
    under=FALSE;
     return(plot + latticeExtra::layer(do.call(lattice::panel.arrows,
                                               list(grid$Var1, grid$Var2,grid$toVar1,grid$toVar2,col=grid$Color,length=grid::unit(0.3*grid$displayLen,"npc"),
                                                    lwd=lwd)),data = as.list(environment()), under = under))
  }






}

#' Plot a direction field for a scalar ODE
#' @param expression an expression giving the right hand side of the ODE, as a function of t and y, in that order.
#' @param tlim the extent of the horizonatal axis
#' @param ylim the extend of the vertical axis
#' @inheritParams plotVectorField
#' @param ics If desired, a list of initial values from which to draw trajectories.
#' @examples
#'  #Consider the ODE
#'  #y'=1/2*y+cos(t)
#'  #Just give the plotter the right hand side with an expression as usual.  You may omit t
#'  #if the equation is autnomous. If you include t, it *must* be the first argument.
#'  plotODEDirectionField(1/2*y+cos(t)~t&y)
#'
#'
#'
#'  #You can set the color, line width, and type as well.
#'  plotODEDirectionField(1/2*y+cos(t)~t&y,col="black",lwd=2)
#'
#'  #You can add it to another plot, or another plot to it.
#'  results=Euler(1/2*y+cos(t)~t&y,tlim=c(0,10),ic=0,stepSize=0.1)
#'  mosaic::plotPoints(y~t,data=results,add=TRUE)
#'
#'  #You can give an autonomous equation and not specify t as an input.
#'  plotODEDirectionField(-y*(y-1)~y,tlim=c(0,4),ylim=c(-0.1,2.1),N=8)
#'
#'  #You can also just plot a vector field and add on trajectories from different initial
#'  #conditions using the y0s option. Add one initial value for every trajectory you want.
#'  plotODEDirectionField(1/2*y+cos(t)~t&y,ics=c(-2,0,-1,2))
#'
#' @export
plotODEDirectionField=function(expression,tlim=c(0,10),ylim=c(-5,5),ics=NA, N=20,
                               col="black",lwd=2,add=FALSE,plot=lattice::trellis.last.object(),...){

  y0s=ics;

  allVars=all.vars(mosaic::rhs(expression))

  if (!("t" %in% allVars)){
    exprstr=as.character(expression)

    expression=stats::as.formula(paste(exprstr[[2]],"~t&",exprstr[[3]],collapse=" "))
    allVars=all.vars(mosaic::rhs(expression))
  }

  if (length(allVars)!=2){
    stop(paste("Vector function must take two variables as input.  You supplied ",as.character(length(allVars)),
               "(",paste0(as.character(allVars),collapse=","),")"))
  }

  exprstr=as.character(expression)
  dydt1=stats::as.formula(paste("c(1,",exprstr[[2]],")~",exprstr[[3]],collapse=" "))

  A=plotVectorField(dydt1,xlim=tlim,ylim=ylim,N=N,col=col,lwd=lwd,add=add,normalize=TRUE,plot=plot,...)

  if (!any(is.na(y0s))){
    dydt=mosaicCore::makeFun(expression);
    odefun=function(t,y,parms){
      return(list(dydt(t,y)))
    }
    for (y0 in y0s){
      #browser()
      soln=as.data.frame(deSolve::ode(y0,seq(tlim[[1]],tlim[[2]],length.out=1000),odefun));
      colnames(soln)[[2]]="y"
      A=mosaic::plotPoints(y~time,data=soln,add=TRUE,type="l",lwd=lwd+2,col="darkorchid",plot=A)
    }
  }
  return(A)

}


#' Plot the phase plane for a vector field.
#' @param ddt A vector-valued function that defines the RHS of your system of ODEs.
#' @param xlim the limits of the horizontal axis.  Can also use a name that matches your variables.
#' @param ylim the limits of the vertical axis. Can also use a name that matches your variables.
#' @param N the number of horizontal and vertical rows of arrows to draw
#' @param lwd arrow line width
#' @param col arrow color
#' @param add T/F variable indicating whether this should be added to plot
#' @param plot the plot to which this should be added, default is the current plot.
#' @param ics -- a list of initial conditions from which to draw trajectories.
#' @param ... additional arguments to pass to xyplot arrows.
#' @examples
#' # example code
#' f=mosaic::makeFun(c(-y,x)~x&y)
#' plotPhasePlane(f(x,y)~x&y, xlim=c(-3,3),ylim=c(-4,4))
#'
#' #Classic predator-prey.
#' ode=mosaic::makeFun(c(2*rabbits-0.5*rabbits*foxes, -1*foxes+0.25*rabbits*foxes)~rabbits&foxes)
#' plotPhasePlane(ode(rabbits,foxes)~rabbits&foxes,rabbitslim=c(0,10),foxeslim=c(0,10),ics=c(4,2))
#' mosaic::plotPoints(c(0,4)~c(0,4),pch=19,cex=1.3,col="magenta",add=TRUE)
#' plotPhasePlane(ode(rabbits,foxes)~rabbits&foxes,rabbitslim=c(0,10),foxeslim=c(0,10),ics=list(c(4,2),c(4,3)))


#' @export

plotPhasePlane<-function(ddt,xlim=c(-5,5),ylim=c(-5,5),add=FALSE, N=20,
                         col="cornflowerblue",lwd=2,plot=lattice::trellis.last.object(),ics=NA,...){


  allVars=all.vars(mosaic::rhs(ddt))
  #If there are three variables, assume the first one is time, and take that variable out.
  if (length(allVars)==3 && "t" %in% allVars){
    exprstr=as.character(ddt)
    ddt=stats::as.formula(paste0(exprstr[[2]],"~",allVars[[2]],"&",allVars[[3]],collapse=" "))
    allVars=allVars[2:3]
  }
  dots=list(...)
  lims <- mosaic::inferArgs(dots = dots, vars = allVars, defaults = list(xlim = xlim,ylim=ylim))
  xlim=lims$xlim
  ylim=lims$ylim
  #browser()
  A=plotVectorField(ddt,xlim=xlim,ylim=ylim,add=add,N=N,col=col,lwd=lwd,plot=plot,normalize=TRUE,...)

  if (!any(is.na(ics))){
    #browser()
    fun1=mosaic::makeFun(ddt)
    odefun=function(t,y,parms){
      return(list(fun1(y[[1]],y[[2]])))
    }

    times=seq(0,100,by=0.1)
    if (is(ics,"numeric")){
      ics=list(ics)
    }

    for (ic in ics){

      if (length(ic)!=2){
        print(paste("Error. Length of initial condition must be 2. Given initial condition (",ic,") is",length(ic),collapse = " "))
      }

      soln=as.data.frame(deSolve::ode(ic,times,odefun))
      colnames(soln)[[2]]="x"
      colnames(soln)[[3]]="y"
      A=mosaic::plotPoints(y~x,data=soln,plot=A,add=TRUE,col="firebrick",type="l",lwd=lwd+1)
    }


  }

  return(A)
}



