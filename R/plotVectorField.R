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

  radius=0.8*max(seqx[[2]]-seqx[[1]],seqy[[2]]-seqx[[1]])

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


  grid$nrm=sqrt(grid$Fx^2+grid$Fy^2)
  maxrootlen=max(sqrt(grid$nrm));

  if (normalize){
    grid$displayLen=radius;
  }
  else{
    grid$displayLen=radius*sqrt(grid$nrm)/maxrootlen
  }

  grid$toVar1=grid$Var1+ifelse(grid$nrm==0,0,grid$displayLen/grid$nrm)*grid$Fx
  grid$toVar2=grid$Var2+ifelse(grid$nrm==0,0,grid$displayLen/grid$nrm)*grid$Fy

  if (!add){
    return(lattice::xyplot(Var2~Var1,panel=function(...){
      lattice::panel.arrows(x0=grid$Var1,y0=grid$Var2,x1=grid$toVar1,y1=grid$toVar2,length=grid::unit(0.3*grid$displayLen,"native"),
                            col=grid$Color,lwd=lwd,xlab=xlab,ylab=ylab)},
                            data=grid,xlab=allVars[[1]],ylab=allVars[[2]],...))
  }

  if (add){
    under=FALSE;
     return(plot + latticeExtra::layer(do.call(lattice::panel.arrows,
                                               list(grid$Var1, grid$Var2,grid$toVar1,grid$toVar2,col=grid$Color,length=grid::unit(0.3*grid$displayLen,"native"),
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
                               col="black",lwd=2,add=FALSE,plot=lattice::trellis.last.object()){

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

  dydt=mosaicCore::makeFun(expression);
  dydt1=mosaicCore::makeFun(c(1,dydt(t,y))~t&y)

  A=plotVectorField(dydt1(t,y)~t&y,xlim=tlim,ylim=ylim,N=N,col=col,lwd=lwd,add=add,normalize=TRUE,plot=plot)

  if (!any(is.na(y0s))){
    for (y0 in y0s){
      em=Euler(expression,tlim,y0,0.01);
      A=mosaic::plotPoints(y~t,data=em,add=TRUE,plot=A)
    }
  }
  return(A)

}

bump <- function(x,x0=0,r=1,h=1){
  y=x*0;
  y=(abs(x-x0)<r-1e-2)*h*exp(1/(((x-x0)/r)^2-1))/exp(-1);
  y[is.na(y)]=0;
  return(y);
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
#' @param ... additional arguments to pass to xyplot arrows.
#' @examples
#' # example code
#' f=mosaic::makeFun(c(-y,x)~x&y)
#' plotPhasePlane(f(x,y)~x&y, xlim=c(-3,3),ylim=c(-4,4))
#'
#' #Classic predator-prey.
#' ode=mosaic::makeFun(c(2*rabbits-0.5*rabbits*foxes, -1*foxes+0.25*rabbits*foxes)~rabbits&foxes)
#' plotPhasePlane(ode(rabbits,foxes)~rabbits&foxes,rabbitslim=c(0,10),foxeslim=c(0,10))
#' mosaic::plotPoints(c(0,4)~c(0,4),pch=19,cex=1.3,col="magenta",add=TRUE)
#' soln=USAFACalc::Euler(ode(rabbits,foxes)~rabbits&foxes,ic=c(4,2),tlim=c(0,10),stepSize=0.001)
#' mosaic::plotPoints(foxes~rabbits,data=soln,pch=".",cex=1.3,add=TRUE,col="forestgreen")
#' @export

plotPhasePlane<-function(ddt,xlim=c(-5,5),ylim=c(-5,5),add=FALSE, N=20,
                         col="cornflowerblue",lwd=2,plot=lattice::trellis.last.object(),...){


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

  plotVectorField(ddt,xlim=xlim,ylim=ylim,add=add,N=N,col=col,lwd=lwd,plot=plot,normalize=TRUE)
}



