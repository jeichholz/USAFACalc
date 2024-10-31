#' Plot a vector field
#' Scales all vectors to have same length
#' @param expression A function of two variables that returns two ouputs to plot.
#' @param xlim the limits of the horizontal axis
#' @param ylim the limits of the vertical axis
#' @param grid.by the spacing between arrow bases
#' @param radius the size of the arrow heads
#' @param lwd arrow line width
#' @param lty arrow line type
#' @param col arrow color
#' @param add T/F variable indicating whether this should be added to plot
#' @param plot the plot to which this should be added, default is the current plot.
#' @examples
#' # example code
#' f=mosaic::makeFun(c(-y,x)~x&y)
#' plotVectorField(f(x,y)~x&y, xlim=c(-3,3),ylim=c(-4,4))
#'
#' #You can set the color, line width, or line type as well.
#' plotVectorField(f(x,y)~x&y, xlim=c(-3,3),ylim=c(-4,4),col="green",lwd=3,lty=1)
#'
#' #You can add it to other plots if you so choose.
#' plotFun(sin(x)~x,xlim=c(0,5))
#' plotVectorField(f(x,y)~x&y, xlim=c(-3,3),ylim=c(-4,4),col="green",lwd=3,lty=1,add=TRUE)
#' @export
plotVectorField = function(expression,xlim=c(-5,5),ylim=c(-5,5),
                            grid.by = (xlim[[2]]-xlim[[1]])/19,radius=grid.by*0.8,col=NA,lwd=1,lty=1,add=FALSE,smooth=FALSE, plot = lattice::trellis.last.object()){

  #expression should be an expression which takes as input two variables and returns a list of length 2 as output.

  allVars=all.vars(mosaic::rhs(expression))

  if (length(allVars)!=2){
    stop(paste("Vector function must take two variables as input.  You supplied ",as.character(length(allVars)),
               "(",paste0(as.character(allVars),collapse=","),")"))
  }

  FUN=mosaic::makeFun(expression);

  #plot.new();
  #frame();
  #plot(xlim,ylim,main="Direction field", ylab = as.character(allVars[[2]]), xlab = as.character(allVars[[1]]), pch = ".")

  # grid points
  seqx = seq(xlim[[1]],xlim[[2]],grid.by)
  seqy = seq(ylim[[1]],ylim[[2]],grid.by)

  grid=expand.grid(seqx,seqy);
  grid$toVar1=0;
  grid$toVar2=0;
  grid$Color=col;

  redfunc= function(theta) {floor(255*bump(theta,x0=0,r=2*pi/3))+floor(255*bump(theta,x0=2*pi,r=2*pi/3))};
  greenfunc=function(theta) {floor(255*bump(theta,x0=2*pi/3,r=2*pi/3))};
  bluefunc=function(theta) {floor(255*bump(theta,x0=4*pi/3,r=2*pi/3))};


  for (i in 1:dim(grid)[[1]]){
    x=grid[i,]$Var1;
    y=grid[i,]$Var2;
    vec=FUN(x,y);
    nrmvec=sqrt(mosaic::dot(vec,vec));
    if (nrmvec>0){
      vecN=vec/nrmvec;
    }
    else{
      vecN=vec;
    }

    slope=vec[[2]]/vec[[1]];


    if (is.na(grid[i,]$Color)){
      if (smooth==TRUE){
        theta=atan(slope);
        if (theta<0){theta=2*pi+theta};
        if (vec[[1]]<0){
          theta=theta-pi;
        }
        cor=ggtern::rgb2hex(redfunc(theta),greenfunc(theta),bluefunc(theta));
        cat(paste(theta,"\n"))
        grid[i,]$Color=cor;
      }
      else if(is.na(slope)){
                cor = "black"
      } else if(slope > 0){
                cor = "blue"
      }else if (slope < 0) {
                cor = "red"
              }
      else if(slope == 0) {
                cor = "green"
      }
      grid[i,]$Color=cor;
    }
    grid[i,]$toVar1=x+radius*vecN[[1]]
    grid[i,]$toVar2=y+radius*vecN[[2]]
  }

  under=FALSE;

  if (!add){
    print(lattice::xyplot(Var2~Var1,panel=function(...){
      lattice::panel.arrows(x0=grid$Var1,y0=grid$Var2,x1=grid$toVar1,y1=grid$toVar2,length=0.2*radius,col=grid$Color,lwd=lwd,lty=lty)

    },data=grid,xlab=allVars[[1]],ylab=allVars[[2]]))
  }

  if (add){
     # print(plot + latticeExtra::layer(do.call(panel.arrows, list(grid$Var1, grid$Var2,grid$toVar1,grid$toVar2)),
    #                             data = as.list(environment()), under = under))

     print(plot + latticeExtra::layer(do.call(lattice::panel.arrows, list(grid$Var1, grid$Var2,grid$toVar1,grid$toVar2,col=grid$Color,length=0.2*radius,lwd=lwd,lty=lty)),
                                 data = as.list(environment()), under = under))
  }






}

#' Plot a direction field for a scalar ODE
#' @param expression an expression giving the right hand side of the ODE, as a function of t and y, in that order.
#' @param tlim the extent of the horizonatal axis
#' @param ylim the extend of the vertical axis
#' @inheritParams plotVectorField
#' @param y0s If desired, a list of initial values from which to draw trajectories.
#' @examples
#'  #Consider the ODE
#'  #y'=1/2*y+cos(t)
#'  #Just give the plotter the right hand side with an expression as usual.  You may omit t
#'  if the equation is autnomous. If you include t, it *must* be the first argument.
#'  plotODEDirectionField(1/2*y+cos(t)~t&y)
#'
#'
#'
#'  #You can set the color, line width, and type as well.
#'  plotODEDirectionField(1/2*y+cos(t)~t&y,col="black",lwd=2)
#'
#'  #You can add it to another plot, or another plot to it.
#'  results=Euler(1/2*y+cos(t)~t&y,tlim=c(0,10),y0=0,stepSize=0.1)
#'  mosaic::plotPoints(y~t,data=results,add=TRUE)
#'
#'  plotODEDirectionField(-y*(y-1)~y)
#'
#'  #You can also just plot a vector field and add on trajectories from different initial
#'  #conditions using the y0s option. Add one initial value for every trajectory you want.
#'  plotODEDirectionField(1/2*y+cos(t)~t&y,ics=c(-2,0,-1,2))
#'
#' @export
plotODEDirectionField=function(expression,tlim=c(0,10),ylim=c(-5,5),ics=NA, grid.by=(tlim[[2]]-tlim[[1]])/19,
                               radius=grid.by*0.8,col=NA,lty=1,lwd=1,add=FALSE,plot=lattice::trellis.last.object()){

  y0s=ics;

  allVars=all.vars(mosaic::rhs(expression))

  if (!("t" %in% allVars)){
    exprstr=as.character(expression)

    expression=as.formula(paste(exprstr[[2]],"~t&",exprstr[[3]],collapse=" "))
    allVars=all.vars(mosaic::rhs(expression))
  }

  if (length(allVars)!=2){
    stop(paste("Vector function must take two variables as input.  You supplied ",as.character(length(allVars)),
               "(",paste0(as.character(allVars),collapse=","),")"))
  }

  dydt=mosaicCore::makeFun(expression);
  dydt1=mosaicCore::makeFun(c(1,dydt(t,y))~t&y)

  plotVectorField(dydt1(t,y)~t&y,xlim=tlim,ylim=ylim,radius=radius,grid.by=grid.by,col=col,lty=lty,lwd=lwd,add=add,plot=plot)

  if (!any(is.na(y0s))){
    for (y0 in y0s){
      em=Euler(expression,tlim,y0,0.01);
      print(mosaic::plotPoints(y~t,data=em,add=TRUE))
    }
  }

}

bump <- function(x,x0=0,r=1,h=1){
  y=x*0;
  y=(abs(x-x0)<r-1e-2)*h*exp(1/(((x-x0)/r)^2-1))/exp(-1);
  y[is.na(y)]=0;
  return(y);
}





