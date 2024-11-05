#' Plots a filled region defined by one or two functions.
#' @param f  an expression defining the first function.
#' @param bottom defined the "bottom" of the shaded region. If you give a
#' number, then the bottom is just the number. You may also give an expression
#' defining a second function to shade the area between two curves.
#' @param col the color of the shaded region
#' @param alpha the transparency of the shaded region
#' @param add whether to add the solid region to current figure or start anew
#' @param xlim the region above which to draw
#' @param addoutline whether or not to graph the one/two curves in a highlight black color to make a nice outline.
#' @examples
#' f=mosaic::makeFun(sin(x)~x)
#' plotFunFill(f(x)~x,xlim=c(-3,3))
#' #mosaic::plotFun(f(x)~x,col="black",lwd=2,add=TRUE)
#' #mosaic::ladd(panel.grid(h=-1,v=-1,col="grey"))
#' @examples
#' f=mosaic::makeFun(sin(x)~x)
#' plotFunFill(f(x)~x,bottom=-2,xlim=c(-3,3))
#' #mosaic::plotFun(f(x)~x,col="black",lwd=2,add=TRUE)
#' #mosaic::ladd(panel.grid(h=-1,v=-1,col="grey"))
#' @examples
#' f=mosaic::makeFun(sin(x)~x)
#' g=mosaic::makeFun(-x-1/2*cos(5*x)~x)
#' plotFunFill(f(x)~x,bottom=g(x)~x,xlim=c(0,3))
#' #mosaic::plotFun(f(x)~x,col="black",lwd=2,add=TRUE)
#' #mosaic::plotFun(g(x)~x,col="black",lwd=2,add=TRUE)
#' #mosaic::ladd(panel.grid(h=-1,v=-1,col="grey"))
#'
#' @export
plotFunFill=function(f,bottom=0,xlim=c(0,1),col="red",alpha=0.5,add=FALSE, addoutline=TRUE, plot=lattice::trellis.last.object(),...){
  dots=list(...)

  #browser()
  f=mosaic::makeFun(f)
  N=100
  xs=seq(xlim[[1]],xlim[[2]],length.out=N);
  ys=f(xs);
  if (is.numeric(bottom)){
    xs=c(xs,xlim[[2]],xlim[[1]]);
    ys=c(ys,bottom,bottom);
  }
  else{
    g=mosaic::makeFun(bottom);
    xxs=seq(xlim[[2]],xlim[[1]],length=N)
    yys=g(xxs)
    xs=c(xs,xxs)
    ys=c(ys,yys)
    xs=c(xs,xlim[[1]])
    ys=c(ys,f(xlim[[1]]))
  }



  #trellis.last.object()+latticeExtra::panel.xyarea(xs,ys,origin=0,data=c(xs=xs,ys=ys))
  d=data.frame(x=xs,y=ys)
  if(add==FALSE){
    A=lattice::xyplot(y ~ x, data = d, panel = lattice::panel.polygon, xlim=xlim, rule = "none", col=col,alpha=alpha, ...)
  }
  else{
    A=plot+lattice::xyplot(y ~ x, data = d, panel = lattice::panel.polygon,  rule = "none", col=col,alpha=alpha, ...)
  }

  if (addoutline){
    xs=seq(xlim[[1]],xlim[[2]],length.out=N);
    A=A+lattice::xyplot(y~x,data=data.frame(x=xs,y=f(xs)),col="black",type="l",lwd=2)
    if (is.numeric(bottom)){
      ys=c(bottom,bottom)
      A=A+lattice::xyplot(y~x,data=data.frame(x=xlim,y=ys),col="black",type="l",lwd=2);
    }
    else{
      A=A+lattice::xyplot(y~x,data=data.frame(x=xs,y=g(xs)),col="black",type="l",lwd=2)
    }
  }


  return(A)
}


#' Places text on a figure.
#' @param text The message to display
#' @param x The x-coordinate
#' @param y The y-coordinate
#' @param col the text color
#' @examples
#' plotFunFill(x^2~x)
#' place.text("x^2",1,0.8)
#'
#' @export
place.text=function(text,x,y,col="black",...){
  dots=list(...)
  return(mosaic::ladd(lattice::panel.text(x=x,y=y,labels=text,col=col,...),data=list(x=x,y=y,text=text,col=col,dots=dots)))
}

#' Turns on a grid for easier viewing.
#' @param col grid color
#' @param lty line type, 1,2,or 3 are good options.
#' @param lwd line width, 2 is about right.
#' @param h # of horizontal grid lines to make, with -1 meaning to just line up to tickmarks.
#' @param v # of vertical grid lines to make, with -1 meaning to just line up to tickmarks.
#' @examples
#' plotFunFill(x^2~x)
#' grid.on()
#'
#' @export
grid.on=function(h=-1,v=-1,lwd=2,col="black",lty=2,...){
  return(mosaic::ladd(lattice::panel.grid(h=h,v=v,lty=lty,col=col),data=list(h=h,v=v,lty=lty,col=col)))
}


#' Plots "traditional" x-y axes.
#' @param xextent the x-range over which the x-axis will be plotted. If it isn't big enough by default, make it bigger.
#' @param yextent the y-range over which the x-axis will be plotted. If it isn't big enough by default, make it bigger.
#' @param col the color of the axes
#' @param lwd the line width of the axes
#' @param xat the y-coordinate where the axis will be drawn.
#' @param yat the x-coordinate where the axis will be drawn.
#' @examples
#' f=mosaic::makeFun(sin(x)~x)
#' plotFunFill(f(x)~x,xlim=c(-3,3))
#' mathaxis.on()
#' @export
mathaxis.on=function(lwd=3,col="black",xat=0,yat=0,xextent=1e30,yextent=1e30,...){
  #browser()
  xlim=c(-xextent,xextent)
  ylim=c(-yextent,yextent)
  d=data.frame(x=xlim,y=c(xat,xat));
  e=data.frame(x=c(yat,yat),y=c(-yextent,yextent));
  return(lattice::trellis.last.object()+lattice::xyplot(y~x,data=d,panel=lattice::panel.polygon,rule="none",lwd=lwd,col=col,...)+lattice::trellis.last.object()+lattice::xyplot(y~x,data=e,panel=lattice::panel.polygon,rule="none",lwd=lwd,col=col,...))
}



