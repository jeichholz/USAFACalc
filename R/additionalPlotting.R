#' Plots a filled region defined by one or two functions.
#' @param expression  an expression defining the first function.
#' @param bottom defined the "bottom" of the shaded region. If you give a
#' number, then the bottom is just the number. You may also give an expression
#' defining a second function to shade the area between two curves.
#' @param col the color of the shaded region
#' @param alpha the transparency of the shaded region
#' @param add whether to add the solid region to current figure or start anew
#' @param xlim the region above which to draw
#' @param addoutline whether or not to graph the one/two curves in a highlight black color to make a nice outline.
#' @param ... additional options which are passed to lattice::xyplot
#' @param plot The figure to which this region should be added.
#' @examples
#' f=mosaic::makeFun(sin(x)~x)
#' plotFunFill(f(x)~x,xlim=c(-3,3))
#' grid.on(col="grey")
#' @examples
#' f=mosaic::makeFun(sin(x)~x)
#' plotFunFill(f(x)~x,bottom=-2,xlim=c(-3,3))
#' grid.on()
#' @examples
#' f=mosaic::makeFun(sin(x)~x)
#' g=mosaic::makeFun(-x-1/2*cos(5*x)~x)
#' plotFunFill(f(x)~x,bottom=g(x)~x,xlim=c(0,3))
#' grid.on()
#'
#' @export
plotFunFill=function(expression,bottom=0,xlim=c(0,1),col="red",alpha=0.5,add=FALSE, addoutline=TRUE, plot=lattice::trellis.last.object(),...){
  dots=list(...)

  #browser()
  f=mosaic::makeFun(expression)
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


  #browser()
  #trellis.last.object()+latticeExtra::panel.xyarea(xs,ys,origin=0,data=c(xs=xs,ys=ys))

  var=deparse(mosaicCore::rhs(expression))
  fnname=deparse(mosaicCore::lhs(expression))

  if(add==FALSE){
    A=lattice::xyplot(y ~ x, data = data.frame(x=xs,y=ys),
                      panel = lattice::panel.polygon, xlim=xlim, rule = "none", col=col,alpha=alpha,xlab=var,
                      ylab=fnname,...)
  }
  else{
    #A=mosaic::ladd(lattice::panel.polygon(xs,ys,rule = "none", col=col,alpha=alpha,dots),
    #               data=list(xs=xs,ys=ys,col=col,alpha=alpha,dots=list(...)),plot=plot)
    A=plot+latticeExtra::layer(lattice::panel.polygon(xs,ys,rule = "none", border="black",col=col,alpha=alpha,dots),data=list(xs=xs,ys=ys,col=col,alpha=alpha,dots=list(...)))
  }

  if (addoutline){
    xs=seq(xlim[[1]],xlim[[2]],length.out=N);
    #browser()
    ys=f(xs)
    A=A+latticeExtra::layer(lattice::panel.xyplot(xs,ys,lwd=2,col="black",type="l"),data=list(xs=xs,ys=ys))
    if (is.numeric(bottom)){
      xs=xlim
      ys=c(bottom,bottom)
    }
    else{
      ys=g(xs)
    }
   A=A+latticeExtra::layer(lattice::panel.xyplot(xs,ys,col="black",type="l",lwd=2),data=list(xs=xs,ys=ys));

  }
  return(A)
}


#' Places text on a figure.
#' @param text The message to display
#' @param x The x-coordinate
#' @param y The y-coordinate
#' @param col the text color
#' @param plot the figure to add the text to.
#' @param zoom the factor by which to expand the text.
#' @param ... additional options that are passed to panel.text.
#' @examples
#' plotFunFill(x^2~x)
#' place.text("x^2",0.5,0.8)
#'
#' @export
place.text=function(text,x,y,col="black",zoom=1,plot=lattice::trellis.last.object(), ...){
  dots=list(...)
  cex=zoom
  return(mosaic::ladd(lattice::panel.text(x=x,y=y,labels=text,alpha=1,cex=cex,col=col),data=list(x=x,y=y,text=text,col=col,cex=cex),plot=plot))
}

#' Turns on a grid for easier viewing.
#' @param col grid color
#' @param lty line type, 1,2,or 3 are good options.
#' @param lwd line width, 2 is about right.
#' @param h # of horizontal grid lines to make, with -1 meaning to just line up to tickmarks.
#' @param v # of vertical grid lines to make, with -1 meaning to just line up to tickmarks.
#' @param plot the figure to which to add the grid.
#' #@examples
#' plotFunFill(x^2~x)
#' grid.on()
#'
#' @export
grid.on=function(h=-1,v=-1,lwd=2,col="black",lty=2,plot=lattice::trellis.last.object()){

  return(mosaic::ladd(lattice::panel.grid(h=h,v=v,lty=lty,col=col),data=list(h=h,v=v,lty=lty,col=col),plot=plot))
}


#' Plots "traditional" x-y axes.
#' @param col the color of the axes
#' @param lwd the line width of the axes
#' @param xat the y-coordinate where the axis will be drawn.
#' @param yat the x-coordinate where the axis will be drawn.
#' @param plot the figure to which to add axes.
#' @param ... additional options that will be passed to
#' @examples
#' f=mosaic::makeFun(sin(x)~x)
#' plotFunFill(f(x)~x,xlim=c(-3,3))
#' mathaxis.on()
#' @export
mathaxis.on=function(lwd=3,col="black",xat=0,yat=0,plot=lattice::trellis.last.object(),...){
  #browser()
  #browser()
  A=mosaic::ladd(lattice::panel.abline(h=xat,col=col,lwd=lwd,dots),data=list(xat=xat,lwd=lwd,col=col,dots=list(...)),plot=plot)
  A=mosaic::ladd(lattice::panel.abline(v=yat,lwd=lwd,col=col,dots),data=list(yat=yat,lwd=lwd,col=col,dots=list(...)),plot=A)
  return(A)
}



