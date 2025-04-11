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


#' Places text on a figure.  Supports some latex symbols, but backslashes must be double backslashes.
#' @param text The message to display
#' @param x The x-coordinate
#' @param y The y-coordinate
#' @param col the text color
#' @param plot the figure to add the text to.
#' @param zoom the factor by which to expand the text.
#' @param tex if TRUE, try to interpret this string as TeX code.
#' @param ... additional options that are passed to panel.text.
#' @examples
#' plotFunFill(x^2~x)
#' place.text("x^2",0.5,0.8)
#'
#' plotFun(pnorm(x)~x,xlim=c(-3,3),lwd=2)
#' mathaxis.on()
#' plotFunFill(pnorm(x)~x,xlim=c(-3,1.1),col='magenta',add=TRUE)
#' place.text("$x$",1.1,-0.04,cex=1.3)
#' place.text("$\\int_{-\\infty}^x \\frac{1}{\\sqrt{2\\pi}} e^{-\\frac{x^2}{2}}\\,dx$",0.5,0.3,cex=1.4)
#' @export
place.text=function(text,x,y,col="black",zoom=1,plot=lattice::trellis.last.object(),tex=TRUE, ...){
  dots=list(...)
  if (tex){
    text=latex2exp::TeX(text)
  }
  cex=zoom
  return(mosaic::ladd(lattice::panel.text(x=x,y=y,labels=text,alpha=1,cex=cex,col=col),
                      data=list(x=x,y=y,text=text,col=col,cex=cex),plot=plot))
}

#' Turns on a grid on your plot.
#' @param col grid color
#' @param lty line type, 1,2,or 3 are good options.
#' @param lwd line width, 2 is about right.
#' @param h number of horizontal grid lines. The default value of -1 means to make a line at every tickmark.
#' -2 means to make a line at every other tickmark.  If h is a vector, then make a horizontal gridline at every y-value
#' in the vector.
#' @param v Same as for h, but for vertical lines.
#' @param plot the figure to which to add the grid.
#' @examples
#' #Make a grid with a gridline at every tickmark.
#' plotFunFill(x^2~x)
#' grid.on()
#'
#'#Make 19 evenly spaced horizontal grid lines, and 9 vertical ones.
#'plotFunFill(x^2~x,xlim=c(-2,2))
#'grid.on(h=19,v=9)
#'
#'#Specify where you want the grid lines.
#'plotFunFill(x^2~x,xlim=c(-2,2))
#'grid.on(h=c(0.5,1,1.5,2,2.5,3,3.5,4),v=c(-1.5,-1,-0.5,0,0.5,1,1.5))
#'
#'#Do the example above, but more concisely using the seq command.
#'plotFunFill(x^2~x,xlim=c(-2,2))
#'grid.on(h=seq(0,4,by=0.5),v=seq(-2,2,by=0.5))
#' @export
grid.on=function(h=-1,v=-1,lwd=2,col="black",lty=2,plot=lattice::trellis.last.object()){
  if (length(h)>1 | length(v)>1){
   return(mosaic::ladd(lattice::panel.abline(h=h,v=v,lty=lty,col=col),data=list(h=h,v=v,lty=lty,col=col),plot=plot))
  }
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

#' Places a vector on an already existing figure.  WILL always add to existing figure!
#' @param offset The vector to plot
#' @param base The base of the vector
#' @param col color of the vector
#' @param lwd width of vector
#' @param plot the plot to add to
#' @param ... additional arguments to panel.arrows.
#' @examples
#' blank.canvas(xlim=c(0,5),ylim=c(0,5))
#' grid.on()
#' place.vector(c(1,2),base=c(1,1))
#'
#' @export
place.vector=function(offset,base=c(0,0),col="black",lwd=2,plot=lattice::trellis.last.object(),...){

  #Get the length of the axes.
  xlim=plot$x.limits;
  ylim=plot$y.limits;

  x_len=xlim[[2]]-xlim[[1]];
  y_len=ylim[[2]]-ylim[[1]];

  #Get the size of the figure, in mm.
  plot_size_mm= lattice::current.panel.limits("mm")

  x_len_mm=plot_size_mm$x[[2]]-plot_size_mm$x[[1]]
  y_len_mm=plot_size_mm$y[[2]]-plot_size_mm$y[[1]]

  vec_len_mm=sqrt((offset[[1]]/x_len*x_len_mm)^2+(offset[[2]]/y_len*y_len_mm)^2)

  return(mosaic::ladd(lattice::panel.arrows(x0 = base[[1]], y0 = base[[2]],
                                            x1 = base[[1]]+offset[[1]], y1 = base[[2]]+offset[[2]],
                                            length=grid::unit(0.15*vec_len_mm,"mm"),col = col, lwd = lwd, dots),
                      data=list(base=base,offset=offset,col=col,lwd=lwd,vec_len_mm=vec_len_mm,dots=list(...)),plot=plot))

}

#' Create a blank canvas with convenient limits for later use.
#' @param xlim the x limits
#' @param ylim the y limits
#' @param xlab the label to put on the x-axis
#' @param ylab the label to put on the y-axis.
#' @param ... additional arguments to pass to lattice::xyplot of particular interest might be asp=1 for setting the aspect ratio of the plot.
#' blank.canvas(xlim=c(0,1),ylim=c(0,1),xlab="This is an axis", ylab="This too")
#' blank.canvas(xlim=c(0,1),ylim=c(0,1),xlab="This is an axis", ylab="This too",asp=1)
#' @export
blank.canvas=function(xlim,ylim,xlab="x",ylab="y", ...){

  xs=c(xlim[[1]],xlim[[1]],xlim[[2]],xlim[[2]]);
  ys=c(ylim[[1]],ylim[[2]],ylim[[1]],ylim[[2]])
  return(lattice::xyplot(y~x,data=data.frame(x=xs,y=ys),xlim=xlim,ylim=ylim,col="white",xlab=xlab,ylab=ylab,...))

}

