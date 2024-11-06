#' Visualize the left-hand rule
#' @param f an expression for the function
#' @param xlim limits of integration
#' @param N number of rectangles.
#' @examples
#' Reimann.LH(exp(-x^2/2)~x,xlim=c(-5,5),N=10)
#'

#' @export
Reimann.LH=function(f,xlim=c(0,5),N=10){

  deltax=(xlim[[2]]-xlim[[1]])/N;
  A=plotFunFill(f,xlim=xlim,col="chocolate1",alpha=0.9)
  ffun=mosaic::makeFun(f)
  for (i in 1:N){
    xi=xlim[[1]]+(i-1)*deltax
    C=ffun(xi)
    val=mosaic::makeFun(C*(x*0+1)~x,C=C)
    A=plotFunFill(val(x)~x,xlim=c(xi,xi+deltax),col="darkmagenta",add=TRUE,addoutline=FALSE,plot=A)
  }
  return(A)

}
#' Visualize the right-hand rule
#' @param f an expression for the function
#' @param xlim limits of integration
#' @param N number of rectangles.
#' @examples
#' Reimann.RH(exp(-x^2/2)~x,xlim=c(-5,5),N=10)
#' @export
Reimann.RH=function(f,xlim=c(0,5),N=10){

  deltax=(xlim[[2]]-xlim[[1]])/N;
  A=plotFunFill(f,xlim=xlim,col="chocolate1",alpha=0.9)
  ffun=mosaic::makeFun(f)
  for (i in 1:N){
    xi=xlim[[1]]+(i-1)*deltax
    C=ffun(xi+deltax)
    val=mosaic::makeFun(C*(x*0+1)~x,C=C)
    A=plotFunFill(val(x)~x,xlim=c(xi,xi+deltax),col="darkmagenta",add=TRUE,addoutline=FALSE,plot=A)
  }
  return(A)

}

#' Visualize the midpoint rule
#' @param f an expression for the function
#' @param xlim limits of integration
#' @param N number of rectangles.
#' @examples
#' Reimann.MP(exp(-x^2/2)~x,xlim=c(-5,5),N=10)
#' @export
Reimann.MP=function(f,xlim=xlim,N=10){

  deltax=(xlim[[2]]-xlim[[1]])/N;

  fig=plotFunFill(f,xlim=xlim,col="chocolate1",alpha=0.9)
  ffun=mosaic::makeFun(f)
  for (i in 1:N){
    xi=xlim[[1]]+(i-1)*deltax
    A=ffun(xi+deltax/2)
    val=mosaic::makeFun(A*(x*0+1)~x,A=A)
    fig=plotFunFill(val(x)~x,xlim=c(xi,xi+deltax),col="darkmagenta",add=TRUE,addoutline=FALSE,plot=fig)
  }

  return(fig)
}


