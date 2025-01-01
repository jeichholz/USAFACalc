#' Animate the numerical solution to an ODE. The expectation is that you have a numerical solution to an ODE coming from Euler, and you'd like to animate some plot.  The time variable *must* be in a column called t.
#' @param relation1 The first relationship you want to plot.  For instance yearth~xearth
#' @param ...  more relationships.  Like ysun~xsun or ymoon~xmoon
#' @param data the data frame containing the numerical solution to the ODE.
#' @param duration the duration, in seconds, that the gif will last.
#' @param fps the number of frames per second of duration to produce.
#' @examples
#' #Simulate the earth orbiting the sun, make an animation.
#'
#' #Mass of sun in KG
#' sun_mass=1.9891e30
#'
#' #Mass of earth in KG
#' earth_mass= 5.97219e24
#'
#' #distance from sun to earth in m.
#' earth_sun_distance=149597870700
#'
#' #velocity of earth relative to sun in m/s
#' earth_vel=29784.8
#'
#' #G in kg m s
#' G=6.6743e-11
#'
#' #distance between two objects
#' dist=mosaic::makeFun(sqrt((x2-x1)^2+(y2-y1)^2)~x1&y1&x2&y2)

#' #Calculate the force between two objects. Directed as toward object 2.
#' Fx=mosaic::makeFun(G*(x2-x1)*m1*m2/dist(x1,y1,x2,y2)^3~m1&x1&y1&m2&x2&y2)
#' Fy=mosaic::makeFun(G*(y2-y1)*m1*m2/dist(x1,y1,x2,y2)^3~m1&x1&y1&m2&x2&y2)

#' #The DE goes in the order
#' #earthx'
#' #earthy'
#' #earthvx'
#' #earthvy'
#' #sunx'
#' #suny'
#' #sunvx'
#' #sunvy'
#'
#' #Define RHS of the ODE.
#' ode=c(earthvx,
#'       earthvy,
#'       1/earth_mass*Fx(earth_mass,earthx,earthy,sun_mass,sunx,suny),
#'       1/earth_mass*Fy(earth_mass,earthx,earthy,sun_mass,sunx,suny),
#'       sunvx,
#'       sunvy,
#'       1/sun_mass*Fx(sun_mass,sunx,suny,earth_mass,earthx,earthy),
#'       1/sun_mass*Fy(sun_mass,sunx,suny,earth_mass,earthx,earthy))~earthx&earthy&earthvx&earthvy&sunx&suny&sunvx&sunvy
#'
#' #There are 3.154e7 seconds in a year.
#' endT=2*3.154e7;
#'
#' #There are 86400 seconds in a day.
#' dt=86400*1/12
#' #simulate!
#' soln=Euler(ode,tlim=c(0,endT), stepSize=dt,ic=c(earthx=earth_sun_distance,earthy=0,
#'                                                 earthvx=0,earthvy=earth_vel,
#'                                                 sunx=0,suny=0,
#'                                                 sunvx=0,sunvy=0))
#'
#' animateODESolution(earthy~earthx,suny~sunx,data=soln)
#' @export
animateODESolution=function(relation1,...,data,duration=5,fps=20){
  dots=list(...)
  l=as.character(mosaic::lhs(relation1))
  r=as.character(mosaic::rhs(relation1))
  toplot=data[,c("t",l,r)]
  toplot=cbind(toplot,rep(deparse(relation1)))
  names(toplot)=c("t","y","x","group")
  ylabs=l
  xlabs=r

  for (rel in dots){
    l=as.character(mosaic::lhs(rel))
    r=as.character(mosaic::rhs(rel))
    tmp=data[,c("t",l,r)]
    tmp=cbind(tmp,rep(deparse(rel)))
    names(tmp)=c("t","y","x","group")
    toplot=rbind(toplot,tmp)
    xlabs=paste0(xlabs,"/",r)
    ylabs=paste0(ylabs,"/",l)
  }

  p<-ggplot2::ggplot(data=toplot,mapping=ggplot2::aes(x=x,y=y,col=group))+ggplot2::geom_point()+
    ggplot2::labs(title="Time: {frame_time}")+
    ggplot2::xlab(xlabs)+ggplot2::ylab(ylabs)+
    ggplot2::theme_bw()+gganimate::transition_time(t)

  gganimate::animate(p,duration=duration,fps=fps)

}
