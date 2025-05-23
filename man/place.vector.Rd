% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additionalPlotting.R
\name{place.vector}
\alias{place.vector}
\title{Places a vector on an already existing figure.  WILL always add to existing figure!}
\usage{
place.vector(
  offset,
  base = c(0, 0),
  col = "black",
  lwd = 2,
  plot = lattice::trellis.last.object(),
  ...
)
}
\arguments{
\item{offset}{The vector to plot}

\item{base}{The base of the vector}

\item{col}{color of the vector}

\item{lwd}{width of vector}

\item{plot}{the plot to add to}

\item{...}{additional arguments to panel.arrows.}
}
\description{
Places a vector on an already existing figure.  WILL always add to existing figure!
}
\examples{
blank.canvas(xlim=c(0,5),ylim=c(0,5))
grid.on()
place.vector(c(1,2),base=c(1,1))

}
