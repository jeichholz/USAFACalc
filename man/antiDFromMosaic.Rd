% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antiDv1.R
\name{antiDFromMosaic}
\alias{antiDFromMosaic}
\title{Calculate an antiderivative, error if no symbolic antiderivative  A wrapper around
mosaicCalc::antiD.}
\usage{
antiDFromMosaic(tilde)
}
\arguments{
\item{tilde}{the expression to integrate}
}
\value{
A function that is the symbolic antiderivative of tilde. Error if no
symbolic antiderivative found.
}
\description{
Calculate an antiderivative, error if no symbolic antiderivative  A wrapper around
mosaicCalc::antiD.
}
\examples{
#Succeeds
#antiD(x^2~x)

#Succeeds
#antiD(cos(2*x)~x)

#This will fail:
#antiD(dnorm(x)~x)

#If you really want to use numerical anti differentiation for an antiderivative
#do this instead:
#mosaic::makeFun(integrate(dnorm(s)~s,slim=c(0,x))~x)
}
