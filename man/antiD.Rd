% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antiDv1.R
\name{antiD}
\alias{antiD}
\title{A replacement for mosaic::antiD that uses sympy.
WARNING: sympy is very powerful and might return the integral
in terms of functions that R doesn't know about, like it might
return things in terms of Ci, the cosine integral.  In that case,
antiD will succeed, but you won't be able to evaluate the resulting function.}
\usage{
antiD(expression, ...)
}
\arguments{
\item{expression}{the expression to integrate.}

\item{...}{additional arguments defining parameters that might exist in expression.}
}
\description{
A replacement for mosaic::antiD that uses sympy.
WARNING: sympy is very powerful and might return the integral
in terms of functions that R doesn't know about, like it might
return things in terms of Ci, the cosine integral.  In that case,
antiD will succeed, but you won't be able to evaluate the resulting function.
}
\examples{
antiD(sin(x)*x~x)
antiD(exp(-z^2)~z)
}
