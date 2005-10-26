\name{abundance.curve}
\alias{abundance.curve}
\title{Abundance curves}
\description{
Given an ecosystem, plots the ranked abundance curve
}
\usage{
abundance.curve(x, show.uncertainty = FALSE, n = 10, ...)
}
\arguments{
  \item{x}{Ecosystem object, coerced to class count}
  \item{show.uncertainty}{Boolean,
    with \code{TRUE} meaning to show bootstrapped estimates for the
    species diversity curve, and default \code{FALSE} meaning to omit
    this.}
  \item{n}{Number of bootstrapped estimates to plot}
  \item{\dots}{Extra parameters passed to \code{untb()}.}
}
\details{
Plots a ranked abundance curve, optionally with parametrically resampled
datasets showing the uncertainties
  
}
\author{Robin K. S. Hankin}
\examples{
data(copepod)
abundance.curve(copepod)

data(butterflies)
abundance.curve(butterflies,show.uncertainty=TRUE)
}
\keyword{math}
