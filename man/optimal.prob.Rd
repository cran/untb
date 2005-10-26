\name{optimal.prob}
\alias{optimal.prob}
\title{Returns estimate for mutation probability}
\description{
Returns a maximum likelihood estimate for the probability of mutation
}
\usage{
optimal.prob(x, interval=NULL, ...)
}
\arguments{
  \item{x}{Ecosystem vector or species count table}
  \item{interval}{Interval for probability with default of \code{NULL}
    meaning to use a wide interval}
  \item{\dots}{Further arguments passed to \code{optimize()}}
}
\author{Robin K. S. Hankin}
\examples{
data(butterflies)
optimal.prob(butterflies)
}
\keyword{math}
