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
\note{
  The fundamental biodiversity parameter \eqn{\theta}{theta} is
  \eqn{2\nu N}{2*nu*N}, where \eqn{\nu}{nu} is the probability of
  mutation (eg, as estimated by \code{optimal.prob()}), and \eqn{N} is
  the size of the ecosystem.
}
\author{Robin K. S. Hankin}
\examples{
data(butterflies)
optimal.prob(butterflies)
}
\keyword{math}
