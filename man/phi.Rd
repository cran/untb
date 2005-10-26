\name{phi}
\alias{phi}
\title{Hubbell's phi}
\description{
Hubbell's phi: counts of species abundances 
}
\usage{
phi(x,addnames=TRUE)
}
\arguments{
  \item{x}{Ecosystem vector; is coerced to class \code{count}}
  \item{addnames}{Boolean with default
    \code{TRUE} meaning to set the name of the \eqn{i}th element
    to the species with abundance \eqn{i} if unique.  Set to
    \code{FALSE} to suppress this, which is useful if the species names
    are long}
}
\details{
  Function \code{phi()} coerces its argument to a \code{count} object and
  by default returns a named vector whose \eqn{i}th element is the
  number of species with \eqn{i} individuals.  The name of the
  \eqn{i}th element is the species with abundance \eqn{i} if unique
  and empty otherwise.  Function \code{phi()} is used by
  \code{theta.prob()}.
}
\references{Hubbell}
\author{Robin K. S. Hankin}
\note{The code for setting the names is a dog's breakfast}
\seealso{\code{\link{preston}}}
\examples{
data(butterflies)
phi(butterflies,add=FALSE)
}
\keyword{math}
