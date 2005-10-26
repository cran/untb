\name{simpson}
\alias{simpson}
\title{Simpson's diversity index}
\description{
Simpson's diversity index
}
\usage{
simpson(x)
}
\arguments{
  \item{x}{Ecosystem vector; is coerced to class \code{count}}
}
\details{
  Returns the Simpson index \eqn{D}: the
  probability that two randomly sampled individuals belonging to
  different species
}
\references{Hubbell}
\author{Robin K. S. Hankin}
\note{The code for setting the names is a dog's breakfast}
\seealso{\code{\link{preston}}}
\examples{
data(butterflies)

D <- simpson(butterflies)
theta <- optimal.prob(butterflies)*2*no.of.ind(butterflies)

#compare theta with D/(1-D) (should be equal):

theta
D/(1-D)
}
\keyword{math}
