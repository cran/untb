\name{theta.prob}
\alias{theta.prob}
\alias{theta.likelihood}
\title{Posterior probabilities for theta}
\description{
  Determines the posterior probability   (\code{theta.prob()}) and
  likelihood (\code{theta.likelihood()}) for theta, given an ecosystem.
}
\usage{
theta.prob(theta, x, give.log = FALSE)
theta.likelihood(theta,x=NULL,S=species.count(x),J=length(x),give.log=FALSE)
}
\arguments{
  \item{theta}{Biodiversity parameter}
  \item{x}{Ecosystem vector}
  \item{give.log}{Boolean, with default \code{FALSE} meaning to return
the probability or likelihood, and \code{TRUE} meaning to return
the logarithm of the probability, or the log-likelihood}
\item{S}{In function \code{theta.likelihood()}, if \code{x} is not
supplied, the number of species in the ecosystem.}
\item{J}{In function \code{theta.likelihood()}, if \code{x} is not
supplied, the number of individuals in the ecosystem.  Arguments
\code{S} and \code{J} are provided so that \code{x} need not be supplied
if \code{S} and \code{J} are known.}
}
\references{
  S. P. Hubbell.  The Unified Neutral Theory of Biodiversity.
  Princeton University Press, 2001.

\url{http://en.wikipedia.org/wiki/Untb}
}
\author{Robin K. S. Hankin}
\note{
  If estimating \code{theta}, use \code{theta.likelihood()} rather than
  \code{theta.probability()}.
}
\seealso{\code{\link{phi}}, \code{\link{optimal.prob}}}
\examples{

gg <- c(rep("a",10),rep("b",3),letters[5:9])
theta.likelihood(theta=2,gg)

optimize(f=theta.likelihood,interval=c(0,100),maximum=TRUE,x=gg)

a <- untb(n=1000,gens=1000,prob=1e-3,start.mono=TRUE)

\dontrun{
optimize(f=theta.likelihood,interval=c(0,100),maximum=TRUE,x=a)
}
## Fails due to numerical overshoot

optimize(f=theta.likelihood,interval=c(0,100),maximum=TRUE,x=a,give.log=TRUE)
## Compare the true value of 2 (=2*1e-3*1000).  Not bad
}
\keyword{math}