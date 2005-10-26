\name{phi}
\alias{phi}
\alias{no.of.spp}
\alias{no.of.ind}
\alias{no.of.singletons}
\alias{singletons}
\alias{simpson}
\title{Ecosystem diagnostics}
\description{
Ecosystem diagnostics such as Hubbell's phi and species count
}
\usage{
phi(x,addnames=TRUE)
no.of.spp(x)
no.of.ind(x)
no.of.singletons(x)
singletons(x)
simpson(x)
}
\arguments{
  \item{x}{Ecosystem vector; is coerced to class \code{count}}
  \item{addnames}{In function \code{phi()}, Boolean with default
    \code{TRUE} meaning to set the name of the \eqn{i}th element
    to the species with abundance \eqn{i} if unique.  Set to
    \code{FALSE} to suppress this, which is useful if the species names
    are long}
}
\details{
  Function \code{phi()} coerces its argument to a count object and
  by default returns a named vector whose \eqn{i}th element is the
  number of species with \eqn{i} individuals.  The name of the
  \eqn{i}th element is the species with abundance \eqn{i} if unique
  and empty otherwise.  Function \code{phi()} is used by
  \code{theta.prob()}.

  Function \code{no.of.spp()} returns the number of species in an
  ecosystem object.

  Function \code{no.of.ind()} returns the number of individuals.

  Function \code{no.of.singletons()} returns the number of singletons.

  Function \code{singletons()} returns the names of the singletons.

  Function \code{simpson()} returns the Simpson index \eqn{D}: the
  probability that two randomly sampled individuals belonging to
  different species
}
\references{Hubbell}
\author{Robin K. S. Hankin}
\note{The code for setting the names is a dog's breakfast}
\seealso{\code{\link{preston}}}
\examples{
data(butterflies)
phi(butterflies,add=FALSE)

no.of.spp(butterflies)
no.of.ind(butterflies)

D <- simpson(butterflies)
theta <- optimal.prob(butterflies)*2*no.of.ind(butterflies)

#compare theta with D/(1-D) (should be equal):

theta
D/(1-D)




}
\keyword{math}
