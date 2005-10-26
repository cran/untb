\name{untb}
\alias{untb}
\alias{select}
\title{Ecological drift simulation under the Unified Neutral
  Theory of Biodiversity}
\description{
Simulates ecological drift under the UNTB.  Function \code{untb()}
carries out the simulation; function  \code{select()}
carries out a single generational step.
}
\usage{
untb(start, prob.of.mutate=0, D=1, gens=150, keep=FALSE, meta=NULL)
select(a, D=length(a), prob.of.mutate=0, meta=NULL)
}
\arguments{
  \item{a, start}{Starting ecosystem; coerced to class census.  Usually,
    pass an object of class count; see examples.  To start
    with a monoculture of size 10, use \code{start=rep(1,10)} and to
    use \code{start=1:10}.}
  \item{prob.of.mutate}{Probability of mutation}
  \item{D}{Number of organisms that die in each timestep}
  \item{gens}{Number of generations to simulate}
  \item{keep}{In function \code{untb()} Boolean with default
    \code{FALSE} meaning to return the system at the end of the
    simulation and \code{TRUE} meaning to return a matrix whose rows are
    the ecosystem at successive times.}
  \item{meta}{Metacommunity, as represented by a vector of integers
    (representing species; thus being a census).  Can be a
    \dQuote{count} or \dQuote{census} object, in which
    case it is coerced to a census and thence to an integer vector.
    Default of \code{NULL} means to use a \dQuote{greedy} system in
    which every mutation gives rise to a new, previously unencountered
    species.  This would correspond to an infinitely large, infinitely
    diverse, Hubbellian ecosystem (which is not too ridiculous an
    assumption for a small island near a large diverse mainland)}
}
\references{
  S. P. Hubbell.  The Unified Neutral Theory of Biodiversity.
  Princeton University Press, 2001.
}
\author{Robin K. S. Hankin}
\examples{

data(butterflies)
untb(start=butterflies, prob=0, gens=100)

a <- untb(start=1:10,prob=0.005, gens=1000,keep=TRUE)
plot(species.count(a),type="b")
matplot(species.table(a),type="l",lty=1)

}
\keyword{math}
