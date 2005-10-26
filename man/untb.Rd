\name{untb}
\alias{untb}
\alias{display.untb}
\alias{select}
\title{Ecological drift simulation under the Unified Neutral
  Theory of Biodiversity}
\description{
Simulates ecological drift under the UNTB.  Function \code{untb()}
carries out the simulation; function \code{display.untb()} displays
an ongoing simulation using nice colours.  Function \code{select()}
carries out a single generational step.
}
\usage{
untb(n = 100, prob.of.mutate = 0.01,  individually = FALSE, gens = 150, keepall = FALSE,
start.mono = FALSE)
display.untb(n = 400, gens = 9e+05, prob.of.mutate = 0.001, cex=3, individually
= TRUE, start.mono = TRUE, flash = FALSE, t1 = 0, flashsleep = 0.1, ...)
select(a, prob.of.mutate = 0)
}
\arguments{
  \item{n}{Size of ecosystem}
  \item{prob.of.mutate}{Probability of mutation}
  \item{individually}{Boolean, indicating whether the ecosystem should
    evolve simultaneously (\code{FALSE}), or individually (\code{TRUE})}
  \item{gens}{Number of generations to simulate}
  \item{cex}{In function \code{display.untb()}, the size of the dots
    used for plotting, defaulting to 3}
  \item{keepall}{In function \code{untb()} Boolean with default
    \code{FALSE} meaning to return the system at the end of the
    simulation and \code{TRUE} meaning to return a matrix whose rows are
    the ecosystem at successive times.}
  \item{start.mono}{Boolean, with default \code{FALSE} meaning to start
    the simulation with a system of maximal diversity and \code{TRUE} to
    start with a monoculture}
  \item{flash}{In function \code{display.untb()},
    Boolean argument with  \code{TRUE} meaning  to indicate the site
    of a death/birth with a flashing ring; and default \code{FALSE}
    meaning to omit the flashing ring and thus simulate \emph{much} faster}
  \item{t1}{In function \code{display.untb()}, time delay between generations}
  \item{flashsleep}{In function \code{display.untb()}, time delay of
    flashing ring}
  \item{a}{Ecosystem start state for use in function \code{select()}}
  \item{...}{In function \code{display.untb()}, further arguments passed
    to \code{plot()} and \code{points()}.}
}
\value{
  Function \code{untb()} returns the final state of the ecosystem.
}
\references{
  S. P. Hubbell.  The Unified Neutral Theory of Biodiversity.
  Princeton University Press, 2001.
}
\author{Robin K. S. Hankin}
\examples{
untb(n=10, prob=0.1, gens=100)

a <- untb(n=10,prob=0.05, gens=1000,keep=TRUE)
plot(species.count(a),type="b")
matplot(species.table(a),type="l",lty=1)

display.untb(prob=0.1, gens=1e2)
}
\keyword{math}
