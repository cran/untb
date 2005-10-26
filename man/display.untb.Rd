\name{display.untb}
\alias{display.untb}
\title{Animation of neutral ecological drift}
\description{
Displays an ongoing simulation of neutral ecological drift using nice
colours and a simple animation technique
}
\usage{
display.untb(start, gens=100, prob.of.mutate = 0.001, cex=3, individually
= TRUE, flash = FALSE, t1 = 0, flashsleep = 0.1, ...)
}
\arguments{
  \item{start}{Starting ecosystem; coerced to class census.  Usually,
    pass an object of class count; see examples.  To start
    with a monoculture of size 10, use \code{start=rep(1,10)} and to
    use \code{start=1:10}.}
  \item{prob.of.mutate}{Probability of mutation.}
  \item{individually}{Boolean, indicating whether the ecosystem should
    evolve simultaneously (\code{FALSE}), or individually (\code{TRUE})}
  \item{gens}{Number of generations to simulate}
  \item{cex}{The size of the dots
    used for plotting, defaulting to 3}
  \item{flash}{Boolean, with \code{TRUE} meaning  to indicate the site
    of a death/birth with a flashing ring; and default \code{FALSE}
    meaning to omit the flashing ring and thus simulate \emph{much} faster}
  \item{t1}{Time delay between generations}
  \item{flashsleep}{In function \code{display.untb()}, time delay of
    flashing ring}
  \item{...}{Further arguments passed
    to \code{plot()} and \code{points()}.}
}
\references{
  S. P. Hubbell.  The Unified Neutral Theory of Biodiversity.
  Princeton University Press, 2001.
}
\author{Robin K. S. Hankin}
\examples{

data(butterflies)
display.untb(start=butterflies,prob=0, gens=1e2)
}
\keyword{math}
