\name{untb-package}
\alias{untb-package}
\docType{package}
\title{
Unified neutral theory of biodiversity
}
\description{
  Numerical simulations, and visualizations, of the unified
  neutral theory of biodiversity
}
\details{
\tabular{ll}{
Package: \tab untb\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2005-10-26\cr
License: \tab GPL\cr
}
The package has two main engines: \code{untb()}, which tracks an
ecosystem through neutral drift; and \code{display.untb()}, which
displays a semi-animated graphic.

Function \code{untb()} produces by default an integer vector of length
\code{n}.  Each element of this vector corresponds to an organism
(member of an ecosystem) and its species is the value of the element.


Functions \code{species.count()} et seq compute ecologically interesting
diagnostics such as species count and so on.
}
\author{
Robin K. S. Hankin

Maintainer: <r.hankin@noc.soton.ac.uk>
}
\references{
  S. P. Hubbell.  The Unified Neutral Theory of Biodiversity.
  Princeton University Press, 2001.
}
\keyword{math}
}
\examples{
a <- untb(1000,prob=0.005,gens=5000,keep=FALSE,start.mono=TRUE)
preston(a)
species.count(a)
species.abundance(a)

display.untb(n=100,prob=0.1,gens=1000)

data(butterflies)
species.curve(butterflies,show=TRUE)

}
