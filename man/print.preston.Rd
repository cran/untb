\name{print.preston}
\alias{print.preston}
\title{Print and plot objects of class Preston}
\description{
Print and plot objects of class Preston}
}
\usage{
\method{print}{preston}(x, ...)
}
\arguments{
  \item{x}{Object of class \dQuote{preston}}
  \item{\dots}{further arguments passed to \code{print()} after class
    reset} 
}
\author{Robin K. S. Hankin}
\note{
Intended to work with the output of function \code{preston()}.
}
\seealso{\code{\link{preston}}}
\examples{
data(butterflies)
print(preston(butterflies))
}
\keyword{math}
