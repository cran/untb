\name{census}
\alias{census}
\alias{as.census}
\alias{is.census}
\title{Construct, coerce, and test for a census object}
\description{
In package untb, ecosystem data is held in one of two preferred forms:
census data and count data.  Function \code{as.census()} coerces to
census format 
}
\usage{
census(a)
as.census(a)
is.census(a)
}
\arguments{
  \item{a}{Ecosystem data.  In function \code{as.census()}, if a table,
    interpret as species count data; otherwise, interpret as census data}
}
\details{
  Semantically, a \dQuote{census} is a list of individuals.  In \R
  idiom, a census is a vector whose elements indicate the individuals'
  species.

  An object of class census is an unordered factor.  The levels are
  always in alphabetical order; the vector itself is sorted so that the
  individuals of the most abundant species appear first, individuals
  of the second most abundant species appear second, and so on: the
  singletons appear last.

  Function \code{census()} takes an object of class \dQuote{count} and
  returns an object of class \dQuote{census}
}
\value{
  Returns an object of class \dQuote{census}.
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{count}}}
\examples{

a <- c(rep("oak",5) ,rep("ash",2),rep("elm",3),rep("xx",4))
# note that "a" is a plain vector here.
as.census(a)
}
\keyword{math}
