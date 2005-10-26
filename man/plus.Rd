\name{+.count}
\alias{+.count}
\title{Add two count objects}
\description{
Adds two count objects
}
\usage{
\method{+}{count}(a,b)
}
\arguments{
  \item{a,b}{objects of class \code{count}}
}
\author{Robin K. S. Hankin, based on a tip from Gabor Grothendiek}
\examples{
a <- count(c(dogs=4,pigs=0,slugs=5))
b <- count(c(slugs=4,hogs=1,frogs=19))

a+b

}
\keyword{math}
