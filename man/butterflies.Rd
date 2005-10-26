\name{butterflies}
\alias{butterflies}
\alias{butterfly}
\docType{data}
\title{abundance data for butterflies}
\description{
A dataset of class \dQuote{count} showing the abundance of several
butterfly species
}
\usage{data(butterflies)}
\format{
A table with names of different butterfly species, and entries
corresponding to the respective numbers of individuals.
}
\source{
  \url{http://texasnaturalist.net/butterflies/bflyct96.html}
}
\references{
Texas Birding and Naturalist Web
}
\examples{
data(butterflies)
plot(butterflies, uncertainty=TRUE)

}
\keyword{datasets}
