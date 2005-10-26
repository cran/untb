\name{count.to.census}
\alias{count.to.census}
\title{Converts count data to census data}
\description{
Given count data in table form, return the corresponding ecosystem.
}
\usage{
count.to.census(a)
}
\arguments{
  \item{a}{Species count data in the form of a dataframe}
}
\author{Robin K. S. Hankin}
\examples{
data(butterflies)
count.to.census(butterflies)
normalize(count.to.census(butterflies))
}
\keyword{math}
