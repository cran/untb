\name{saunders}
\alias{saunders}
\docType{data}
\title{Dataset due to Saunders}
\description{
A list of 40 holdfasts and a species inventory for each
}
\usage{data(saunders)}
\format{
  A data frame with 40 observations on 177 variables.  The observations
  are 40 holdfasts and the variables are, with the exception of the
  first one, the numbers of each of 176 species living there.  The first
  variable is a Boolean, indicating whether or not that holdfast was
  exposed (\code{TRUE}) or sheltered (\code{FALSE}).
}
\details{
Extract rows for a species table of a particular holdfast
}
\source{
  Data supplied by Justine Saunders
}
\examples{
data(saunders)
preston(as.count(saunders[1,-1]))

abundance.curve(count(apply(saunders,2,sum)[-1]))

}
\keyword{datasets}
