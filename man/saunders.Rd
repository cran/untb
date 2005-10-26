\name{saunders}
\alias{saunders}
\alias{saunders.tot}
\alias{saunders.exposed}
\alias{saunders.exposed.tot}
\alias{saunders.sheltered}
\alias{saunders.sheltered.tot}
\docType{data}
\title{Dataset due to Saunders}
\description{
A dataframe of 40 holdfasts and a species inventory for each
(\code{saunders}); also two \code{data.frame.count} objects:
one for the 20 exposed holdfasts (\code{saunders.exposed}) and one for
the 20 sheltered holdfasts (\code{saunders.sheltered}); and three
\code{count} objects: \code{saunders.tot}, \code{saunders.exposed}, and
\code{saunders.sheltered}.
}
\usage{data(saunders)}
\format{
  Dataset \code{saunders} is a data frame with 40 observations on 177
  variables.  The observations are 40 holdfasts and the variables are,
  with the exception of the first one, the numbers of each of 176
  species living there.  The first variable is a Boolean, indicating
  whether or not that holdfast was exposed (\code{TRUE}) or sheltered
  (\code{FALSE}).

  Summary datasets \code{saunders}, \code{saunders.exposed}, and
  \code{saunders.sheltered} are objects of class \code{count} that
  are the species abundance for sheltered holdfasts, exposed
  holdfasts, and the entire dataset.
}
\note{
  It is not entirely obvious how to derive the summary datasets from the
  \code{saunders} dataframe.  Use function \code{extractor()} for this.
}

\source{
  Data supplied by Justine Saunders
}
\seealso{\code{\link{extractor}}}
\examples{
data(saunders)
plot(saunders.sheltered.tot, uncertainty=TRUE, n=1)

preston(saunders.tot)
}
\keyword{datasets}
