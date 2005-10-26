\name{count}
\alias{count}
\alias{as.count}
\alias{is.count}
\title{Construct, coerce, and test for a count object}
\description{
  In package untb, ecosystem data is held in one of two preferred forms:
  census data and count data.  Function \code{count} creates an object
  of class \dQuote{count}, and \code{as.count()} coerces to this class.
}
\usage{
as.count(a,add="")
count(a)
is.count(a)
}
\arguments{
  \item{a}{Ecosystem data.  In function \code{as.count()}, if a table,
    interpret as species count data; otherwise, interpret as census
    data.  Special dispensation is made for single rows of a dataframe}
  \item{add}{Character, with default \code{""} (empty string) meaning to
    leave the species names unchanged.  A non-empty string is prepended
    to the species names.  This is useful if the species names are
    integers because the display can become confusing}
}
\details{
  Semantically, a \dQuote{count} is a list of species together with
  their abundance.  In~\R\ idiom, a count  is a table of species
  abundances.
  
  An object of class count is a table sorted from most to least abundant
  species.  The singletons are thus tabulated last.

  Function \code{count()} takes a named vector, the elements of which
  are interpreted as abundances and the names of which are interpreted
  as species names.
}
\value{
  Returns an object of class \dQuote{count}.
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{census}}}
\examples{ 
a <- c(rep("oak",5) ,rep("ash",2),rep("elm",3),rep("xx",4))
as.count(a)

data(saunders)
as.count(saunders[1,-(1:150)])

jj <- sample(1:5,5,replace=TRUE)
as.count(jj)
as.count(jj,add="spp.")

}
\keyword{math}
