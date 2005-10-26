\name{no.of.ind}
\alias{no.of.ind}
\alias{no.of.spp}
\alias{no.of.singletons}
\alias{no.of.extinct}
\alias{singletons}
\alias{extinct}
\alias{extant}
\title{Ecosystem diagnostics}
\description{
Ecosystem diagnostics such as species count, individual count, number of
singletons, etc
}
\usage{
no.of.ind(x)
no.of.spp(x, include.extinct=FALSE)
no.of.singletons(x)
singletons(x)
extinct(x)
extant(x)
}
\arguments{
  \item{x}{Ecosystem vector; is coerced to class \code{count}}
  \item{include.extinct}{In function \code{no.of.spp()}, Boolean
    argument with \code{TRUE} meaning to include extinct species (ie
    species with an abundance of zero), and default \code{FALSE}
    meaning to return the number of extant species.}
}
\details{
  Function \code{no.of.spp()} returns the number of species in an
  ecosystem object, treating exctinct species in line with argument
  \code{include.extinct}.

  Function \code{no.of.ind()} returns the number of individuals.

  Function \code{no.of.singletons()} returns the number of singletons.

  Function \code{no.of.extinct()} returns the number of extinct species.

  Function \code{singletons()} returns a \code{count} object containing
  only the singletons: each abundance is one.

  Function \code{extinct()} returns a \code{count} object containing
  only the extinct species: each abundance is zero.

  Function \code{extant()} returns a \code{count} object containing
  only the extant species: each abundance is greater than zero.

}
\references{Hubbell}
\author{Robin K. S. Hankin}
\note{It is sometimes useful to include species with an abundance of
  zero when, for example, taking a single row of the Saunders dataframe.

  The default for \code{include.extinct} is \code{FALSE} because this is
  required for (eg) \code{optimal.theta()}.
  
  In function \code{singletons()}, the code for setting the names
  is a dog's breakfast}
\examples{
data(butterflies)
no.of.spp(butterflies)
no.of.ind(butterflies)
}
\keyword{math}
