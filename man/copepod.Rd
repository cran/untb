\name{copepod}
\alias{copepod}
\alias{ostracod}
\docType{data}
\title{Copepod data supplied by Phil Pugh}
\description{
A dataset of copepod abundances supplied by Dr Phil Pugh of the National
Oceanography Centre, Soutampton
}
\usage{data(copepod)}
\format{
A table with names of different copepod species, and entries
corresponding to the numbers of individuals of each species.
}
\source{
\url{http://www.soc.soton.ac.uk/}
}
\examples{
data(copepod)
optimize(f=theta.likelihood,interval=c(10,100), maximum=TRUE,
S=nrow(copepod), J=sum(copepod$count), give.log=TRUE)

data(ostracod)
preston(ostracod)

}
\keyword{datasets}
