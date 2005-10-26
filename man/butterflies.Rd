\name{butterflies}
\alias{butterflies}
\alias{butterfly}
\docType{data}
\title{abundance data for butterflies}
\description{
A dataset giving 
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

# First, estimate the probability of mutation:
nu <-
  optimize(f=theta.likelihood,interval=c(1,100),maximum=TRUE,give.log=TRUE,x=butterflies)$maximum
prob <- nu/(2*sum(butterflies))

#Now plot the dominance diversity curve:
species.curve(butterflies,pch=16,type="o",cex=3,col="red",lwd=3)

#Now plot simulated dominance diversity curves using the above estimation:
for(i in 1:10){
  jj <- untb(n=376,gens=1000,start.mono=TRUE, prob.of.mutate=prob)

  jj <- sort(table(jj),decreasing=TRUE)
  points(1:length(jj),jitter(jj),type="l")
}
}
\keyword{datasets}
