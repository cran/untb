\name{simpson}
\alias{simpson}
\title{Simpson's diversity index}
\description{
Simpson's diversity index
}
\usage{
simpson(x)
}
\arguments{
  \item{x}{Ecosystem vector; coerced to class \code{count}}
}
\details{
  Returns the Simpson index \eqn{D}: the
  probability that two randomly sampled individuals belong to
  different species.  This is 
  \deqn{1-\sum_{i=1}^S\left(p_i\right)^2}{1-sum((p_i)^2)}

  where \eqn{p_i} is the relative frequency of the \eqn{i^{\rm th}}{i-th}
  species  
}
\references{
  \itemize{
    \item S. P. Hubbell 2001.
    \dQuote{The Unified Neutral Theory of Biodiversity}.
    Princeton University Press.
  
    \item F. He and X.-S. Hu 2005.
    \dQuote{Hubbell's Fundamental Biodiversity
      Parameter and the Simpson Diversity Index}.  Ecology Letters, volume 8,
    pp386-390. doi: \code{10.1111/j.1461-0248.2005.00729.x}
    \item E. H. Simpson 1949. \dQuote{Measurement of diversity}.
    Nature, volume 163, p688
  }  
}
\author{Robin K. S. Hankin}
\note{
  Early versions of the \pkg{untb} package effectively defined \eqn{D}
  as the probability that two randomly selected \emph{different}
  individuals belong to different species
}
\seealso{\code{\link{preston}}}
\examples{
data(butterflies)

D <- simpson(butterflies)
theta <- optimal.prob(butterflies)*2*no.of.ind(butterflies)

#compare theta with D/(1-D) (should be roughly equal):

theta
D/(1-D)
}
\keyword{math}
