\name{logkda}
\alias{logkda}
\alias{logkda.a11}
\title{Etienne's K(D,A)}
\description{
Calculates Etienne's \eqn{K(D,A)}
}
\usage{
logkda(a)
logkda.a11(a)
}
\arguments{
  \item{a}{Count object}
}
\details{
Function \code{logkda()} calculates \eqn{K(D,A)} using the method
appearing in Etienne (2005), supplementary online material.  Function
\code{logkda.a11} is a direct implementation of formula A11 in Etienne
(2005).
The formula is
\deqn{
  K(D,A)=
  \sum_{\left\{a_1,\ldots,a_S|\sum a_i=A\right\}}
  \prod_{i=1}^S\frac{
    \overline{s}\left(n_i, a_i\right)
    \overline{s}\left(a_i,   1\right) }{
  \overline{s}\left(n_i,1\right)}}{ommitted...see PDF}
   where \eqn{\overline{s}\left(n_i,a_i\right)} are Stirling numbers of
   the first kind (see \code{logS1}).
 }
\references{
R. S. Etienne 2005. \dQuote{A New Sampling Formula for Neutral
Biodiversity}.  Ecology Letters, volume {8}, pp253--260.
\code{doi: 10.111/j.1461-0248.2004.00717.x}
}
\author{Robin K. S. Hankin; \code{logkda()} is an \R\ transliteration of
  \code{pari/gp} code appearing in Etienne 2005 (supplementary online
  material) due to Chave.
  }
\note{There are two functions here because \code{logkda()}, although
  faster, is based on a non-obvious method.  Function \code{logkda.a11()}
  is much slower (especially for large datasets) but the computational
  method is a direct transcription of formula A11.

  The two functions should produce numerically near-identical output.
  
}
\seealso{\code{\link{etienne}},\code{\link{logS1}}}
\examples{
a <- count(c(dogs=7,pigs=3,crabs=1,hogs=1,slugs=1))
logkda(a)
logkda.a11(a)
}
\keyword{math}
