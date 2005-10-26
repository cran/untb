\name{species.count}
\alias{species.count}
\alias{species.table}
\alias{species.abundance}
\alias{species.curve}
\alias{normalize}
\alias{phi}
\title{Ecosystem diagnostics}
\description{
Provides ecosystem diagnostics of species count datasets (species
counts and species tables), mainly  for the output of \code{untb()}
}
\usage{
species.count(x)
species.table(x)
species.abundance(x)
normalize(x,keepnames=FALSE)
species.curve(x, show.uncertainty = FALSE, n=10, ...)
phi(x)
}
\arguments{
  \item{x}{In functions \code{species.count()} and
    \code{species.table()}, an integer matrix whose rows are integers
    representing the individuals' species.  In function
    \code{species.abundance()} and \code{species.curve()},
    \code{x} is assumed to be an integer
    vector (ie one row of the previous functions' argument) but
    \code{species.curve()} will take a table as well.}
  \item{keepnames}{In function \code{normalize}, Boolean with default
    \code{FALSE} meaning to delete the names and \code{TRUE} meaning to
    retain them.  Retaining names can be confusing.}
  \item{show.uncertainty}{In function \code{species.curve()}, Boolean
    with \code{TRUE} meaning to show bootstrapped estimates for the
    species diversity curve, and default \code{FALSE} meaning to omit
    this.}
  \item{n}{In function \code{species.curve()}, the number of
    bootstrapped estimates to plot}
  \item{...}{In function \code{species.curve()}, extra parameters passed
    to \code{plot()}.}
}
\details{
  Function \code{species.count()} returns the total number of species
  present as a function of time

  Function \code{species.table()} returns a matrix \eqn{M} where
  \code{M[i,j]}   column of the matrix is the abundance of species \eqn{j} 
  at time \code{i}.

  Function \code{species.abundance()} takes an integer vector and
  returns a table of abundances of species, ordered by abundance.
  (Actually it's not a table but a two-row matrix because a matrix can
  have labelled rows and a table can't).
 
  Function \code{normalize()} takes an ecosystem vector and renames the
  species so that the most abundant one is \dQuote{1}, the second most
  abundant is \dQuote{2}, and so on.

  Function \code{phi()} takes an ecosystem vector and returns a vector
  whose \eqn{i}th element is the number of species with \eqn{i}
  individuals; this is used by function \code{theta.prob()}.

  The  functions work well with the output of function \code{untb()}.

 Function \code{species.curve()} plots a dominance diversity curve.
}
\value{
  Returns a vector
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{preston}}}
\examples{
a <- untb(n=50, prob=0.01, gens=2000, keepall=TRUE, start.mono=TRUE)

plot(species.count(a),type="b")
matplot(species.table(a),type="l",lty=1)

jj <- a[2000,]
print(jj)
species.abundance(jj)

gg <- c(rep("a",10),rep("b",3),letters[5:9])
normalize(gg)
phi(gg)

## Thus there are 5 species with one individual (ie singletons) , zero
## species with two individuals, one species with three individuals, and
## so on.


data(butterflies)
phi(butterflies)

}
\keyword{math}
