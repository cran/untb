\name{preston}
\alias{preston}
\alias{Preston}
\title{Preston diagram of an ecosystem}
\description{
  Gives a standard Preston diagram for an ecosystem.
}
\usage{
preston(x,n=8)
}
\arguments{
  \item{x}{Either a table, or an Integer vector whose entries are the
    species identifier of each individual}
  \item{n}{An integer specifying the number of species abundance classes
      to use.  Must be greater than 2.}
}

\details{
  The Preston diagram is a table showing the number of species having
  abundances in specified abundance classes.  Consider the following
  Preston diagram:
\preformatted{
                    (0,1] (1,2] (2,4] (4,8] (8,16] (16, Inf)
number of species    73    33    24    36     39      55
}

This shows that there are 73 species withe abundance 1 (that is,
singletons), 33 species with abundance 2-4, 24 species with abudances
4-8 and so on.

Note that some authors count any species with an abundance on the
boundary between two adjacent abundance classes as being split 50-50
between the classes.  This is not done here; note the type of brackets
used in the diagram above that show whether the interval is open or
closed (here, open on the left and closed on the right, apart from the
most abundant class).

}
\references{Preston}
\author{Robin K. S. Hankin}
\examples{
preston(untb(n=1000, prob=0.01, gens=1000, keep=FALSE))

data(butterflies)
preston(butterflies)


data(copepod)
preston(copepod)
preston(copepod, n=21)

}
\keyword{math}
