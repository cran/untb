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
  A dataframe with each row representing species inventories for a kelp
  holdfast (\code{saunders}) together with a Boolean flag indicating
  whether the holdfast was in a sheltered or exposed location.  Also two
  data frames, one for the 20 exposed holdfasts (\code{saunders.exposed})
  and one for the 20 sheltered holdfasts (\code{saunders.sheltered}).

  Also three \code{count} objects, giving counts for all organisms
  (\code{saunders.tot}), all those from exposed locations
  (\code{saunders.exposed.tot}), and all those from sheltered locations
  only (\code{saunders.sheltered.tot}).
}
\usage{data(saunders)}
\format{
  Dataset \code{saunders} is a data frame with 40 observations on 177
  variables.  The observations (rows) are 40 holdfasts and the variables
  are, with the exception of the first one, the numbers of each of 176
  species living there.  The first variable is a Boolean, indicating
  whether or not that holdfast was exposed (\code{TRUE}) or sheltered
  (\code{FALSE}).
  
  Summary datasets \code{saunders}, \code{saunders.exposed}, and
  \code{saunders.sheltered} are objects of class \code{count} that
  are the species abundance for sheltered holdfasts, exposed
  holdfasts, and the entire dataset.
}
\details{

  Kelp are large seaweeds classified in kingdom Chromista.  Kelp grows
  in shallow oceans in kelp forests.
    
  The \dfn{holdfast} is a root-like structure that anchors the kelp to
  the ocean floor.  Fauna inhabiting kelp holdfasts, being
  \dQuote{incredibly diverse} (Anderson et al 2005), are often used as
  indicators of environmental change.

  The data was collected in New Zealand, from eight sites along the
  Leigh coastline from north of Leigh Harbour down to the southern end
  of Kawau Island (a stretch of roughly 20 km).  Four sites were
  wave-exposed, four were sheltered (although two of the latter were
  arguably quite tidally-dominated).  Each site had a spatial extent of
  roughly one hectare.  They were collected from 5 - 10 November, 2003
  
}  
\note{
  It is not entirely obvious how to derive the summary datasets from the
  \code{saunders} dataframe.  Use function \code{extractor()} for this.
}
\source{
  Data supplied by Justine Saunders
}
\seealso{\code{\link{extractor}}}
\references{
  \itemize{
    \item J. Saunders 2007.  \dQuote{Biodiversity of kelp holdfasts}
      (provisional title).  PhD thesis (in preparation); School of Geography
      and Environmental Sciences, The University of Auckland
      \item M. J. Anderson and others 2005.  \dQuote{Consistency and
	variation in kelp holdfast assemblages: Spatial patterns of
	biodiversity for the major phyla at different taxonomic resolutions}.
      Journal of Experimental Marine Biology and Ecology. volume 320, pages
      35-56
    }
}
\examples{
data(saunders)
plot(saunders.sheltered.tot, uncertainty=TRUE, n=1)

preston(saunders.tot)

plot(extant(extractor(saunders,1)))

}
\keyword{datasets}
