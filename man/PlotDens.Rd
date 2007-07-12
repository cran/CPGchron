\name{PlotDens}
\alias{PlotDens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot age density at various depths }
\description{
  Computes plots of the density for ages at various (specified) depths. Best used via \code{CPGmenu()} which
  fills in all arguments automatically.
}
\usage{
PlotDens(datesfile, depthfile, fullname,outfile,vers)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datesfile}{ A character string detailing the location and name of the output dates }
  \item{depthfile}{ A character string detailing the location and name of the desired depths }
  \item{fullname}{ The full name of the core }
  \item{outfile}{ A character string detailing the desired location and name of the outputfile }
  \item{vers}{ CPGchron version number used in plots}
}

\value{
  A density plot with highest posterior density ranges (HDR), as well as a file (specified in 'outfile') which
  contains HDR regions.
}
\references{Haslett and Parnell (2007; submitted to JRSSC) }
\author{Andrew Parnell }
\seealso{
\code{\link{CPGmenu}}
\code{\link{hdr2.den}}
}
\examples{
    ## Given in CPGmenu    
}
\keyword{dplot}
\keyword{programming}