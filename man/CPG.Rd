\name{CPG}
\alias{CPG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ An MCMC run of the CPG model }
\description{
  This function will perform a run of the MCMC for the compound Poisson Gamma (CPG) model for radiocarbon dated
  depth chronologies. It is best accessed through the \code{\link{CPGmenu}}. 
}
\usage{
CPG(CALPATH, INFILE, OUTFILE, ndet, iterations, burnin, howmany, thinby)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CALPATH}{ A character string detailing the location and name of the calibration curve file.
  Note that CPGchron uses a specific version of the Northern hemisphere Intcal calibration curve. Other calibration
  curves are not supported. }
  \item{INFILE}{ A character string detailing the location and name of the determinations, as required by
  CPGchron package. See \code{\link{CPGmenu}} for more details. }
  \item{OUTFILE}{ A character string detailing the desired location and name of the output file. }
  \item{ndet}{ The number of radiocarbon determinations for this core. }
  \item{iterations}{ The number of iterations required for the MCMC. }
  \item{burnin}{ The length of the MCMC burn-in. }
  \item{howmany}{ How often you wish the program to report updates in iterations. }
  \item{thinby}{ How much thinning should occur to the MCMC output. }
}
\value{
  No output produced but the file specified in \code{OUTFILE} is populated by parameter values.
}
\references{Haslett and Parnell (2007; submitted to JRSSC) }
\author{Andrew Parnell }
\examples{
    ## Given in CPGmenu    
}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}