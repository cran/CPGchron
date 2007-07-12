\name{CPGmenu}
\alias{CPGmenu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Main menu for CPGchron }
\description{
  A walk-through menu for accessing and running the entire CPGchron package. 
}
\usage{
CPGmenu()
}
\details{
The running of the CPGchron package essentially occurs in 4 stages:
 
1. Load data

2. MCMC model run

3. Prediction

4. Plotting


The package also includes an option to calibrate radiocarbon dates ignoring any relationship between them
but this is not a necessity for a run of the model, it is simply available for producing better graphics.

In particular loading data needs to occur at the start of every session so option 1 in CPG.menu should 
always be chosen first. Other instructions and a walkthrough is in \code{\link{CPG}}.

}
\value{
  No value is produced, though plots and file output is created through the use of the menus.
}
\references{Haslett and Parnell (2007; submitted to JRSSC) }
\author{Andrew Parnell }

\examples{
    ## CPGmenu()    
}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
