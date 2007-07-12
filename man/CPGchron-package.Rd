\name{CPGchron-package}
\alias{CPGchron-package}
\alias{CPGchron}
\docType{package}
\title{
Compound Poisson-Gamma chronologies
}
\description{
Runs MCMC and predictive distributions for the monotonic piece-wise linear compound Poisson Gamma process.
Before a run of the CPG model can commence, the file structures must be in the correct format. Included in the package
is an example set of data from Glendalough, Ireland. All other cores must be in a similar format for the routines
to work.
}
\details{
\tabular{ll}{
Package: \tab CPGchron\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2007-07-12\cr
License: \tab GPL (version 2 or later) \cr
}
Installation instructions for basic and first time users:

1. Type library(CPGchron)

2. Now create a directory on your C drive called CPGchron (it needs to be directly on the c drive).

3. Within this directory, create three more called Input, Output and CalCurve.

4. Now navigate to C:\\ program files\\ R\\ R-XXXXX\\ library\\ CPGChron\\ Data. In here there should be a file called Rdata.zip

5. Double click on it and move the files Glendalough.dat, GlendaloughEventDepths.txt and Glendaloughddepths.txt to the input directory.

6. Move the BigCal.txt file to the CalCurve directory.

Everything is now set up.

Example full model run (with Glendalough.dat):

1. At the command prompt in R, type library(CPGchron)

2. Type CPGmenu() and choose option 1.

3. If you've followed the steps above you should not need to change the default path, you just need to tell it that the file name is Glendalough.

4. Now choose option 2 to calibrate the radiocarbon dates.

5. Choose option 3 and short to do a run of the CPG model.

6. Choose option 4 to run the prediction stage.

7. Choose option 6 to produce a plot of the chronology.


Once a satisfactory short run has been obtained, a long run should be undertaken. The long run will take much longer
than the short run, but will only be required once.


Example event prediction stage (with GlendaloughEventDepths.txt)

1. Type CPGmenu() and choose Option 1.

2. Follow step 3 as above to enter the data.

3. Assuming a run of the CPG model has already been done (as above) and that the file coreEventDepths.txt is in the input directory, choose option 5

4. Choose option 7 to produce appropriate plots of the data.

5. Check the output directory for CoreEvents.txt which will contain 95\% HDR age intervals for the depths of interest.


Other cores can be created in a spirit similar to that found in the Glendalough.dat file. This example input file has
5 radiocarbon dates (and the top of the core). The columns are tab
delimited and represent the lab code of the sample, the radiocarbon age, the sample standard error, the depth (in cm)
at which it was found, the thickness of the sample in cm (NOTE: this data is not used in the current version of CPGchron
but it must be included so that the file is read properly.) If the thickness is unknown, zero is acceptable. The last
two columns contain the probabilities of being an outlier. The first probability identifies censored outliers as proposed
by Christen (1994), the second indicates the probability an outlier is ignored completely by the CPGchron model.

}

\author{
Andrew Parnell <Andrew.Parnell@tcd.ie>
}

\references{
Haslett and Parnell (2007, submitted to JRSSC).

Christen, J. A. (1994). Summarizing a set of radiocarbon determinations: A robust approach.
Applied Statistics 43 (3), 489 -- 503.
}
\seealso{
\code{\link{CPGmenu}}
}
\examples{
## See CPGmenu()
}
\keyword{dplot}
\keyword{programming}
\keyword{models}
\keyword{nonparametric}
\keyword{smooth}
\keyword{package}
