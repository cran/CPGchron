`CPGmenu` <-
function() {

CPGversion <-"1.0"

cat("------------------------------- \n")
cat(paste("Welcome to CPGchron version", CPGversion, "\n"))
cat(paste("Author: Andrew Parnell, Trinity College Dublin\n"))
cat(paste("Please report bugs to: Andrew.Parnell@tcd.ie\n"))
cat("------------------------------- \n")

PATH <- NULL

EXIT <- FALSE
while(EXIT==FALSE)
{

choices <- c("Read in a CPGchron data file (RUN THIS FIRST)","Calibrate a set of radiocarbon dates (30 seconds - 5 minutes)","Run CPGchron (20 minutes - 12 hours)","Predict ages for the entire core (~30 seconds)","Predict ages for certain depths in the core (~30 seconds)","Produce plots of the entire chronology (~2 mins)","Produce plots of the age for certain depths in the core (~30 seconds)","FIRST TIME USERS AND HELP SYSTEM","Exit")
title <- "The available options are:"
choose <- menu(choices,title = title)


#####################################################################################################

# Section 1
if(choose == 1) {

cat("You have chosen to read in a CPGchron data file. \n")
cat("For this you will need to have created an appropriate file as outlined in the help section. \n")
cat("\n")

cat("Please enter the PATH to the file that contains the data file \n")
cat("eg C:\\cores\\mycore \n")
if(R.Version()$arch=="i386")
{
    cat("Leave blank for the default of C:\\CPGchron \n")
} else {
    cat("Leave blank for the default of ~/CPGchron \n")
}
cat("\n")
PATH <- scan(what="",nlines=1,quiet=TRUE)
if(length(PATH)==0 && R.Version()$arch=="i386") PATH <- "C:\\CPGchron"
if(length(PATH)==0 && R.Version()$arch!="i386") PATH <- "~/CPGchron"
cat("Now please enter the name of the .dat file you wish to run CPG chron from. \n")
name <- scan(what="",nlines=1,quiet=TRUE)
cat("Please enter the full name of the core (leave blank if the same as above). \n")
fullname <- scan(what="",nlines=1,quiet=TRUE,sep="\t")
if(length(fullname)==0) fullname <- name
cat("Thank you.\n\n")

# CHECKS
if(R.Version()$arch=="i386")
{
    if(!file.exists(paste(PATH,"\\CalCurve\\BigCal.txt",sep=""))) 
    {
        cat(paste("Check the path, you have specified it as ",PATH,"\\CalCurve\\BigCal.txt\n",sep=""))
        stop("Calibration curve cannot be found",call.=FALSE)
    } 
    if(!file.exists(paste(PATH,"\\Input\\",name,".dat",sep=""))) 
    {
        cat(paste("Check the path, you have specified it as ",PATH,"\\Input\\",name,".dat\n",sep=""))
        stop("Data input file cannot be read",call.=FALSE)
    } 
} else {
    if(!file.exists(paste(PATH,"/CalCurve/BigCal.txt",sep=""))) 
    {
        cat(paste("Check the path, you have specified it as ",PATH,"/CalCurve/BigCal.txt\n",sep=""))
        stop("Calibration curve cannot be found",call.=FALSE)
    } 
    if(!file.exists(paste(PATH,"/Input/",name,".dat",sep=""))) 
    {
        cat(paste("Check the path, you have specified it as ",PATH,"/Input/",name,".dat\n",sep=""))
        stop("Data input file cannot be read",call.=FALSE)
    } 
}

pause()
}

#####################################################################################################

if(choose == 2) {

cat("Now calibrating dates... \n")

if(length(PATH)==0) stop("You have not entered any data. Make sure you choose option 1 with every run of CPGmenu()",call.=FALSE)

# Read in the data file to get the number of determinations
if(R.Version()$arch=="i386")
{
    Temp <- read.table(paste(PATH,"\\Input\\",name,".dat",sep=""),header=TRUE)
} else {
    Temp <- read.table(paste(PATH,"/Input/",name,".dat",sep=""),header=TRUE)
}
ndet <- nrow(Temp)

# Put in the path for the calibration curve
if(R.Version()$arch=="i386")
{
    calpath <- paste(PATH,"\\CalCurve\\BigCal.txt",sep="")
    infile <- paste(PATH,"\\Input\\",name,".dat",sep="")
    outfile <- paste(PATH,"\\Output\\",name,"TrueDates.txt",sep="")
} else {
    calpath <- paste(PATH,"/CalCurve/BigCal.txt",sep="")
    infile <- paste(PATH,"/Input/",name,".dat",sep="")
    outfile <- paste(PATH,"/Output/",name,"TrueDates.txt",sep="")
}
calibrate(calpath,infile,outfile,ndet)


cat("\n \n \n")
pause()

}

#####################################################################################################

if(choose == 3) {

cat("Setting up CPG model... \n")

if(length(PATH)==0) stop("You have not entered any data. Make sure you choose option 1 with every run of CPGmenu()",call.=FALSE)

# Read in the data file to get the number of determinations
if(R.Version()$arch=="i386")
{
    Temp <- read.table(paste(PATH,"\\Input\\",name,".dat",sep=""),header=TRUE)
} else {
    Temp <- read.table(paste(PATH,"/Input/",name,".dat",sep=""),header=TRUE)
}
ndet <- nrow(Temp)


# Put in the path for the calibration curve
if(R.Version()$arch=="i386")
{
    calpath <- paste(PATH,"\\CalCurve\\BigCal.txt",sep="")
    infile <- paste(PATH,"\\Input\\",name,".dat",sep="")
    outfile <- paste(PATH,"\\Output\\",name,"pars.txt",sep="")
} else {
    calpath <- paste(PATH,"/CalCurve/BigCal.txt",sep="")
    infile <- paste(PATH,"/Input/",name,".dat",sep="")
    outfile <- paste(PATH,"/Output/",name,"pars.txt",sep="")
}
choices2 <- c("short","long","super long")

choose2 <- menu(choices2,title="Would you like a short or long run of the model?")

iterations <- 100000
burnin <- 10000
howmany <- 2000
thinby <- 8
if(choose2==2) {
    iterations <- 1000000
    burnin <- 200000
    howmany <- 20000
    thinby <- 75
}
if(choose2==3) {
    iterations <- 10000000
    burnin <- 2000000
    howmany <- 20000
    thinby <- 750
}


CPG(calpath,infile,outfile,ndet,iterations,burnin,howmany,thinby)

if(choose2 > 1) {
cat("\n Now checking convergence of CPGchron parameters... \n")
pars <- read.table(outfile)
pars <-  pars[,c(seq(1,ndet),c(ncol(pars)-1,ncol(pars)))]

pvals <- boa.geweke(pars,0.1,0.5)[,2]
bad <- pvals<0.01
vbad <- pvals<0.001

if(sum(vbad[!is.na(vbad)])>0) {
cat("SEVERE WARNING: it looks like",sum(vbad),"of the parameters have not converged. \n")
} else if(sum(bad[!is.na(bad)])>0) {
cat("WARNING:",sum(bad),"parameters may not have converged. \n")
cat("Problem does not appear to be fatal. \n")
}
if(sum(bad[!is.na(bad)])>0) {
    cat("Plotting trace plots and densities. If these look unsatisfactory, do a longer run. \n")
    bad[is.na(bad)] <- FALSE
    badpars <- as.matrix(pars[,bad==TRUE])
    for(i in 1:ncol(badpars)) {
        if(i!=1) windows()
        par(mfrow=c(2,1))
        plot(badpars[,i],main="Trace plot")
        plot(density(badpars[,i]),main="Density plot")
    }
    
    par(mfrow=c(1,1))

}

}

cat("\n \n \n")
pause()

}

#####################################################################################################

if(choose == 4) {

cat("Now predicting ages for the entire core \n")

if(length(PATH)==0) stop("You have not entered any data. Make sure you choose option 1 with every run of CPGmenu()",call.=FALSE)

# Read in the data file to get the number of determinations
if(R.Version()$arch=="i386")
{
    Temp <- read.table(paste(PATH,"\\Input\\",name,".dat",sep=""),header=TRUE)
} else {
    Temp <- read.table(paste(PATH,"/Input/",name,".dat",sep=""),header=TRUE)
}
ndet <- nrow(Temp)

# Read in the ddepths file to get the number of ddepths
if(R.Version()$arch=="i386")
{
    Temp2 <- read.table(paste(PATH,"\\Input\\",name,"ddepths.txt",sep=""),header=FALSE)
} else {
    Temp2 <- read.table(paste(PATH,"/Input/",name,"ddepths.txt",sep=""),header=FALSE)
}
nddepths <- nrow(Temp2)

# Put in the paths for the various files
if(R.Version()$arch=="i386")
{
    infile <- paste(PATH,"\\Input\\",name,".dat",sep="")
    parsfile <- paste(PATH,"\\Output\\",name,"pars.txt",sep="")
    outfile <- paste(PATH,"\\Output\\",name,"chrons.txt",sep="")
    ddepthfile <- paste(PATH,"\\Input\\",name,"ddepths.txt",sep="")
    outlierfile <- paste(PATH,"\\Output\\",name,"Outliers.txt",sep="")
} else {
    infile <- paste(PATH,"/Input/",name,".dat",sep="")
    parsfile <- paste(PATH,"/Output/",name,"pars.txt",sep="")
    outfile <- paste(PATH,"/Output/",name,"chrons.txt",sep="")
    ddepthfile <- paste(PATH,"/Input/",name,"ddepths.txt",sep="")
    outlierfile <- paste(PATH,"/Output/",name,"Outliers.txt",sep="")
}

extract <- (1950-as.numeric(substr(date(),21,24)))/1000
if(Temp[1,4]!=0) {
    cat("Input date of core extraction in k cal yrs BP (leave blank for default = ",(1950-as.numeric(substr(date(),21,24))),"):\n",sep="")
    extract <- scan(what="",nlines=1,quiet=TRUE)
    if(length(extract)==0) extract <- (1950-as.numeric(substr(date(),21,24)))/1000
}

numchron <- 10000

predictCPG(parsfile,infile,outfile,ndet,ddepthfile,nddepths,numchron,extract,outlierfile)

cat("Completed!\n")
cat("\n \n \n")
pause()

}

#####################################################################################################

if(choose == 5) {

cat("Please create a file of the desired depths you wish to create ages for in the input directory. \n")
cat("Make sure it is called MyCoreEventDepths.txt where 'MyCore' is the name of the core \n")
pause()

if(length(PATH)==0) stop("You have not entered any data. Make sure you choose option 1 with every run of CPGmenu()",call.=FALSE)

# Read in the data file to get the number of determinations
if(R.Version()$arch=="i386")
{
    Temp <- read.table(paste(PATH,"\\Input\\",name,".dat",sep=""),header=TRUE)
} else {
    Temp <- read.table(paste(PATH,"/Input/",name,".dat",sep=""),header=TRUE)
}
ndet <- nrow(Temp)

# Read in the ddepths file to get the number of ddepths
if(R.Version()$arch=="i386")
{
    Temp2 <- read.table(paste(PATH,"\\Input\\",name,"EventDepths.txt",sep=""),header=FALSE)
} else {
    Temp2 <- read.table(paste(PATH,"/Input/",name,"EventDepths.txt",sep=""),header=FALSE)
}
nevents <- nrow(Temp2)

# Put in the paths for the various files
if(R.Version()$arch=="i386")
{
    infile <- paste(PATH,"\\Input\\",name,".dat",sep="")
    parsfile <- paste(PATH,"\\Output\\",name,"pars.txt",sep="")
    outfile <- paste(PATH,"\\Output\\",name,"eventages.txt",sep="")
    evfile <- paste(PATH,"\\Input\\",name,"EventDepths.txt",sep="")
    outlierfile <- paste(PATH,"\\Output\\",name,"Outliers.txt",sep="")
} else {
    infile <- paste(PATH,"/Input/",name,".dat",sep="")
    parsfile <- paste(PATH,"/Output/",name,"pars.txt",sep="")
    outfile <- paste(PATH,"/Output/",name,"eventages.txt",sep="")
    evfile <- paste(PATH,"/Input/",name,"EventDepths.txt",sep="")
    outlierfile <- paste(PATH,"/Output/",name,"Outliers.txt",sep="")
}

extract <- -57/1000

numchron <- 10000

predictCPG(parsfile,infile,outfile,ndet,evfile,nevents,numchron,extract,outlierfile)

cat("Completed!\n")
cat("\n \n \n")
pause()

}

#####################################################################################################

if(choose == 6) {

cat("Now plotting chronology. \n")

if(length(PATH)==0) stop("You have not entered any data. Make sure you choose option 1 with every run of CPGmenu()",call.=FALSE)

# Create arguments
if(R.Version()$arch=="i386")
{
    infile <- paste(PATH,"\\Output\\",name,"chrons.txt",sep="")
    ddepthfile <- paste(PATH,"\\Input\\",name,"ddepths.txt",sep="")
    detsfile <- paste(PATH,"\\Input\\",name,".dat",sep="")
    datesfile <- paste(PATH,"\\Output\\",name,"TrueDates.txt",sep="")
    rangesfile <- paste(PATH,"\\Output\\",name,"Ranges.txt",sep="")
} else {
    infile <- paste(PATH,"/Output/",name,"chrons.txt",sep="")
    ddepthfile <- paste(PATH,"/Input/",name,"ddepths.txt",sep="")
    detsfile <- paste(PATH,"/Input/",name,".dat",sep="")
    datesfile <- paste(PATH,"/Output/",name,"TrueDates.txt",sep="")
    rangesfile <- paste(PATH,"/Output/",name,"Ranges.txt",sep="")
}

choices3 <- c("colour","black and white")
choose3 <- menu(choices3,title="How would you like the chronology to be drawn?")

cols <- ifelse(choose3==1,TRUE,FALSE)

PlotCPG(fullname,infile,ddepthfile,detsfile,datesfile,rangesfile,cols,CPGversion)
cat("Completed!\n")
cat("\n \n \n")
pause()

}

#####################################################################################################

if(choose == 7) {

cat("To use this option, you must have previously run option 5 \n")
cat("\n")
cat(paste("Now plotting date estimates for each of the depths in ",name,"eventages.txt \n",sep=""))

if(length(PATH)==0) stop("You have not entered any data. Make sure you choose option 1 with every run of CPGmenu()",call.=FALSE)

# Create arguments
if(R.Version()$arch=="i386")
{
    infile <- paste(PATH,"\\Output\\",name,"eventages.txt",sep="")
    depthfile <- paste(PATH,"\\Input\\",name,"EventDepths.txt",sep="")
    outfile <- paste(PATH,"\\Output\\",name,"EventHDRs.txt",sep="")
} else {
    infile <- paste(PATH,"/Output/",name,"eventages.txt",sep="")
    depthfile <- paste(PATH,"/Input/",name,"EventDepths.txt",sep="")
    outfile <- paste(PATH,"/Output/",name,"EventHDRs.txt",sep="")
}
PlotDens(infile,depthfile,fullname,outfile,CPGversion)

cat("\n")
cat("Done! \n")
pause()

}

#####################################################################################################

if(choose == 8) {

cat("Thank you for downloading CPGchron, a (reasonably) user-friendly program for producing \n")
cat("reliably dated radiocarbon depth chronologies in Windows. If you are ever stuck whilst using this \n")
cat("program, you can always type help(CPGmenu) at the command prompt for more information.\n")
cat("\n")
cat("WARNING: Some of these tasks, especially a long run of the CPG model, may take several hours \n")
cat("Be aware before you start a long run that this may occur. It is often sensible to run these \n")
cat("overnight. The good thing is you only have to do it once for each core. \n")
cat("\n")
cat("To start, it is recommended for basic users to create a folder called CPGchron at the \n")
cat("root directory of the hard disk. Within, you should create three directories, called inputs, \n")
cat("outputs and calcurve. The inputs directory should contain all the information about the cores for \n")
cat("which you require chronologies. This should include:\n")
cat("1. 'core.dat', a tab-delimited file with columns for the laboratory code of the sample, the \n")
cat("   14C age, the error, the depth (in cm), the thickness (in cm), the probability of \n")
cat("   being a standard outlier, the probability of being an extreme outlier and the type of date information. \n")
cat("   Note that there are three allowed types: 1) a standard radiocarbon date, 2) a calendar age with a Normally \n")
cat("   distributed error, 3) a calendar age with a Uniformly distributed error. For types 1) and 2), the error given \n")
cat("   in column 3 is a 1-sigma standard error. For type 3), the error is the distance to the upper or lower bound \n")
cat("   of the desired Uniform distibution. The format of the .dat file can \n")
cat("   be copied from the supplied Glendalough.dat file. Note that the two outlier columns can \n")
cat("   generally be left as is, though you may want to set some of them to 0 if the supplied dates are \n")
cat("   not radiocarbon ages.\n")
cat("2. 'coreddepths.txt', a single column file which contains each of the desired depths, for \n")
cat("   example, where pollen was collected.\n")
cat("3. (optional) 'coreeventdepths.txt', a single column file containing depths at which \n")
cat("   events of specific interest may have occurred. CPGchron allows the ages of these to \n")
cat("   be looked at individually.\n")
cat("\n")
cat("Finally, the supplied 'BigCal.txt' file should be placed in the CalCurve directory. \n")
cat("Note that CPGchron only uses (at present) the Northern hemisphere 2004 Intcal calibration \n")
cat("curve. Other calibration curves are not supported.\n")
cat("\n")
cat("With all this setup, it should be possible to follow the instructions after typing CPGmenu() \n")
cat("at the command prompt. If you find any bugs, or wish to suggest enhancements, please contact \n")
cat("the author at Andrew.Parnell@tcd.ie \n")
pause()
}

if(choose == 9)
{
cat("Thank you. Exiting... \n")
EXIT=TRUE
}


}

}
