// This function runs the main CPGchron malarkey

#include<R.h>
#include<Rmath.h>
#include<stdio.h>
#include<time.h>
#include"use.h"

////////////////////// Some other functions //////////////////////////////////

double Max2 (double a, double b)
{
	// find the max of 2 numbers
   double larger;
   if (a > b)
      larger = a;
   else
      larger = b;
   return larger;
}

//rtruncn function:
double rtruncn (double a, double b)
{
    double A, B;
    double maxA, maxB, maxR, r2, r, th, u, v, x, accept=1.0;
    
    A = atan(a);
    B = atan(b);
    
    maxA = exp(-pow(a,2)/4)/cos(A);
    maxB = exp(-pow(b,2)/4)/cos(B);
    maxR = Max2(maxA, maxB);

    if((a<1) && (b>-1)) maxR = exp(-0.25)*sqrt(2.0);

    while (accept!=0)
    {
        r2 = runif(0.0,1.0);
        r = sqrt(r2)*maxR;
        th = runif(A,B);
        u = r*cos(th);
        v = r*sin(th);
        x = tan(th);
        accept = ((pow(x,2)) < (log(u)*-4));
    }
    return x;

}        

//truncated normal function:
double truncatedwalk (double old, double sd, double low, double high)
{
    double lowlimold, upplimold, y, newvalue;
    lowlimold = (low - old)/sd;
    upplimold = (high - old)/sd;
    y = rtruncn(lowlimold, upplimold);
    newvalue = old + sd*y;
           
    return newvalue;
}

//truncated normal ratio function:
double truncatedrat (double old, double sd, double low, double high, double newvalue)
{
    double lowlimold, upplimold, lowlimnew, upplimnew, plowold, puppold, plownew, puppnew, ratio;
    
    lowlimold = (low - old)/sd;
    upplimold = (high - old)/sd;
    lowlimnew = (low - newvalue)/sd;
    upplimnew = (high - newvalue)/sd;
    plowold = pnorm(lowlimold,0.0,1.0,1,0);
    puppold = pnorm(upplimold,0.0,1.0,1,0);
    plownew = pnorm(lowlimnew,0.0,1.0,1,0);
    puppnew = pnorm(upplimnew,0.0,1.0,1,0);
    ratio = (puppold - plowold)/(puppnew - plownew);
    return ratio;        
}

double Max(double *Numbers, int Count)
{
	// Find the maximum of a sequence of numbers
	double Maximum;
	Maximum = Numbers[0];

	for(int i = 0; i < Count; i++)
		if( Maximum < Numbers[i] )
			Maximum = Numbers[i];

	return Maximum;
}

int seq(double from,double to,double len,double sequence[])
{
	// Create a sequence of numbers from 'from' to 'to' of length 'len'
	// Simple huh?
	
	double by = (to-from)/(len-1);
	int i;
	for(i=0;i<len;i++) 
		sequence[i] = from + i*by;

	return(0);

}

double dtweedielogwsmallp(double y, double phi, double power)
{
	// Matches the R function of the same name. Oh, actually it's called 

	double p,a,a1,r,drop=37,logz,jmax,j,cc,wmax,estlogw,oldestlogw;
	int hij,lowj;

    if (power < 1) 
        exit(-99);
	if (power > 2) 
		exit(-99);
    if (phi <= 0)
		exit(-99);
    if (y <= 0)
		exit(-99);
    p = power;
    a = (2 - p)/(1 - p);
    a1 = 1 - a;
    r = -a * log(y) + a * log(p - 1) - a1 * log(phi) - log(2 - p);
    logz = r;
	
    jmax = (pow(y,(2 - p)))/(phi * (2 - p));
    j = Max2(1, jmax);
    cc = logz + a1 + a * log(-a);
    wmax = a1 * jmax;
    estlogw = wmax;
    while (estlogw > (wmax - drop)) 
	{
        j = j + 2;
        estlogw = j * (cc - a1 * log(j));
    }
	
    hij = (int)ceil(j);
    logz = r;
    jmax = pow(y,(2 - power))/(phi * (2 - power));
    j = Max2(1, jmax);
    wmax = a1 * jmax;
    estlogw = wmax;
    while ((estlogw > (wmax - drop)) && (j >= 2)) 
	{
        j = Max2(1, j - 2);
        oldestlogw = estlogw;
        estlogw = j * (cc - a1 * log(j));
    }
    lowj = (int)Max2(1, floor(j));

	double newj[hij-lowj+1];
    seq(lowj, hij,(hij-lowj+1),newj);
    
	double g[hij-lowj+1]; 
	int k;
	for(k=0;k<hij-lowj+1;k++) g[k] = lgamma(newj[k]+1)+lgamma(-a*newj[k]);
	
	double A[hij-lowj+1];
	for(k=0;k<hij-lowj+1;k++) A[k] = r*(double)newj[k]-g[k];
	
	double m=Max(A,hij-lowj+1);
    double we[hij-lowj+1];
	for(k=0;k<hij-lowj+1;k++) we[k] = exp(A[k]-m);
	double sumwe=0;
	for(k=0;k<hij-lowj+1;k++) sumwe+=we[k];
	double logw=log(sumwe)+m;

	return(logw);

}

double dtweedieseriessmallp(double power,double y, double mu, double phi)
{

// This function matches the R function of the same name (with a few dots in it though)

double logw = dtweedielogwsmallp(y,phi,power);
double tau = phi*(power-1)*pow(mu,power-1);
double lambda = pow(mu,2-power)/(phi*(2-power));
double logf = -y/tau-lambda-log(y)+logw;
double f = exp(logf);

return(f);

}

double dtweediep1(double y, double power, double mu, double phi)
{
// Same as my R function
// Calculates the density of a tweedie plus one random variable

double eps = 0.00000001;
double lambda2 = pow(mu,2-power)/(phi*(2-power))-eps;
double alpha = (2-power)/(power-1);
double beta = 1/(phi*(power-1)*pow(mu,power-1));

double mu2 = alpha*lambda2/beta;
double phi2 = (alpha+1)/(pow(lambda2*alpha,(1/(alpha+1)))*pow(beta,(alpha/(alpha+1))));

double fTplus = dtweedieseriessmallp(power,y,mu,phi)
	+(1/eps)*(dtweedieseriessmallp(power,y,mu,phi)
	-dtweedieseriessmallp(power,y,mu2,phi2));

return(fTplus);

}

double UpdateMCMC(double newloglik,double oldloglik,double newval,double oldval,double rat)
{
// Function to update MCMC when given log likelihoods

double u,mh;
u = runif(0.0,1.0);
mh = exp(newloglik - oldloglik)*rat;
if (u < mh)
{
	return(newval);
} else 
{
	return(oldval);
}

}

int fact(int number) {
	int temp;

	if(number <= 1) return 1;

	temp = number * fact(number - 1);
	return temp;
}

double dlogbinom(int x,int n, double p)
{
double dlogbinom = log(fact(n))-log(fact(x))-log(fact(n-x))+x*log(p)+(n-x)*log(1-p);
return(dlogbinom);

}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
////////////////////// Main CPG function //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void cpg(char**CALPATH,char**INFILE,char**OUTFILE,int*ndets,int*m,int*burnin,int*howmany,int*thinby)
{

/////////////////////////////////// CONSTANTS ////////////////////////////////////
// Some constants for reading stuff in
// len is calcurve length
int BigCalSize=26006;
// Parameters to use on outlier variance
double beta1=2.0,beta2=100.0;

///////////////////////////// READ IN CALIBRATION CURVE /////////////////////////////

double BigC14[BigCalSize],BigSigma[BigCalSize];

FILE *CalFile;
int i;

CalFile = fopen(*CALPATH,"r");

if(CalFile==NULL) {
    error("Error: can't open calibration file.\n");
    exit(1);
} else {
    Rprintf("Big calibration file opened successfully.\n");

    // Now read in
    for(i=0;i<BigCalSize;i++)
    {
       fscanf(CalFile,"%lf",&BigC14[i]);                       
       fscanf(CalFile,"%lf",&BigSigma[i]);                                           
    }

    fclose(CalFile);
}	


///////////////////////////// READ IN DETERMINATIONS /////////////////////////////

//enter determinations and their errors - create them as dynamic arrays and enter them
//from a separate file using same method as cal curve:
char labcode[*ndets][50];
double cage[*ndets],sd[*ndets],depth[*ndets],thick[*ndets],outprob1[*ndets],outprob2[*ndets];
 
FILE *dets;

double numb1[*ndets],numb2[*ndets],numb3[*ndets],numb4[*ndets],numb5[*ndets],numb6[*ndets];

dets = fopen(*INFILE,"r");

if(dets==NULL) {
    error("Error: can't open determinations file.\n");
    exit(0);
} else {
    Rprintf("Determinations file opened successfully.\n");

    // First get rid of header
    char temp[100];
    fgets(temp,100,dets);

    // Now read in as ints and then loop again to convert to double
    for(i=0;i<*ndets;i++)
    {
       fscanf(dets,"%s",labcode[i]);   
       fscanf(dets,"%lf",&numb1[i]);                       
       fscanf(dets,"%lf",&numb2[i]);                       
       fscanf(dets,"%lf",&numb3[i]);                       
       fscanf(dets,"%lf",&numb4[i]);                       
       fscanf(dets,"%lf",&numb5[i]);
       fscanf(dets,"%lf",&numb6[i]);
    }

    for(i=0;i<*ndets;i++)
    {
       cage[i] = (double)numb1[i]/1000;                       
       sd[i] = (double)numb2[i]/1000;
       depth[i] = (double)numb3[i]/100;
       thick[i] = (double)numb4[i]/100;
       outprob1[i] = (double)numb5[i];
       outprob2[i] = (double)numb6[i];
    }
    
    Rprintf("Determinations read successfully.\n");
    
    fclose(dets);
}

///////////////////// STARTING VALUES ////////////////////////////

// Create starting values for dates
double thetaall[*ndets],shift1[*ndets],shift1new,priorp1[*ndets],shift2[*ndets];
double shift2new,priorp2[*ndets];
double thetanew[*ndets],thetadiff[*ndets-1],thetanewdiff[*ndets-1];
double thetanewrat,shift1newrat,psinewrat,meannewrat,shift2newrat;
int flag1[*ndets],flag1new,flag2[*ndets],flag2new;
double p=1.2,mean=5.0,meannew,psi=2.0,psinew; // Starting values
double hi = 10000000; // big number to represent infinity

for(i=0;i<*ndets;i++) 
{
    // Estimate starting values for theta, could do better here wiht linearinterp
    thetaall[i] = cage[i];
    flag1[i] = 0;
	shift1[i] = 0;
    flag2[i] = 0;
	shift2[i] = 0;
	priorp1[i] = outprob1[i];
    priorp2[i] = outprob2[i];
}
// depth differences
double depthdiff[*ndets-1];
diff(depth,ndets,depthdiff);

int iter;		// iterations loop int;
int q,q1;		// determinations loop int;

// somewhere to store likelihoods
double piytheta[*ndets],pixtheta[*ndets];
double piyflag1[*ndets],pixflag1[*ndets];
double piyshift1[*ndets],pixshift1[*ndets];
double piyflag2[*ndets],pixflag2[*ndets];
double piyshift2[*ndets],pixshift2[*ndets];
double piytwmean=0.0,pixtwmean=0.0;
double piytwpsi,pixtwpsi=0.0;

// First sort thetas and get differences of all ages
int t;
int thetaall2[*ndets];

for (t=0; t<*ndets; t++) 
{
    thetaall[t]=thetaall[t]*1000000;
    thetaall2[t]=(int)thetaall[t];
}   
qsort (thetaall2, *ndets, sizeof(int),compare);
for (t=0; t<*ndets; t++) 
{
    thetaall[t]=(double)thetaall2[t]/1000000;
}   

diff(thetaall,ndets,thetadiff);

///////////////////////////// MCMC PART /////////////////////////////

// Open up the output file
FILE *parameterfile;
parameterfile = fopen(*OUTFILE,"w");

Rprintf("Total number of iterations required: %i \n",*m);
Rprintf("Burn in: %i \n",*burnin);

// Do some timing
clock_t c0, c1, c2;
c0 = clock();



// Start iterations here 
for (iter=0;iter<*m;iter++)
{  

    // Get a new seed
    GetRNGstate();

    // print out some of the values of iter
    if((iter % *howmany == 0) & (iter>0)) Rprintf("%i \n",iter);

    // Give some update and some estimated time to finish
    if(iter == *howmany) 
    {
        c2 = clock();
        Rprintf("Estimated time to finish is %5.2f minutes or %5.2f seconds \n",(float) (*m/iter)*(c2 - c0)/(60*CLOCKS_PER_SEC),(float) (*m/iter)*(c2 - c0)/(CLOCKS_PER_SEC));
        Rprintf("Iterations so far ... \n");
    }

    // Write everything to files
	for(q=0; q<*ndets; q++)	
		if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%lf ", thetaall[q]);
	for(q=0; q<*ndets; q++)	
		if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%i ", flag1[q]);
    for(q=0; q<*ndets; q++)	
        if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%lf ", shift1[q]);
    for(q=0; q<*ndets; q++)	
		if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%i ", flag2[q]);
    for(q=0; q<*ndets; q++)	
        if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%lf ", shift2[q]);
    if(iter % *thinby == 0 && iter > *burnin) fprintf(parameterfile,"%lf %lf \n", mean,psi);

    ///////////////////////////////// THETAS ////////////////////////////////////////



	// Update thetas by looping through each determination
	for(q1=0; q1<*ndets; q1++)
    {       

        // Need to randomise the order in which they're updated
        // But need to update them all at least once (to initialise pixtheta[q])
        // for each q, so wait 10 iterations before randomising.
        if(iter > 10) {
            q = (int)floor(runif(0.0,(double)*ndets));
        } else { 
            q = q1;
        }
        
		//sample a new value using a truncated random walk:
		for(i=0;i<*ndets;i++) thetanew[i]= thetaall[i];

		if(q==0)
		{
			thetanew[0] = truncatedwalk(thetaall[0],0.1,0.0,thetaall[1]);
			thetanewrat = truncatedrat(thetaall[0],0.1,0.0,thetaall[1],thetanew[0]);
		} else if(q==*ndets-1)
		{
			thetanew[*ndets-1] = truncatedwalk(thetaall[*ndets-1],0.1,thetaall[*ndets-2],26);
			thetanewrat = truncatedrat(thetaall[*ndets-1],0.1,thetaall[*ndets-2],26,thetanew[*ndets-1]);
		} else {
			thetanew[q] = truncatedwalk(thetaall[q],0.1,thetaall[q-1],thetaall[q+1]);		
			thetanewrat = truncatedrat(thetaall[q],0.1,thetaall[q-1],thetaall[q+1],thetanew[q]);
		}
	

		// Difference the new thetas
        diff(thetanew,ndets,thetanewdiff);

		//calculate old likelihood on first iteration:
        if(iter==0) 
        {
           pixtheta[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)+5]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)+5],2)),1);
        
            for(i=0;i<*ndets-1;i++) 
              pixtheta[q] += log(dtweediep1(thetadiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1)));
        }

        //calculate new likelihood:
        piytheta[q] = dnorm(cage[q],BigC14[(int)(thetanew[q]*1000+0.5)+5]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetanew[q]*1000+0.5)+5],2)),1);
       
        for(i=0;i<*ndets-1;i++) 
            piytheta[q] += log(dtweediep1(thetanewdiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1))); 

		//Update the thetas
		thetaall[q] = UpdateMCMC(piytheta[q],pixtheta[q],thetanew[q],thetaall[q],thetanewrat);
		if(thetaall[q] == thetanew[q]) 
        {
           pixtheta[q] = piytheta[q];
           // Difference the thetas again before re-looping
           diff(thetaall,ndets,thetadiff);
        }

	}		


    for(q1=0; q1<*ndets; q1++)
	{	

        // Need to randomise the order in which they're updated
        // But need to update them all at least once (to initialise pixtheta[q])
        // for each q, so wait 10 iterations before randomising.
        if(iter > 10) {
            q = (int)floor(runif(0.0,(double)*ndets));
        } else { 
            q = q1;
        }

		// Sample a new flag
		if(flag1[q]==0)
		{
			flag1new = 1;
		} else {
			flag1new = 0;
		}
		
		if(iter==0) 
		{     
		    pixflag1[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)+5]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)+5],2)),1)
                  +dlogbinom(flag1[q],1,priorp1[q]);		
        }
		
        piyflag1[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)+5]+flag1new*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)+5],2)),1)
              +dlogbinom(flag1new,1,priorp1[q]);		
		
		flag1[q] = (int)UpdateMCMC(piyflag1[q],pixflag1[q],flag1new,flag1[q],1.0);
		if(flag1[q] == flag1new) pixflag1[q] = piyflag1[q];
	}


    for(q1=0;q1<*ndets;q1++)
    {

        // Need to randomise the order in which they're updated
        // But need to update them all at least once (to initialise pixtheta[q])
        // for each q, so wait 10 iterations before randomising.
        if(iter > 10) {
            q = (int)floor(runif(0.0,(double)*ndets));
        } else { 
            q = q1;
        }

        // Get a new shift
        if(iter==0) 
		{
		       pixshift1[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)+5]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)+5],2)),1)
                    +dnorm(shift1[q],0,sqrt(beta1)*sd[q],1);		
        }
		if(flag1[q]==1)
		{
            shift1new = truncatedwalk(shift1[q],sqrt(beta1)*sd[q],-thetaall[q],26.0-thetaall[q]);
			shift1newrat = truncatedrat(shift1[q],sqrt(beta1)*sd[q],-thetaall[q],26.0-thetaall[q],shift1new);
			
			piyshift1[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)+5]+flag1[q]*shift1new+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)+5],2)),1)
       			  +dnorm(shift1new,0,sqrt(beta1)*sd[q],1);
			
			shift1[q] = UpdateMCMC(piyshift1[q],pixshift1[q],shift1new,shift1[q],shift1newrat);
     		if(shift1[q] == shift1new) pixshift1[q] = piyshift1[q];
		}
        
	}
    
    for(q1=0; q1<*ndets; q1++)
	{	

        // Need to randomise the order in which they're updated
        // But need to update them all at least once (to initialise pixtheta[q])
        // for each q, so wait 10 iterations before randomising.
        if(iter > 10) {
            q = (int)floor(runif(0.0,(double)*ndets));
        } else { 
            q = q1;
        }

		// Sample a new big outlier flag
		if(flag2[q]==0)
		{
			flag2new = 1;
		} else {
			flag2new = 0;
		}
		
		if(iter==0) 
		{     
		    pixflag2[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)+5]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)+5],2)),1)
                  +dlogbinom(flag2[q],1,priorp2[q]);		
        }
		
        piyflag2[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)+5]+flag1[q]*shift1[q]+flag2new*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)+5],2)),1)
              +dlogbinom(flag2new,1,priorp2[q]);		
		
		flag2[q] = (int)UpdateMCMC(piyflag2[q],pixflag2[q],flag2new,flag2[q],1.0);
		if(flag2[q] == flag2new) pixflag2[q] = piyflag2[q];
	}


    for(q1=0;q1<*ndets;q1++)
    {

        // Need to randomise the order in which they're updated
        // But need to update them all at least once (to initialise pixtheta[q])
        // for each q, so wait 10 iterations before randomising.
        if(iter > 10) {
            q = (int)floor(runif(0.0,(double)*ndets));
        } else { 
            q = q1;
        }

        // Get a new big outlier shift
        if(iter==0) 
		{
		       pixshift2[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)+5]+flag1[q]*shift1[q]+flag2[q]*shift2[q],sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)+5],2)),1)
                    +dnorm(shift2[q],0,sqrt(beta2)*sd[q],1);		
        }
		if(flag2[q]==1)
		{
            shift2new = truncatedwalk(shift2[q],sqrt(beta2)*sd[q],-thetaall[q],26.0-thetaall[q]);
			shift2newrat = truncatedrat(shift2[q],sqrt(beta2)*sd[q],-thetaall[q],26.0-thetaall[q],shift2new);
			
			piyshift2[q] = dnorm(cage[q],BigC14[(int)(thetaall[q]*1000+0.5)+5]+flag1[q]*shift1[q]+flag2[q]*shift2new,sqrt(pow(sd[q],2)+pow(BigSigma[(int)(thetaall[q]*1000+0.5)+5],2)),1)
       			  +dnorm(shift2new,0,sqrt(beta2)*sd[q],1);
			
			shift2[q] = UpdateMCMC(piyshift2[q],pixshift2[q],shift2new,shift2[q],shift2newrat);
     		if(shift2[q] == shift2new) pixshift2[q] = piyshift2[q];
		}
        
	}

    ///////////////////////////////// Tweedie ////////////////////////////////////////
    
    // Update tweedie parameters
	// First do mean
	meannew = truncatedwalk(mean,0.5,0.0,hi);
	meannewrat = truncatedrat(mean,0.5,0.0,hi,meannew);
	piytwmean = 0.0;
	
	if(iter==0) {
       pixtwmean=0.0;
	   for(i=0;i<*ndets-1;i++) 
          pixtwmean += log(dtweediep1(thetadiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1)));   
       pixtwmean += dgamma(1/mean,0.01,1/0.01,1);
    }
	
	for(i=0;i<*ndets-1;i++) 
		piytwmean += log(dtweediep1(thetadiff[i],p,meannew*depthdiff[i],psi/pow(depthdiff[i],p-1)));
	piytwmean += dgamma(1/meannew,0.01,1/0.01,1);
	
	mean = UpdateMCMC(piytwmean,pixtwmean,meannew,mean,meannewrat);
	if(mean == meannew) pixtwmean = piytwmean;
	
	// Now update psi
	psinew = truncatedwalk(psi,0.01,0.0,hi);
	psinewrat = truncatedrat(psi,0.01,0.0,hi,psinew);
	piytwpsi = 0.0;

    if(iter==0) {
       pixtwpsi=0.0; 
	   for(i=0;i<*ndets-1;i++) 
          pixtwpsi += log(dtweediep1(thetadiff[i],p,mean*depthdiff[i],psi/pow(depthdiff[i],p-1)));   
       pixtwpsi += dgamma(1/psi,0.01,1/0.01,1);
    }
    
	for(i=0;i<*ndets-1;i++) 
		piytwpsi += log(dtweediep1(thetadiff[i],p,mean*depthdiff[i],psinew/pow(depthdiff[i],p-1)));

	piytwpsi += dgamma(1/psinew,0.01,1/0.01,1);
	
	psi = UpdateMCMC(piytwpsi,pixtwpsi,psinew,psi,psinewrat);
	if(psi == psinew) pixtwpsi = piytwpsi;
	

    // Sort out the RNG state
    PutRNGstate();

    }




fclose(parameterfile);

c1 = clock();
Rprintf("Completed!\n");
Rprintf("Elapsed time in sec: %5.2f\n",(float) (c1 - c0)/CLOCKS_PER_SEC,2);
Rprintf("Elapsed time in minutes: %5.2f\n",(float) (c1 - c0)/(60*CLOCKS_PER_SEC));    
Rprintf("Elapsed time in hours: %5.2f\n",(float) (c1 - c0)/(60*60*CLOCKS_PER_SEC));

}
