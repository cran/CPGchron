// A couple of files that more than one thingy uses
#include<R.h>
#include<Rmath.h>

int diff(double *arr,int *len,double *retarr)
{
// this function takes a one-dimensional array arr and its length len, and returns the differenced
// vector retarr of length len-1

int i;
for(i=0;i<*len-1;i++)
{
	retarr[i] = arr[i+1]-arr[i];
}
return(0);

}

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

////////////////////// Some other functions //////////////////////////////////


int GetLengthCurrentDepths(double depthlow,double depthhigh,double ddepth[],int length)
{
// This function finds the sizes of the array required to store the current depths

// Need to loop through the depths to get the number of things in ddepths which are
// between depthlow and depthhigh

int i,count=0;
for(i=0;i<length;i++)
{
    if((ddepth[i]<=depthhigh) & (ddepth[i]>=depthlow)) count = count+1;

}

return(count);
}

int GetCurrentDepths(double depthlow,double depthhigh,double ddepth[],int length,double current[])
{
// This function takes currentdepths and fills it with all the ddepths between depthlow and
// depthhigh

int i,count=0;
for(i=0;i<length;i++)
{
    if((ddepth[i]<=depthhigh) & (ddepth[i]>=depthlow))
    {
      current[count] = ddepth[i];
      count = count+1;  
    }
    
}

return(0);
}

int GetCurrentDepthRows(double depthlow,double depthhigh,double ddepth[],int length,int rows[])
{
// This function takes currentdepths and fills it with all the ddepths between depthlow and
// depthhigh

int i,count=0;
for(i=0;i<length;i++)
{
    if((ddepth[i]<=depthhigh) & (ddepth[i]>=depthlow))
    {
      rows[count] = i;
      count = count+1;  
    }
    
}

return(0);
}

double linearinterp(int n, double newx, double *a, double *b)
{
    double newvalue;
    int i;

//condition is where y lies between the two closests approximations to 
//it in the cal curve
   
    for(i=0; i<n-1; i++)
    {
        if (((newx >= a[i]) & (newx <= a[i+1])) | ((newx <= a[i]) & (newx >= a[i+1])))
        {
                newvalue = b[i] + ((newx-a[i])/(a[i+1]-a[i]))*(b[i+1]-b[i]);
                if(newx==a[i]) newvalue = b[i];
                return(newvalue);
                //break;
        }        
    }
  
  return(-999.0);
}
