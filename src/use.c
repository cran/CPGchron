// A couple of files that more than one thingy uses

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
