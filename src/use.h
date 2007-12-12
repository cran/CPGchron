// Header file for use.c

int diff(double *arr,int *len,double *retarr);
int compare (const void * a, const void * b);
int GetLengthCurrentDepths(double depthlow,double depthhigh,double ddepth[],int length);
int GetCurrentDepths(double depthlow,double depthhigh,double ddepth[],int length,double current[]);
int GetCurrentDepthRows(double depthlow,double depthhigh,double ddepth[],int length,int rows[]);
double linearinterp(int n, double newx, double *a, double *b);
