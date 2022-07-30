#include <stdlib.h>
#include <stdio.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"


/*
-------------------------- Function epp1 -----------------------------

             x= sign(v) max( |v|- rho, 0)
 
Usage (in matlab)
 x=epp1(v, rho); 

-------------------------- Function epp1 -----------------------------
 */
 
void  epp1(double *x, double *v, int n, double rho)
{
	int i;
	
	for(i=0;i<n;i++)
	{ 
		if (fabs(v[i])<=rho)
			x[i]=0;
    else if (v[i]< -rho)
      x[i]=v[i]+rho;
    else
		  x[i]=v[i]-rho;
  }
}


void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

    double *v, *x;
    int m, n, numel;
    double rho;
    
    if (nrhs != 2 && nlhs != 1)
    {
 			mexErrMsgTxt("Invalid arguments.");
    }
     
    /*set up input arguments */
   	v = mxGetPr(prhs[0]);
    numel = mxGetNumberOfElements(prhs[0]);
    rho = mxGetScalar(prhs[1]);
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);    
    x=mxGetPr(plhs[0]);
    
    epp1(x, v, numel, rho);
}