/* 
MEX interface for TQLB. Matlab calling sequence:  
  [lambda,top,bot,err] = tqlb(alpha,beta)   
*/


#include <string.h>
#include "mex.h"

/* Template for tqlb: */
void tqlb_(int *n, double *d__, double *e, double *bnd, 
	   double *bnd2, int *ierr);

/* Here comes the gateway function to be called by Matlab: */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
  int m, n,i, ierr;
  double x, *tmp;

  if (nrhs != 2)
     mexErrMsgTxt("tqlb requires two input arguments");
  else if  (nlhs != 4)
     mexErrMsgTxt("tqlb requires four output arguments");

  for (i=0; i<2; i++) { 
    m = mxGetM(prhs[i]); /* get the dimensions of the input */
    n = mxGetN(prhs[i]);
    
    /* make sure input is m x 1 */
    if (n != 1) 
      mexErrMsgTxt("Input must be a m x 1 vectors");
  }

  /* Create/allocate return argument, a 1x1 real-valued Matrix */
  for (i=0; i<3; i++) { 
    plhs[i]=mxCreateDoubleMatrix(m,1,mxREAL); 
  }
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);   
  tmp = mxCalloc(m,sizeof(double));
  
  memcpy(mxGetPr(plhs[0]), mxGetPr(prhs[0]),m*sizeof(double));
  memcpy(tmp,mxGetPr(prhs[1]), m*sizeof(double));
  tqlb_(&m,mxGetPr(plhs[0]),tmp,mxGetPr(plhs[1]),
	mxGetPr(plhs[2]),&ierr);
  
  *(mxGetPr(plhs[3])) = (double) ierr;
  mxFree(tmp);
}
