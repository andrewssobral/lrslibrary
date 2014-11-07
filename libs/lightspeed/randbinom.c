/* compile with: 
      Windows: mex randbinom.c util.obj
      Others:  cmex randbinom.c util.o -lm
 */
#include "mexutil.h"
#include "util.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  mwSize ndims, *dims, len, n, i;
  double *indata, *outdata;

  if(nrhs != 2)
    mexErrMsgTxt("Usage: r = bino_sample(p, n)");

  /* prhs[0] is first argument.
   * mxGetPr returns double*  (data, col-major)
   */
  ndims = mxGetNumberOfDimensions(prhs[0]);
  dims = (mwSize*)mxGetDimensions(prhs[0]);
  indata = mxGetPr(prhs[0]);
  len = mxGetNumberOfElements(prhs[0]);

  if(mxGetNumberOfElements(prhs[1]) != 1)
    mexErrMsgTxt("n is not scalar");
  n = (int)*mxGetPr(prhs[1]);

  if(mxIsSparse(prhs[0]))
    mexErrMsgTxt("Cannot handle sparse matrices.  Sorry.");

  /* plhs[0] is first output */
  plhs[0] = mxCreateNumericArrayE(ndims, dims, mxDOUBLE_CLASS, mxREAL);
  outdata = mxGetPr(plhs[0]);
  for(i=0;i<len;i++) {
    *outdata++ = BinoRand(*indata++, n);
  }
}

