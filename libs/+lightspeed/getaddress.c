/* Written by Tom Minka
 * (c) Microsoft Corporation. All rights reserved.
 */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
  double *p;
  if(nrhs != 1) {
    mexErrMsgTxt("usage: getaddress(a)");
  }
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  p = mxGetPr(plhs[0]);
	*p = (double)(long)mxGetData(prhs[0]);
}
