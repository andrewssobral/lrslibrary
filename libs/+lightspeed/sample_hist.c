/* compile with: cmex sample_hist.c util.o -lm
 *
 * Each column of p is a distribution which is sampled from n times.
 * Requires all(col_sum(p) == 1).
 */
#include "mex.h"
#include "util.h"
#include <string.h>

void MultiRand(double *p, mwSize len, mwSize n, double *result)
{
  double z = 1;
  mwSize i;
  for(i=1;i<len;i++) {
    int r = BinoRand(*p/z, n);
    *result++ = r;
    n -= r;
    z -= *p;
    if(z == 0) {
      memset(result, 0, (len-i)*sizeof(double));
      return;
    }
    p++;
  }
  *result = n;
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  mwSize rows, cols, n, i;
  double *indata, *outdata;

  if(nrhs != 2)
    mexErrMsgTxt("Usage: h = sample_hist(p, n)");

  /* prhs[0] is first argument.
   * mxGetPr returns double*  (data, col-major)
   * mxGetM returns int  (rows)
   * mxGetN returns int  (cols)
   */
  rows = mxGetM(prhs[0]);
  cols = mxGetN(prhs[0]);
  indata = mxGetPr(prhs[0]);
  if(mxGetNumberOfElements(prhs[1]) != 1)
    mexErrMsgTxt("n is not scalar");
  n = (int)*mxGetPr(prhs[1]);

  if(mxIsSparse(prhs[0]))
    mexErrMsgTxt("Cannot handle sparse matrices.  Sorry.");

  /* plhs[0] is first output */
  plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
  outdata = mxGetPr(plhs[0]);
  for(i=0;i<cols;i++) {
    MultiRand(indata, rows, n, outdata);
    indata += rows;
    outdata += rows;
  }
}

