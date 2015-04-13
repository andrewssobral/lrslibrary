/* setnonzeros(s,v) where s is sparse and numel(v)==nnz(s)
 * returns a new sparse matrix having the structure of s but nonzero values
 * equal to the elements of v(:), in order.
 * Written by Tom Minka
 */
#include "mex.h"
#include "string.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
	const mxArray *s, *v;
	mwSize m,n;
	mwIndex nnz;
  if(nrhs != 2) {
		mexErrMsgTxt("Usage: s = setnonzeros(s,v) where numel(v) == nnz(s)");
  } 
  if(nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }
  if(!(mxIsSparse(prhs[0]))) {
    mexErrMsgTxt("Input argument 1 must be sparse.");
  }
  if(!(mxIsDouble(prhs[0]))) {
    mexErrMsgTxt("Input argument 1 must be of type double.");
  }
  if(mxIsSparse(prhs[1])) {
    mexErrMsgTxt("Input argument 2 must not be sparse.");
  }
  if(!(mxIsDouble(prhs[1]))) {
    mexErrMsgTxt("Input argument 2 must be of type double.");
  }
	s = prhs[0];
	v = prhs[1];
	m = mxGetM(s);
	n = mxGetN(s);
	nnz = mxGetJc(s)[n];
	if(mxGetNumberOfElements(v) != nnz) {
    mexErrMsgTxt("numel(v) != nnz(s)");
	}
	plhs[0] = mxCreateSparse(m,n,nnz,mxIsComplex(v)?mxCOMPLEX:mxREAL);
	memcpy(mxGetPr(plhs[0]), mxGetPr(v), nnz*sizeof(double));
	if(mxIsComplex(v)) memcpy(mxGetPi(plhs[0]), mxGetPi(v), nnz*sizeof(double));
	memcpy(mxGetIr(plhs[0]), mxGetIr(s), nnz*sizeof(mwIndex));
	memcpy(mxGetJc(plhs[0]), mxGetJc(s), (n+1)*sizeof(mwIndex));
}
