/* Adapted from partXY.c of Wotao Yin*/
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *U, *V;
  int *I, *J;
  double *out;
  mwIndex m, n, r, l, p, ir, jr, k;

  if (nrhs != 4 || nlhs > 1)
    mexErrMsgTxt("Usage: out = maskmult(U,V,I,J)");
  if (mxGetClassID(prhs[0])!=mxDOUBLE_CLASS || mxGetClassID(prhs[1])!=mxDOUBLE_CLASS) {
    mexErrMsgTxt("First and second inputs must be double type");
  }
  if (!mxIsInt32(prhs[2]) || !mxIsInt32(prhs[3])) {
    mexErrMsgTxt("Third and fourth inputs must be int32 type");
  }

  U = mxGetPr(prhs[0]);
  V = mxGetPr(prhs[1]);
  m = mxGetN(prhs[0]);
  r = mxGetM(prhs[0]);
  n = mxGetN(prhs[1]);

  if (r!=mxGetM(prhs[1])) {
    mexErrMsgTxt("Rows of U must be equal to rows of V");
  }

  I = (int*) mxGetData(prhs[2]);
  J = (int*) mxGetData(prhs[3]);

  if (mxGetN(prhs[2])!=1 || mxGetN(prhs[3])!=1) {
    mexErrMsgTxt("I and J must be column vectors");
  }
  if (mxGetM(prhs[2])!=mxGetM(prhs[3])) {
    mexErrMsgTxt("I and J must be of the same length");
  }

  l = mxGetM(prhs[2]);
  
  plhs[0] = mxCreateDoubleMatrix(1, l, mxREAL);
  out = mxGetPr(plhs[0]);
  
  for (p=0; p<l; p++) {
    ir = (I[p]-1)*r;
    jr = (J[p]-1)*r;
    out[p] = 0.0;
    for (k=0; k<r; k++) {
      out[p] += U[ir+k]*V[jr+k];
    }
  }
}

