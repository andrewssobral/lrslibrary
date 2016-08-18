#include "mex.h"
/* Computes sparse matrix (Y) vector (x) multiplication. 
   Sparse matrix is specified by non-zero elements Ynz and their indices I, J.
   Usage: Yx=compOmegaYx(m,n,Ynz,I,J,x).
   Written by: Prateek Jain (pjain@cs.utexas.edu) and Raghu Meka (raghu@cs.utexas.edu)
   Last updated: 10/29/09
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  mwSize nz_size;
  mwIndex nz_idx;
  int m,n,i,j;
  double *Ynz, *Yx, *x, *I, *J;
  m=*((double *)mxGetPr(prhs[0]));
  n=*((double *)mxGetPr(prhs[1]));
  nz_size=mxGetM(prhs[2]);
  Ynz=mxGetPr(prhs[2]);
  I=mxGetPr(prhs[3]);
  J=mxGetPr(prhs[4]);
  x=mxGetPr(prhs[5]);
  plhs[0]=mxCreateDoubleMatrix(m, 1, mxREAL);
  Yx=mxGetPr(plhs[0]);
  for(i=0;i<m;i++){
    Yx[i]=0;
  }

  for(nz_idx=0;nz_idx<nz_size;nz_idx++){
    Yx[(int)(I[nz_idx])-1]+=Ynz[nz_idx]*x[(int)(J[nz_idx])-1];
  }
}
