/* solve_triu.c
 * Written by Tom Minka
 */
#include "mex.h"
#include <string.h>

/* How to write MEX files that call LAPACK:
 * http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_external/f13120.html#f45091
 */
/* In R2009a, the blas signatures have been changed to take ptrdiff_t instead of int. We call this BLAS64. */
/* various ways to detect R2009a:
 * new definitions in matrix.h:
 * mxCreateScalarDouble
 * MATHWORKS_MATRIX_MATRIX_PUB_FWD_H
 * MATHWORKS_MATRIX_MXARRAY_PUB_FWD_H
 * MATHWORKS_MATRIX_VAGUE_MXARRAY_HPP
 * existence of blascompat32.h
 */
#ifdef BLAS64
#include "blas.h"
#else
#ifdef UNDERSCORE_LAPACK_CALL
/* Thanks to Ruben Martinez-Cantin */
extern int dtrsm_(char *side, char *uplo, char *transa, char *diag, 
		  int *m, int *n, double *alpha, double *a, int *lda, 
		  double *b, int *ldb);
#else
extern int dtrsm(char *side, char *uplo, char *transa, char *diag, 
		  int *m, int *n, double *alpha, double *a, int *lda, 
		  double *b, int *ldb);
#endif
#endif

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  mwSize m,n;
	ptrdiff_t m64, n64;
#ifndef BLAS64
	int im,in;
#endif
  double *T,*b,*x;
  char side='L',uplo='U',trans='N',diag='N';
  double one = 1;

  if(nrhs != 2 || nlhs > 1)
    mexErrMsgTxt("Usage: x = solve_triu(T,b)");

  /* prhs[0] is first argument.
   * mxGetPr returns double*  (data, col-major)
   * mxGetM returns int  (rows)
   * mxGetN returns int  (cols)
   */
  /* m = rows(T) */
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  if(m != n) mexErrMsgTxt("matrix must be square");
  /* n = cols(b) */
  n = mxGetN(prhs[1]);
  T = mxGetPr(prhs[0]);
  b = mxGetPr(prhs[1]);

  if(mxIsSparse(prhs[0]) || mxIsSparse(prhs[1])) {
    mexErrMsgTxt("Sorry, can't handle sparse matrices yet.");
  }
  if(mxGetNumberOfDimensions(prhs[0]) != 2) {
    mexErrMsgTxt("Arguments must be matrices.");
  }
  if(mxGetNumberOfDimensions(prhs[1]) != 2) {
    mexErrMsgTxt("Arguments must be matrices.");
  }

  /* plhs[0] is first output */
  /* x is same size as b */
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  x = mxGetPr(plhs[0]);
  /* copy b into x to speed up memory access */
  memcpy(x,b,m*n*sizeof(double));
  b = x;
#ifdef BLAS64
	m64 = m;
	n64 = n;
  dtrsm(&side,&uplo,&trans,&diag,&m64,&n64,&one,T,&m64,x,&m64);
#else
	im = (int)m;
	in = (int)n;
#ifdef UNDERSCORE_LAPACK_CALL
  dtrsm_(&side,&uplo,&trans,&diag,&im,&in,&one,T,&im,x,&im);
#else
  dtrsm(&side,&uplo,&trans,&diag,&im,&in,&one,T,&im,x,&im);
#endif
#endif
}

#if 0
  /* Upper triangular */
  for(j=0;j<n;j++) x[m-1 + m*j] = b[m-1 + m*j]/T[m*m - 1];
  for(i=m-2;i>=0;i--) {
    for(j=0;j<n;j++) {
      double s = 0;
      for(k=i+1;k<m;k++) {
	s += T[i + m*k]*x[k + m*j];
      }
      x[i + m*j] = (b[i + m*j] - s)/T[i + m*i];
    }
  }    
  /* Lower triangular */
  for(j=0;j<n;j++) x[m*j] = b[m*j]/T[0];
  for(i=1;i<m;i++) {
    for(j=0;j<n;j++) {
      double s = 0;
      for(k=0;k<i;k++) {
	s += T[i + m*k]*x[k + m*j];
      }
      x[i + m*j] = (b[i + m*j] - s)/T[i + m*i];
    }
  }
#endif
