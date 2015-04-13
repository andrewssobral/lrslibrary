/* 
MEX interface for LAPACK routine bdsqr.
Matlab calling sequence:
  [sigma,bnd] = bdsqr(alpha,beta)   

  Part of PROPACK by R. Larsen
*/

/* Stephen Becker, 11/12/08
 * Now, dbdqr is in a C file.  So, no underscore for the name.
 * Also, for Windows, the Lapack libraries don't need underscores
 * either.  These are controlled by pre-processor definitions
 * that I put in; either edit this source file, or, better, pass
 * in a -D[name of variable].
 * e.g. mex ... -DWINDOWS ...
 *
 * 11/6/09
 * Noticed major problems with 64-bit systems
 * Download dbdsqr.f from netlib and use in place of MATLAB's LAPACK version
 * Compile as follows (linux)
 * mex bdsqr_mex.c dbdqr.c dbdsqr.f -DDBDQR_IN_C -output bdsqr -llapack -largeArrayDims -lblas
 * or...
 * mex bdsqr_mex.c dbdqr.f dbdsqr.f -UDBDQR_IN_C -output bdsqr -llapack -largeArrayDims -lblas
 *
 * on Windows, to compile dbdsqr.f, need a Fortran compiler
 * See http://www.codingday.com/compile-lapack-and-blas-as-dll-on-windows/
 * and the compiler:
 *  http://sourceforge.net/projects/mingw/files/
 * (once compiled to dbdsqr.o, can distribute)
 * 
 * 11/9/09
 * Defining DBDQR_IN_C by default
 * Switching function prototypes to use ptrdiff_t instead of int
 * Adding "extern" to definitions
 * Couldn't compile fortran on windows (mingw works well, but mangles names
 * in a way that's not compatible with the MS Visual C++ compilter)
 * However, by switching to ptrdiff_t in the definitions, and using mwlapack library,
 * it now works on it's own, so no need to use dbdsqr.f explicitly. Problem solved!
 * On Linux, compile as:
 *  mex -v bdsqr_mex.c dbdqr.c -output bdsqr -llapack -largeArrayDims -lblas
 *
 *  12/4/09
 *  use blas.h if its available
 *
 * */


#include <stdio.h>
#include <string.h>
#include "mex.h"

#ifdef _WIN32
    #define WINDOWS
#endif
#ifndef FORTRAN_WRAPPER
    #if defined(WINDOWS) || defined(__hpux)
    #define FORTRAN_WRAPPER(x) x
    #else
    #define FORTRAN_WRAPPER(x) x ## _
    #endif
#endif

/* SRB, Oct'09  -- actually, is this necessary?*/
#ifndef NO_MATRIX_H
    #include "matrix.h"
#endif

/* if we include lapack.h, it redefines lapack symbols to have the appropriate underscore! nice */
#include "lapack.h"

/* if lapack.h is not included, we need to worry about underscores ourselves */
#ifndef dbdsqr
    #define dbdsqr FORTRAN_WRAPPER(dbdsqr)
#endif
/* if lapack.h is not included, here is the prototype for the lapack routinen we need: */
/*extern void dbdsqr(char *uplo, ptrdiff_t *n, ptrdiff_t *ncvt, ptrdiff_t *nru, ptrdiff_t *ncc,
	     double *d, double *e, double *vt, ptrdiff_t *ldt, double *u,
	     ptrdiff_t *ldu, double *c, ptrdiff_t *ldc, double *work, ptrdiff_t *info);
*/

/* this is a helper file */
extern void dbdqr(int *n, double *d, double *e, double *c1, double *c2);

/* the gateway function to be called by Matlab: */
void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
  int m, n, i, info, zero=0, one=1;
  ptrdiff_t M, Zero = 0, One = 1, Info;
  double dummy=1, *wrk, *bnd, *d, *e;

  if (nrhs != 2)
     mexErrMsgTxt("bdsqr requires two input arguments");
  else if  (nlhs != 2)
     mexErrMsgTxt("bdsqr requires two output arguments");

  m = mxGetM(prhs[0]); /* get the dimensions of the input */
  n = mxGetN(prhs[0]);
  M = m; 

  /* make sure input input vectors are same length */
  if (m != mxGetM(prhs[1]) )
    mexErrMsgTxt("alpha and beta must have the same size");

  /* make sure input is m x 1 */
  if ( n != 1 || mxGetN(prhs[1]) != 1 || n != mxGetN(prhs[1])) 
    mexErrMsgTxt("alpha and beta must be a m x 1 vectors");
    
  /* Create/allocate return arguments */
  for (i=0; i<2; i++) {  
    plhs[i]=mxCreateDoubleMatrix(m,1,mxREAL);  
    if ( plhs[i] == NULL )
        mexErrMsgTxt("unable to allocate output");
  } 

  e = mxCalloc(m,sizeof(double)); 
  wrk = mxCalloc(4*m,sizeof(double)); 
  if ( (e==NULL) || (wrk==NULL) ) /* SRB adding Oct;09 */
      mexErrMsgTxt("Failed to allocate memory");
  d = mxGetPr(plhs[0]); 

  memcpy(d,mxGetPr(prhs[0]), m*sizeof(double));
  memcpy(e,mxGetPr(prhs[1]), m*sizeof(double));
  bnd = mxGetPr(plhs[1]); /* automatically zeroed out */

  /* Reduce to upper m-by-m upper bidiagonal */
  dbdqr(&m, d, e, &bnd[m-1],&dummy);

  /* Compute singular values and last row of U */
  dbdsqr("U",&M,&Zero,&One,&Zero, d, e, &dummy, &One, 
       bnd, &One, &dummy, &One, wrk, &Info);  
  info = (int) Info;

  /* Check exit status of dbdsqr */
  if ( info < 0 )
    mexErrMsgTxt("DBDSQR was called with illegal arguments");
  else if ( info > 0)
    mexWarnMsgTxt("DBDSQR: singular values did not converge");

  /* Free work arrays */
  mxFree(e); 
  mxFree(wrk);
}
