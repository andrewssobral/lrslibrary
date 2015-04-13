#include "mex.h"
#include "matrix.h"

/* compile with -largeArrayDims flag */

/* to use "mwSize", need matrix.h, but doesn't work for R2006a */
/* So, use the following definitions instead: */
#ifndef mwSize
    #define mwSize size_t
#endif
#ifndef mwIndex
    #define mwIndex size_t  /* should make it compatible w/ 64-bit systems */
#endif


void printUsage() {
    mexPrintf("usage:\tsmvp(A,b), where A is a sparse m x n matrix\n");
    mexPrintf("\tand b is a dense n x 1 vector.  Handles real or complex data\n");
}


void mexFunction( int nlhs,       mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] 
                 )
                
{   /* Sparse matrix-dense vector product 
     *       
     *  Inputs:
     *           
     *           A - m x n sparse matrix
     *           b - n x 1 dense vector
     *           
     *  Outputs:
     *   
     *           x = A*b
     *
     *  Darren Engwirda 2006.
     *
     *  Modified by Stephen Becker, April 2009
     *  Notes on compiling:
     *  compile with       
     *      mex -O -v -largeArrayDims smvp.c
     *
     *  Modified by Stephen Becker, May 2009, to allow complex multiplies
     */

    double *s, *b, *x, rhs;
    double *si, *bi, *xi, rhsi;  /* the imaginary parts, if data is complex */
    /*int *ir, *jc, i, ncol, k;    */
    mwIndex *ir, *jc, i, k, stop;
    mwSize ncol;
    int COMPLEX=0, B_REAL=1;
    
    
    /* Check I/O number */
    if (nlhs>1) {
        printUsage();
        mexErrMsgTxt("SMVP: Too many outputs");
    }
    if (nrhs!=2) {
        printUsage();
        mexErrMsgTxt("SMVP: Incorrect number of inputs");
    }
    
    /* Brief error checking */
    
    ncol = mxGetN(prhs[0]);
    
    if ((ncol!=mxGetM(prhs[1])) || (mxGetN(prhs[1])!=1)) {
        printUsage();
        mexErrMsgTxt("Wrong input dimensions");
    }
    if (!mxIsSparse(prhs[0])) {
        printUsage();
        mexErrMsgTxt("Matrix must be sparse");
    }
    
    /* handle complex data -- for now, both A and b must be the same type */
    if ( mxIsComplex( prhs[0] ) ){
        COMPLEX = 1; B_REAL =0;
        if (!mxIsComplex( prhs[1] )) {
            /*
            printUsage();
            mexErrMsgTxt("If A is complex, b must be also");
             **/
            /* this situation arise a lot, so let's handle it */
            B_REAL = 1;
        }
    } else if (mxIsComplex( prhs[1] )) {
        printUsage();
        mexErrMsgTxt("If b is complex, A must be also");
    }
    
    
    
    /* Allocate output */
    if (COMPLEX)
        plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), 1, mxCOMPLEX);
    else
        plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), 1, mxREAL);
    
    /* I/O pointers */
    ir = mxGetIr(prhs[0]);      /* Row indexing      */
    jc = mxGetJc(prhs[0]);      /* Column count      */
    s  = mxGetPr(prhs[0]);      /* Non-zero elements */
    b  = mxGetPr(prhs[1]);      /* Rhs vector        */
    x  = mxGetPr(plhs[0]);      /* Output vector     */
    if (COMPLEX) {
        si  = mxGetPi(prhs[0]);    
        xi  = mxGetPi(plhs[0]);    
        if (!B_REAL) {
            bi = mxGetPi( prhs[1] );
        }
    }
    
    
    /* Multiplication */
    if (COMPLEX) {
        for (i=0; i<ncol; i++) {            /* Loop through columns */
            stop = jc[i+1];
            rhs = b[i];
            if ( B_REAL )
                rhsi = 0.0;
            else
                rhsi = bi[i];
            for (k=jc[i]; k<stop; k++) {    /* Loop through non-zeros in ith column */
                x[ir[k]] += s[k] * rhs - si[k] * rhsi;
                xi[ir[k]]+= si[k]* rhs + s[k] * rhsi;
            }
        }
    } else {
        for (i=0; i<ncol; i++) {            /* Loop through columns */
            stop = jc[i+1];
            rhs = b[i];
            for (k=jc[i]; k<stop; k++) {    /* Loop through non-zeros in ith column */
                x[ir[k]] += s[k] * rhs;
            }
        }
    }

}
