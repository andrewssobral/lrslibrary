/* 
 * Stephen Becker, 11/10/08
 * Updates a sparse vector very quickly
 * calling format:
 *      updateSparse(Y,b)
 * which updates the values of Y to be b
 *
 * Modified 11/12/08 to allow unsorted omega
 * (omega is the implicit index: in Matlab, what
 *  we are doing is Y(omega) = b. So, if omega
 *  is unsorted, then b must be re-ordered appropriately 
 *
 * Modified May '09 to allow complex entries
 * Need to make sure this works on 64-bit systems! use mwSize, mwIndex
 *
 * Modified July '09 to check that input sparse matrix contains
 * numeric entries, not doubles.  Compiled versions of this are
 * not yet distributed.
 *
 * */

#include "mex.h"


/* to use "mwSize", need matrix.h, but doesn't work for R2006a */
/* #include "matrix.h" */
/* So, use the following definitions instead: */
#ifndef mwSize
    #define mwSize size_t
#endif
#ifndef mwIndex
    #define mwIndex size_t  /* should make it compatible w/ 64-bit systems */
#endif


#ifndef true
    #define true 1
#endif
#ifndef false
    #define false 0
#endif

void printUsage() {
    mexPrintf("usage:\tupdateSparse(Y,b)\nchanges the sparse Y matrix");
    mexPrintf(" to have values b\non its nonzero elements.  Be careful:\n\t");
    mexPrintf("this assumes b is sorted in the appropriate order!\n");
    mexPrintf("If b (i.e. the index omega, where we want to perform Y(omega)=b)\n");
    mexPrintf("  is unsorted, then call the command as follows:\n");
    mexPrintf("\tupdateSparse(Y,b,omegaIndx)\n");
    mexPrintf("where [temp,omegaIndx] = sort(omega)\n");
}

void mexFunction(
         int nlhs,       mxArray *plhs[],
         int nrhs, const mxArray *prhs[]
         )
{
    /* Declare variable */
    mwSize M, N;
    mwIndex i, j, m, n;
    double *b, *S, *omega;
    int SORTED = true;
    int COMPLEX;
    
    /* Check for proper number of input and output arguments */    
    if ( (nrhs < 2) || (nrhs > 3) )  {
        printUsage();
        mexErrMsgTxt("Needs 2 or 3 input arguments");
    } 
    if ( nrhs == 3 ) SORTED = false;
    if(nlhs > 0){
        printUsage();
        mexErrMsgTxt("No output arguments!");
    }
    
    /* Check data type of input argument  */
    if (!(mxIsSparse(prhs[0])) || !((mxIsDouble(prhs[1]))) ){
        printUsage();
        mexErrMsgTxt("Input arguments wrong data-type (must be sparse, double).");
    }   
    
    /* check if "b" is complex */
    if ( mxIsComplex(prhs[1]) ) {
        COMPLEX = 1;
       /*  mexPrintf("Input is complex\n");  */
        if (!mxIsComplex( prhs[0] )) {
            printUsage();
            mexErrMsgTxt("if second input is complex, first input must be also");
        }
    }else{
        COMPLEX = 0;
        if (mxIsComplex( prhs[0] )) {
            printUsage();
            mexErrMsgTxt("if first input is complex, second input must be also");
        }
    }

    /* Get the size and pointers to input data */
    /* Check second input */
    N = mxGetN( prhs[1] );
    M = mxGetM( prhs[1] );
    if ( (M>1) && (N>1) ) {
        printUsage();
        mexErrMsgTxt("Second argument must be a vector");
    }
    N = (N>M) ? N : M;

    
    /* Check first input */
    M = mxGetNzmax( prhs[0] );
    if ( M != N ) {
        printUsage();
        mexErrMsgTxt("Length of second argument must match nnz of first argument");
    }

    /* Check that the sparse matrix contains doubles, not logicals
     * or something else */
    if (!mxIsDouble( prhs[0] )) {
        printUsage();
        mexErrMsgTxt("First input must contain numeric entries (double floats), not logicals");
    }

    /* if 3rd argument provided, check that it agrees with 2nd argument */
    if (!SORTED) {
       m = mxGetM( prhs[2] );
       n = mxGetN( prhs[2] );
       if ( (m>1) && (n>1) ) {
           printUsage();
           mexErrMsgTxt("Third argument must be a vector");
       }
       n = (n>m) ? n : m;
       if ( n != N ) {
           printUsage();
           mexErrMsgTxt("Third argument must be same length as second argument");
       }
       omega = mxGetPr( prhs[2] );
    }


    b = mxGetPr( prhs[1] );
    S = mxGetPr( prhs[0] );

    if (SORTED) {
        /* And here's the really fancy part:  */
        for ( i=0 ; i < N ; i++ )
            S[i] = b[i];
    } else {
        for ( i=0 ; i < N ; i++ ) {
            /* this is a little slow, but I should really check
             * to make sure the index is not out-of-bounds, otherwise
             * Matlab could crash */
            j = (int)omega[i]-1; /* the -1 because Matlab is 1-based */
            if ((j >= N)||(j<0)){
                printUsage();
                mexErrMsgTxt("Third argument must have values < length of 2nd argument");
            }
/*             S[ j ] = b[i]; */  /* this is incorrect */
            S[ i ] = b[j];  /* this is the correct form */
        }
    }
    
    if (COMPLEX){
        b = mxGetPi( prhs[1] );
        S = mxGetPi( prhs[0] );
        if (SORTED) {
            for ( i=0 ; i < N ; i++ )
                S[i] = b[i];
        } else {
            for ( i=0 ; i < N ; i++ ) {
                j = (int)omega[i]-1; 
                S[ i ] = b[j]; 
            }
        }
    }
     
}
