/*************************************************************************
 * function B = inplacecolumnmex(A, k)
 * Return the inplace-column A(:,k)
 * Important notes: 
 * - use MEX function releaseinplace(B) to release properly shared-data
 *   pointer before clear/reuse B.
 * - All inplace variables shared data with A must be released before
 *   the original array A is cleared/reused.
 * Thanks to James Tursa
 ************************************************************************/
#include "mex.h"
#include "matrix.h"

/* Uncomment this on older Matlab version where size_t has not been 
 defined */
/*
#define mwSize int
#define size_t int
 */

/* The following file defines the internal representation of mxArray,
 * inspired from mxArray_tag declared in the header <matrix.h>.
 * This file is built by calling the MATLAB function
 * buildInternal_mxArrayDef.m */
#include "Internal_mxArray.h"

/* Gateway of inplacecolumnmex */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    mwSize k, N, M;
    double *Pr;
    
    /* Check arguments */
    if (nrhs!=2)
        mexErrMsgTxt("INPLACECOLUMN: Two input arguments required.");
    
    if (!mxIsNumeric(prhs[0]))
        mexErrMsgTxt("INPLACECOLUMN: First input A argument must be numeric.");
    
    if (!mxIsNumeric(prhs[1]))                              
        mexErrMsgTxt("INPLACECOLUMN: Second input K must be numeric.");

    /* Get the size */
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    
    /* Get the column number k from the second input */
    if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1)
        mexErrMsgTxt("INPLACECOLUMN: Second input K must be a scalar.");
    
    if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
        mexErrMsgTxt("INPLACECOLUMN: Second input K must be a double.");
    
    k = (mwSize)(*mxGetPr(prhs[1]));
    /* Make sure k is valid */
    if (k<1 || k>N)
        mexErrMsgTxt("INPLACECOLUMN: K is not valid.");
    
    /* Create the Matrix result (first output) */
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxFree(mxGetPr(plhs[0])); /* Free the data, normally Pr is NULL and
                               * this call does nothing */

    /* Set the dimension as one column */
    mxSetM(plhs[0], M);
    mxSetN(plhs[0], 1);
    
    /* Inplace data pointer of A */
    Pr = mxGetPr(prhs[0]);
    Pr += (k-1)*M; /* Point to the column #k */
    /* Equivalent to doing this: mxSetPr(plhs[0], Pr); 
       but access directly to data pointer in order to by pass Matlab
       checking */
    ((Internal_mxArray*)(plhs[0]))->data.number_array.pdata = Pr;
      
    return;

} /* Gateway of INPLACECOLUMNMEX.c */

