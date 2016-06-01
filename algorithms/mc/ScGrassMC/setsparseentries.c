/*=================================================================
% function setsparseentries(M, X)
% M is a sparse matrix (real double).
% X is a vector (real double).
% After this function returns, the nonzero entries of M will be
% equal to those of X.
%
% Compile with: mex setsparseentries.c -largeArrayDims
%
% May 20, 2011, Nicolas Boumal, UCLouvain
 *=================================================================*/

/* #include <math.h> */
#include "string.h"
#include "mex.h"
#include "matrix.h"

/* Input Arguments */

#define	pM	(prhs[0])
#define	pX	(prhs[1])

void mexFunction(
          int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray* prhs[] )
{
    if(!mxIsDouble(pM) || !mxIsDouble(pX))
        mexErrMsgTxt("M and X must be of type DOUBLE.");
    if(!mxIsSparse(pM))
        mexErrMsgTxt("M must be sparse.");
    if(mxGetNzmax(pM) < mxGetNumberOfElements(pX))
        mexErrMsgTxt("M must be able to contain at least as many elements as X.");
    
    memcpy(mxGetPr(pM), mxGetPr(pX), mxGetNumberOfElements(pX) * sizeof(double));
    
    return;
}
