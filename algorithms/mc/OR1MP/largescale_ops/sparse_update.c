/* -------------------------------------------------------------------------- */
/* mexFunction: sparse_update      */
/* set values to sparse matrix S   */
/* -------------------------------------------------------------------------- */

#include "mex.h"

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin  [ ]
)
{
    /* Declare variables */ 
    size_t elements;
    double *Sval, *v;
    unsigned long int k;

    if (nargin != 2 || nargout > 0)
        mexErrMsgTxt ("Usage: sparse_update ( S, v )") ;

    /* ---------------------------------------------------------------- */
    /* inputs */
    /* ---------------------------------------------------------------- */
    
    Sval = mxGetPr( pargin [0] );
    v    = mxGetPr( pargin [1] );
    elements=mxGetNumberOfElements(pargin [1]);
    
    /* TODO: Validate the length */
   
    for (k = 0; k < elements; k++) 
        Sval[k] = v[k]; 
    return;
}
