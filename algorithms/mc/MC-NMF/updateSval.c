/* -------------------------------------------------------------------------- */
/* SetSpv mexFunction */
/* -------------------------------------------------------------------------- */

#include "mex.h"
/* #include "blas.h" */
/*#include "matrix.h"*/

/* set values to sparse matrix S */
void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin  [ ]
)
{
    double *Sval, *v, LL;
    unsigned long int L, k;
    ptrdiff_t inc = 1;

    if (nargin != 3 || nargout > 0)
        mexErrMsgTxt ("Usage:  SetSpv ( S, v, L )") ;

    /* ---------------------------------------------------------------- */
    /* inputs */
    /* ---------------------------------------------------------------- */
    
    Sval = mxGetPr( pargin [0] );
    v   = mxGetPr( pargin [1] );
    LL = mxGetScalar( pargin [2] );  L = (unsigned long int) LL;
    
    /* ---------------------------------------------------------------- */
    /* output */
    /* ---------------------------------------------------------------- */

    for (k = 0; k < L; k++) Sval[k] = v[k]; 
    /*dcopy(&L,v, &inc, Sval, &inc);*/
    return;
}
