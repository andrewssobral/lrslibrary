/* -------------------------------------------------------------------------- */
/* partXY_mex mexFunction */
/* -------------------------------------------------------------------------- */

#include "mex.h"
/*#include "blas.h"*/

/* compute a part of X*Y */
void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin  [ ]
)
{
    double *Xt, *Y, *Z, *I, *J, *v, LL;
    ptrdiff_t m, n, r, L, p, ir, jr, k;
    ptrdiff_t inc = 1;

    if (nargin != 5 || nargout > 1)
        mexErrMsgTxt ("Usage: v = partXY (Xt, Y, I, J, L)") ;

    /* ---------------------------------------------------------------- */
    /* inputs */
    /* ---------------------------------------------------------------- */
    
    Xt = mxGetPr( pargin [0] );
    Y  = mxGetPr( pargin [1] );
    I  = mxGetPr( pargin [2] );
    J  = mxGetPr( pargin [3] );
    LL = mxGetScalar( pargin [4] ); L = (ptrdiff_t) LL;
    m  = mxGetN( pargin [0] );
    n  = mxGetN( pargin [1] );
    r  = mxGetM( pargin [0] ); 
    if ( r != mxGetM( pargin [1] ))
        mexErrMsgTxt ("rows of Xt must be equal to rows of Y") ;
    if ( r > m || r > n )
        mexErrMsgTxt ("rank must be r <= min(m,n)") ;
    
    /* ---------------------------------------------------------------- */
    /* output */
    /* ---------------------------------------------------------------- */

    pargout [0] = mxCreateDoubleMatrix(1, L, mxREAL);
    v = mxGetPr( pargout [0] );
    
    /* C array indices start from 0 */
    for (p = 0; p < L; p++) {
        ir = ( I[p] - 1 ) * r;
        jr = ( J[p] - 1 ) * r;
         v[p] = 0.0;
        for (k = 0; k < r; k++)
            v[p] += Xt[ ir + k ] * Y[ jr + k ];
        /*v[p] = ddot(&r, Xt+ir, &inc, Y+jr, &inc);*/
    }
    
    return;
}

