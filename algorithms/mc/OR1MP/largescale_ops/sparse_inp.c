/* ---------------------------------------------------------- */
/* mexFunction: sparse_inp_native                             */
/*                                                            */
/* inner-product and apply sparse selection                   */
/* for array simplicity, we follow LMaFit to use transpose U  */
/*                                                            */
/* Jiayu Zhou (jiayu.zhou@asu.edu)                            */
/* ---------------------------------------------------------- */

#include "mex.h"

void mexFunction
(
    int nargout,                  /*number of output variables*/
    mxArray *pargout [ ],         /*pointer of output variables*/
    int nargin,                   /*number of input variables*/
    const mxArray *pargin  [ ]    /*pointer of input variables*/
)
{
    /* define variables */
    double *Ut, *V, *iR, *iC, *v, LL;
    ptrdiff_t m, n, r;       /*matrix information*/
    size_t   len_r,  len_c;  
    size_t   ii, ir, jr, k;  /*local variables for computation*/

    if (nargin != 4 || nargout > 1)
        mexErrMsgTxt ("Usage: out = sparse_inp_native (U', V, iR, iC)") ;

    /* get inputs from mex */
    Ut  = mxGetPr( pargin [0] );
    V   = mxGetPr( pargin [1] );
    iR  = mxGetPr( pargin [2] );
    iC  = mxGetPr( pargin [3] );
    len_r = mxGetNumberOfElements(pargin [2]);
    len_c = mxGetNumberOfElements(pargin [3]);
    
    /* get matrix information */
    m   = mxGetN( pargin [0] );
    n   = mxGetN( pargin [1] );
    r   = mxGetM( pargin [0] ); 
    
    /* validation input */
    if ( r != mxGetM( pargin [1] ))
        mexErrMsgTxt ("the columns of U must the same as the rows of V") ;
    if ( r > m || r > n )
        mexErrMsgTxt ("r must be less or equal than min(m,n)") ;
    if (len_r != len_c)
        mexErrMsgTxt ("iR and iC do not have the same length") ;
    
    /* allocate memory for output*/
    pargout [0] = mxCreateDoubleMatrix(1, len_r, mxREAL);
    v = mxGetPr( pargout [0] );
    
    /* compute innner product on selected indices */
    for (ii = 0; ii < len_r; ii++) {
        ir = ( iR[ii] - 1 ) * r;
        jr = ( iC[ii] - 1 ) * r;
        v[ii] = 0.0;
        for (k = 0; k < r; k++)
            v[ii] += Ut[ ir + k ] * V[ jr + k ];
    }
    return;
}

