/*************************************************************************
 * function releaseinplace(b)
 * Release the data from an inplace column mxArray that was created
 * with the inplacecolumnmex function.
 *  Author Bruno Luong <brunoluong@yahoo.com>
 *  Last update: 27/June/2009
 ************************************************************************/

#include "mex.h"
#include "matrix.h"

/* Uncomment this on older Matlab version where size_t has not been 
 defined */
/*#define size_t int*/

/* The following file defines the internal representation of mxArray,
 * inspired from mxArray_tag declared in the header <matrix.h>.
 * This file is built by calling the MATLAB function
 * buildInternal_mxArrayDef.m */
#include "Internal_mxArray.h"

/* Gateway of releaseinplace */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    /* Check arguments */
    if (nrhs!=1)
        mexErrMsgTxt("RELEASEINPLACE: One input argument required.");
    
    if( nlhs != 0 ) {
        mexErrMsgTxt("RELEASEINPLACE: Zero output arguments required.");
    }
    
    mxSetM((mxArray *)prhs[0], 0);
    mxSetN((mxArray *)prhs[0], 0);
    /* Equivalent to doing this:  mxSetPr(prhs[0], NULL); 
       but access directly to data pointer in order to by pass Matlab
       checking - Thanks to James Tursa */
    ((Internal_mxArray*)(prhs[0]))->data.number_array.pdata = NULL;
    ((Internal_mxArray*)(prhs[0]))->data.number_array.pimag_data = NULL; 
      
    return;

} /* Gateway of releaseinplace.c */ 