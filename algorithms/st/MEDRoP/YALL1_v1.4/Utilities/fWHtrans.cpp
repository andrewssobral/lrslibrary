/*=================================================================
 *
 * \file fWHtrans.cpp
 *
 *
 * This code computes the (real) fast discrete Walsh-Hadamard transform with sequency order according to the K.G. Beauchamp's book -- Applications of Walsh and Related Functions.
 *
 *
 * This file is written by Chengbo Li from Computational and Applied Mathematics Department of Rice University.
 *
 *
 * This is a MEX-file for MATLAB.
 *
 * 02/15/2010
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"
#include "matrix.h"
//#include <stdio.h>
//#include <stdlib.h>
// #include <malloc.h>
// #include <stack>


//! Matlab entry function
/*!
 * \param nlhs number of left-hand-side output arguments
 * \param plhs mxArray of output arguments
 * \param nrhs number of right-hand-side input arguments
 * \param prhs mxArray of input arguments
 */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
        
{
    int p, nStage, L, clm;
    int i, j, m, n, N, J, K, M;
    double *v_in, *v_out, *v_ext, *v_final, *temp;
    
    /* Check for proper number of arguments */
    if (nrhs != 1) {
        mexErrMsgTxt("Only one input arguments required.");
    }
    else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Get the size and pointers to input data. */
    m  = mxGetM(prhs[0]);
    n  = mxGetN(prhs[0]);
    
    
    /* Make sure that both input and output vectors have the length with 2^p where p is some integer. */
    p = (int)ceil(log2(m));
    N = 1<<p; /* pow(2, p); */
    
    plhs[0] = mxCreateDoubleMatrix(N, n, mxREAL);
    
    v_final = mxGetPr(plhs[0]);
    v_in = mxGetPr(prhs[0]);
    
    /* Extend the input vector if necessary. */
    v_ext = (double*) mxCalloc(N, sizeof(double));
    v_out = (double*) mxCalloc(N, sizeof(double));
    
    for (clm=0; clm<n; clm++) {
        
        /* C is row major while Matlab is column major */
        for (j=0; j<m; j++){
            v_ext[j] = (v_in+clm*m)[j];
        }
        for (j=m; j<N; j++){
            v_ext[j] = 0;
        }
        
        
        for (i=0; i<N-1; i = i+2) {
            v_ext[i] = v_ext[i] + v_ext[i+1];
            v_ext[i+1] = v_ext[i] - v_ext[i+1] * 2;
        }
        L = 1;
        
        /* main loop */
        for (nStage = 2; nStage<=p; ++nStage){
            M = 1<<L; /* pow(2, L); */
            J = 0;
            K = 0; // difference between Matlab and C++
            while (K<N-1){
                for (j = J;j < J+M-1; j = j+2){
                    /* sequency order  */
                    v_out[K] = v_ext[j] + v_ext[j+M];
                    v_out[K+1] = v_ext[j] - v_ext[j+M];
                    v_out[K+2] = v_ext[j+1] - v_ext[j+1+M];
                    v_out[K+3] = v_ext[j+1] + v_ext[j+1+M];
                    
                    K = K+4;
                }
                
                J = J+2*M;
            }
            
            temp = v_ext;
            v_ext = v_out;
            v_out = temp;
            
            L = L+1;
        }
        
        /* Perform scaling of coefficients. */
        for ( i =0; i<N; ++i){
            (v_final+clm*N)[i] = v_ext[i]/N;
        }
    }
    
    /* Set free the memory. */
    mxFree(v_out);
    mxFree(v_ext);
    
    return;
}
