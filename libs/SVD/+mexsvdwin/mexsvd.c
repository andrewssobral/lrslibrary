/***********************************************************************
* mexsvd.c : compute singular value decomposition of a matrix
*            based on a divide-and-conquer method implemented in the
*            LAPACK routine dgesdd.f
*
* Compilation in 
* unix: mex -O -largeArrayDims -lmwlapack -lmwblas mexsvd.c
*
* [U,S,V,flag] = mexsvd(A,options); A = U*S*VT 
* or
* [S] = mexsvd(A); compute only singular values
* 
* options = 0, (default)
*              economical SVD: 
*              dim(S) = minMN x minMN, dim(U) = M x minMN, dim(VT) = minMN x N.
*         = 1, full SVD: 
*              dim(S) = M x N, dim(U) = M x M, dim(VT) = N x N, 
*
*  flag   = 0:  successful exit.
*         < 0:  if flag = -i, the i-th argument had an illegal value.
*         > 0:  LAPACK DBDSDC did not converge, updating process failed.
*
* NNLS, version 0: 
* Copyright (c) 2009 by
* Kim-Chuan Toh and Sangwoon Yun 
***********************************************************************/

#include <math.h>
#include <mex.h>
#include <matrix.h>
#include <string.h> /* needed for memcpy() */

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define  MIN(A, B)   ((A) < (B) ? (A) : (B))
#endif

/**********************************************************
* 
***********************************************************/
void mexFunction(
      int nlhs,   mxArray  *plhs[], 
      int nrhs,   const mxArray  *prhs[] )

{    double   *A, *U, *V, *VT, *S, *flag, *work, *AA;  

     mwIndex  subs[2];
     mwSize   nsubs=2; 
     mwIndex  *irS, *jcS, *iwork; 
     mwSize   M, N, lwork, info, options, minMN, maxMN, k, j; 
     mwSize   LDA, LDU, LDVT, nU, mVT, mS, nS; 
     char     *jobz;

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs > 2){
      mexErrMsgTxt("mexsvd: requires at most 2 input arguments."); }
   if (nlhs > 4){ 
      mexErrMsgTxt("mexsvd: requires at most 4 output argument."); }   

/* CHECK THE DIMENSIONS */

    M = mxGetM(prhs[0]); 
    N = mxGetN(prhs[0]); 
    if (mxIsSparse(prhs[0])) {
       mexErrMsgTxt("mexeig: sparse matrix not allowed."); }   
    A = mxGetPr(prhs[0]); 
    LDA = M; 
    minMN = MIN(M,N); 
    maxMN = MAX(M,N); 
    options = 0; 
    if (nrhs==2) { options = (int)*mxGetPr(prhs[1]); } 
    if (options==0) { 
        /***** economical SVD ******/
       jobz="S"; nU = minMN; mVT = minMN; mS = minMN; nS = minMN; 
    } else {
        /***** full SVD ******/
       jobz="A"; nU = M; mVT = N; mS = M; nS = N; 
    }   
    LDU  = M; 
    LDVT = mVT; 
    /***** create return argument *****/
    if (nlhs >=3) {
       /***** compute singular values and vectors ******/
       plhs[0] = mxCreateDoubleMatrix(M,nU,mxREAL); 
       U = mxGetPr(plhs[0]);  
       plhs[1] = mxCreateSparse(mS,nS,minMN,mxREAL); 
       S   = mxGetPr(plhs[1]); 
       irS = mxGetIr(plhs[1]); 
       jcS = mxGetJc(plhs[1]);
       plhs[2] = mxCreateDoubleMatrix(N,mVT,mxREAL);     
       V = mxGetPr(plhs[2]); 
       plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);     
       flag = mxGetPr(plhs[3]); 
    } else { 
       /***** compute only singular values ******/
       plhs[0]= mxCreateDoubleMatrix(minMN,1,mxREAL); 
       S = mxGetPr(plhs[0]); 
       plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);     
       flag = mxGetPr(plhs[1]); 
       U = mxCalloc(M*nU,sizeof(double)); 
       V = mxCalloc(N*mVT,sizeof(double)); 
       jobz="N";
    }
    /***** Do the computations in a subroutine *****/   
    lwork = 4*minMN*minMN + MAX(maxMN,5*minMN*minMN+4*minMN);  
    work  = mxCalloc(lwork,sizeof(double)); 
    iwork = mxCalloc(8*minMN,sizeof(int)); 
    VT = mxCalloc(mVT*N,sizeof(double)); 
    AA = mxCalloc(M*N,sizeof(double)); 
    memcpy(AA,mxGetPr(prhs[0]),(M*N)*sizeof(double));

    dgesdd(jobz,&M,&N, AA,&LDA,S,U,&LDU,VT,&LDVT,work,&lwork,iwork, &info); 

    flag[0] = (double)info; 
    if (nlhs >= 3) { 
       for (k=0; k<minMN; k++) { irS[k] = k; }
       jcS[0] = 0;
       for (k=1; k<=nS; k++) { 
         if (k<minMN) { jcS[k] = k; } else { jcS[k] = minMN; }
       }  
       for (k=0; k<mVT; k++) { 
          for (j=0; j<N; j++) { V[j+k*N] = VT[k+j*mVT]; }
       }
    }
    return;
 }
/**********************************************************/
