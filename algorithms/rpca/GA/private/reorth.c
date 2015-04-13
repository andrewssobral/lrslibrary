/***********************************************************************
%% mexreorth.c : C mex file 
%%
%% unix: 
%% mex -O -largeArrayDims -lmwlapack -lmwblas mexreorth.c
%%
%%   [R_NEW,NORMR_NEW,NRE] = reorth(Q,R,NORMR,INDEX,ALPHA,METHOD)
%%   reorthogonalizes R against the subset of columns of Q given by INDEX. 
%%   If INDEX==[] then R is reorthogonalized all columns of Q.
%%   If the result R_NEW has a small norm, i.e. if norm(R_NEW) < ALPHA*NORMR,
%%   then a second reorthogonalization is performed. If the norm of R_NEW
%%   is once more decreased by  more than a factor of ALPHA then R is 
%%   numerically in span(Q(:,INDEX)) and a zero-vector is returned for R_NEW.
%%
%%   If method==0 then iterated modified Gram-Schmidt is used.
%%   If method==1 then iterated classical Gram-Schmidt is used.
%%
%%   The default value for ALPHA is 0.5. 
%%   NRE is the number of reorthogonalizations performed (1 or 2).
%% 
%% NNLS, version 0: 
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun 
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

/********************************
* realdot: x dense vector,  
*          y dense vector
*********************************/
double realdot(const double *x, const double *y, const int n)

{  int i;
   double r; 

   r=0.0;
   for (i=0; i<n-3; i++) {             /* LEVEL 4 */
      r += x[i] * y[i]; i++;
      r += x[i] * y[i]; i++;
      r += x[i] * y[i]; i++;
      r += x[i] * y[i]; }
   if (i<n-1) {                        /* LEVEL 2 */
      r += x[i] * y[i]; i++;
      r += x[i] * y[i]; i++; }
   if (i<n) {                          /* LEVEL 1 */
      r += x[i] * y[i]; }
   return r; 
}
/**********************************************************
* modified Gram-Schmidt: 
* Orthogalizes R against the k vectors in Q by the 
* iterative process 
***********************************************************/

void mgs(mwSize n, double *Q, double *R, mwIndex *index, mwSize k)  

{ mwSize i,j,idx,idxn;
  double s;   

    for (i=0; i<k; i++) { 
      idx = index[i]; 
      idxn = idx*n; 
      s = 0.0; 
      for (j=0; j<n; j++) { s += Q[j+idxn]*R[j]; }
      for (j=0; j<n; j++) { R[j] -= s*Q[j+idxn]; } 
    }
}

/**********************************************************
* 
***********************************************************/
void mexFunction(
      int nlhs,   mxArray  *plhs[], 
      int nrhs,   const mxArray  *prhs[] )

{    double   *Q, *R, *Rout, *indextmp, *normrout, *nreout;  
    
     mwIndex  subs[2];    
     mwSize   nsubs=2; 
     mwIndex  *index; 
     mwSize   n, k1, k, j, nre, method;
     double   alpha, normr, normrold; 

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs < 3){
      mexErrMsgTxt("mexreorth: requires at least 2 input arguments."); }
   if (nlhs > 4){ 
      mexErrMsgTxt("mexreorth: requires at most 4 output argument."); }   

/* CHECK THE DIMENSIONS */

    n  = mxGetM(prhs[0]); 
    k1 = mxGetN(prhs[0]); 
    if (mxIsSparse(prhs[0])) {
       mexErrMsgTxt("mexreorth: sparse Q not allowed."); }   
    Q = mxGetPr(prhs[0]);     
    if (mxIsSparse(prhs[1])) {
       mexErrMsgTxt("mexreorth: sparse R not allowed."); }   
    if (mxGetM(prhs[1])!=n || mxGetN(prhs[1])!=1) { 
       mexErrMsgTxt("mexreorth: R is not compatible."); }  
    R = mxGetPr(prhs[1]);     
    if (nrhs < 3) { 
       normr = sqrt(realdot(R,R,n)); 
    } else {
       normr = *mxGetPr(prhs[2]); 
    }
    if (nrhs < 4) {
       k = k1;  
       index = mxCalloc(k,sizeof(mwSize)); 
       for (j=0; j<k; j++) { index[j]=j; } 
    } else {
       indextmp = mxGetPr(prhs[3]);  
       if (indextmp == NULL) { 
          mexErrMsgTxt("mexreorth: index is empty"); }        
       k = MAX(mxGetM(prhs[3]),mxGetN(prhs[3]));
       index = mxCalloc(k,sizeof(mwSize));  
       for (j=0; j<k; j++) { index[j] = (mwSize) (indextmp[j]-1); }
    }
    if (nrhs< 5) { alpha = 0.5; } else { alpha = *mxGetPr(prhs[4]); }
    if (nrhs==5) { method = (mwSize)*mxGetPr(prhs[5]); } else { method = 0; } 

    /***** create return argument *****/
    plhs[0]  = mxCreateDoubleMatrix(n,1,mxREAL); 
    Rout     = mxGetPr(plhs[0]);  
    plhs[1]  = mxCreateDoubleMatrix(1,1,mxREAL); 
    normrout = mxGetPr(plhs[1]);  
    plhs[2]  = mxCreateDoubleMatrix(1,1,mxREAL); 
    nreout   = mxGetPr(plhs[2]);  

    memcpy(mxGetPr(plhs[0]),mxGetPr(prhs[1]),(n)*sizeof(double));

    /***** main  *****/    
    normrold = 0; 
    nre = 0; 
    while ((normr<alpha*normrold) || (nre==0)) {
       mgs(n,Q,Rout,index,k);       
       normrold = normr;       
       normr = sqrt(realdot(Rout,Rout,n)); 
       nre = nre+1;
       if (nre > 4) { 
	 /** printf("R is in span(Q)"); **/
         for (j=0; j<n; j++) { Rout[j] = 0; }  
	 normr = 0; 
         break; 
       }
    }   
    normrout[0] = normr; 
    nreout[0] = nre; 
    return;
 }
/**********************************************************/
