/* 
 * Stephen Becker, 11/10/08
 * Computes A(omega), where A = U*V' is never explicitly computed 
 * srbecker@caltech.edu
 *
 * Partial Changelog:
 *
 *  3/24/09
 *      adding in "mxIsComplex" and "mxGetPi" routines to handle complex data seamlessly
 *      Not fully implemented yet, though; need to define .' vs ' convention (i.e. transpose vs adjoint)
 *  5/13/09
 *      Finished implementing the complex stuff.
 *      The convention is A = U*V', not A = U*V.'
 *      Having problems with complex data on 64-bit linux though: zdot isn't working
 *      see
 *      http://matlab.izmiran.ru/help/techdoc/matlab_external/ch04cr17.html
 *      and the fort.h file in this directory, for alternative methods
 *  5/19/09
 *      try mex -largeArrayDims -lblas -DWINDOWS XonOmega.c  ??
 *      it only affects zdotu, not the level-3 blas
 *      tried -fPIC (as opposed to -fpic) in LDFLAGS, no luck
 *      So, my current solution is to implement the blas function myself for 64-bit processors
 *
 * 10/29/09
 *      see XonOmegaTranspose.c for a much faster implementation
 *      The idea is that we pass in U' and V' (instead of U and V)
 *      Very little modification of the code, but benefit is enormous:
 *      instead of stride lengths of M and N, stride lengths are now 1
 *      (because C and MATLAB use column-major order).
 *      Also, adding in new level-3 BLAS feature: if omega isn't sorted,
 *         and even if length(omega) isn't quite N*M, we still do level-3 BLAS
 *         call, and THEN sort it or subsample it.  Also, extend to omegaX,omegaY case
 *      Also, why were omega, omegaI and omegaJ all doubles, not ints, or similar??
 *          Because MATLAB passes them in as doubles.  Do NOT cast as integers!!
 *  11/9/09
 *      Consider including blas.h and using ptrdiff_t as the types, not ints
 *      Added _WIN32 definition (implies WINDOWS) definition
 *      Switched malloc and free to mxMalloc and mxFree
 *
 * 
 * Note to the user:
 *      to really take advantage of this function, you should recompile it on
 *      your computer so that you can take advantage of any good BLAS libraries you might have
 *      However, it is unlikely that this will be the speed bottleneck for most applications...
 *
 *      There are 3 ways to call this program:
 *
 *      (1) b = XonOmega(U,V,omega)
 *              implicity forms the matrix A = U*V', then returns b = A(omega)
 *              where omega are linear indices.
 *              This is perhaps 20% faster if omega is sorted
 *
 *      (2) b = XonOmega(U,V,omegaX,omegaY)
 *              does the same, but omegaX and omegaY are subscripts
 *              (i.e. omega = ind2sub(... OmegaX, OmegaY), or [OmegaX,OmegaY] = sub2ind(...,omega) )
 *              This is sometimes faster than method (1).  If in addition, omegaX and omegaY
 *              are created from the indices of a sorted linear index
 *              (so that omegaY is sorted), then it will be about 20% faster
 *
 *      (3) b = XonOmega(U,V,Y)
 *              gets information about omega from the sparsity pattern of the sparse matrix Y
 *              This is is actually the fastest method usually!
 *
 *      (note: the relative speed of the different calling sequences
 *          will depend on your computer)
 *
 *      Update: all of the above methods to call the program can be modified by adding
 *      an extra argument, cuttoff, to the end of the input sequence, e.g.
 *      (1) b = XonOmega(U,V,omega,cutoff)
 *      (2) b = XonOmega(U,V,omegaX,omegaY,cutoff)
 *      (3) b = XonOmega(U,V,Y,cutoff)
 *      where cutoff is a number between 0 and 1 that affects how the internal
 *      computation is done; specifically, if length(omega) > cutoff*size(U,1)*size(V,1),
 *      then the intermediate matrix A = U*V' is actually computed (using level-3 BLAS,
 *      which is usually very optimized).  The default value of cutoff is .25,
 *      which is roughly the level where you begin to see benefit from not explicitly constructing A.
 *
 * Compilation notes:
 *   For efficieny, this program makes BLAS calls whenever possible;
 *   it uses level-3 BLAS if the entire matrix A = U*V' is needed,
 *   and if only a small subset is needed, then it makes many level-1 BLAS calls.
 *
 *   On Windows computers, BLAS functions have no underscore, but for linux/unix/darwin, they do 
 *   So, for Windows, either uncomment the line that says #define WINDOWS, or alternatively,
 *   when compiling, use mex -L(location) -lmwblas -largeArrayDims -DWINDOWS XonOmega.c
 *      (the -DWINDOWS is the same as putting "#define WINDOWS" in the code;
 *      "location<S-Del>) is fullfile(matlabroot,'extern','lib','win32',cc)
 *       where cc is usually 'lcc' (default) or 'microsoft' (for visual studio);
 *       The location string should be enclosed in quotes)
 *    ex. on my Windows machine, I do this
mex -O -L"C:\Program Files\MATLAB\R2008b\extern\lib\win32\microsoft" ...
        -largeArrayDims -lmwblas -DWINDOWS XonOmega.c
 *    
 * 
 *
 *   For linux, etc., you can specifically undefine the WINDOWS symbol
 *   by passing in the -UWINDOWS option to mex, but that's not
 *   necessary unless you uncommented the #define WINDOWS line below.
 *   You can use the mathwork's blas, -lmwblas, or, recommended,
 *   use a good blas that's already installed, and probably available
 *   with just -lblas.  So, try this:
 *      mex -O -lblas -largeArrayDims updateSparse.c
 *   the -O is for optimized
 *   the -lblas may actually not even be necessary
 *
 *   For both windows and linux, the -largeArrayDims flag is important if you have
 *   a 64-bit system, otherwise you can skip it.
 *
 *   If you've never compiled with mex before, please see the MATHWORKS help page.
 *   You probably want to run "mex -setup" the first time.
 *
 *
 * Dec '09 notes:
 *  zdotu -- if using libmwblas and matlab/extern/.../blas.h exists, use this for prototype and use blas call on 64-bit
 *  also, use fortran_wrapper stuff
 * Source: zdotu.f
#define zdotu FORTRAN_WRAPPER(zdotu)
 #ifndef FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID
 extern doublecomplex zdotu(
 #else
 extern void zdotu(
     doublecomplex* retval,
     #endif
         ptrdiff_t *n,
             double *zx,
                 ptrdiff_t *incx,
                     double *zy,
                         ptrdiff_t *incy
                         );

 *
  extern complex16 zdotu_( int *K, complex16 *x, int *incx, complex16 *y, int *incy );
 * */

#include "mex.h"
/* #include "limits.h"  in mex.h already */
/* #include "blas.h" */  /* might want to try */


/* #define WINDOWS */
/* I use the "WINDOWS" symbol to change between underscores and no underscores
 * _WIN32 *should* be automatically defined on WINDOWS machines, for most compilers */
#if defined(_WIN32)
    #define WINDOWS
#endif


/* to use "mwSize", need matrix.h, but doesn't work for R2006a */
#include "matrix.h" 
/* So, use the following definitions instead: */
#ifndef mwSize
    #define mwSize size_t
#endif
#ifndef mwIndex
    #define mwIndex size_t  /* should make it compatible w/ 64-bit systems */
#endif
/* R2009a: BLAS is now 64-bit on 64-bit systems, so use mwSignedIndex
 * For R2006a and earlier, when matrix.h doesn't exist, the following
 * definitions can be used: */
#ifndef mwSignedIndex
    #define mwSignedIndex ptrdiff_t
#endif
/* However, I'm not changing this just yet...
 * http://www.mathworks.com/access/helpdesk/help/techdoc/index.html?/access/helpdesk/help/techdoc/apiref/mwsignedindex.html&http://www.mathworks.com/access/helpdesk/help/techdoc/rn/brvak9c-1.html  
 *
 *http://www.mathworks.com/access/helpdesk/help/techdoc/index.html?/access/helpdesk/help/techdoc/rn/brvak9c-1.html&http://www.mathworks.com/access/helpdesk/help/techdoc/rn/bryg9vd-1.html#br02jub-1
 * */

#ifndef true
    #define true 1
#endif
#ifndef false
    #define false 0
#endif

typedef struct{ double re; double im; } complex16; /* a hack */
/* this is what mathworks does in lapack.h: */
/*
#ifndef COMPLEX_TYPES
#define COMPLEX_TYPES
  typedef struct{float r,i;} complex;
  typedef struct{double r,i;} doublecomplex;
#endif
  */


#ifdef WINDOWS
  /* level-1 BLAS */
  extern double ddot( int *K, double *x, int *incx, double *y, int *incy );
  /* level-3 BLAS */
  extern void dgemm( char *transA, char *transB, int *M, int *N, int *K, 
          double *alpha, double *A, int *LDA, double *R, int *LDB, 
          double *beta, double *C, int *LDC );
  extern void zgemm( char *transA, char *transB, int *M, int *N, int *K,
          complex16 *alpha, complex16 *A, int *LDA, complex16 *R, int *LDB, 
          complex16 *beta, complex16 *C, int *LDC );
  extern complex16 zdotu( int *K, complex16 *x, int *incx, complex16 *y, int *incy );
#else
  extern double ddot_( int *K, double *x, int *incx, double *y, int *incy );
  extern void dgemm_( char *transA, char *transB, int *M, int *N, int *K,
          double *alpha, double *A, int *LDA, double *R, int *LDB, 
          double *beta, double *C, int *LDC );
  extern void zgemm_( char *transA, char *transB, int *M, int *N, int *K,
          complex16 *alpha, complex16 *A, int *LDA, complex16 *R, int *LDB, 
          complex16 *beta, complex16 *C, int *LDC );
  extern complex16 zdotu_( int *K, complex16 *x, int *incx, complex16 *y, int *incy );
#endif


void printUsage() {
    mexPrintf("XonOmega.c: usage is\n\t b = XonOmega(U,V,omega)\n");
    mexPrintf("where A = U*V' and b = A(omega)\n");
    mexPrintf("\nAlternative usage is\n\t b = XonOmega(U,V,omegaI,omegaJ),\n where [omegaI,omegaJ] = ind2sub(...,omega)\n");
    mexPrintf("\nAnother alternative usage is\n\t b = XonOmega(U,V, OMEGA),\n where OMEGA is a sparse matrix with nonzeros on omega.\nThis will agree with the other forms of the command if omega is sorted\n\n");
    mexPrintf("If U and V are complex, then make sure A=U*V' and not A=U*V.'\n\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
{

/* deal with the underscores now so that they can be ignored after this point */
#ifdef WINDOWS
    double (*my_ddot)( int *, double *, int *, double *, int *) = ddot;
    complex16 (*my_zdot)( int *, complex16 *, int *, complex16 *, int *) = zdotu;
    void (*my_dgemm)(char *, char *, int *, int *, int *, double *, double *,
            int *, double *, int *, double *, double *, int *) = dgemm;
    void (*my_zgemm)(char *, char *, int *, int *, int *, complex16 *, complex16 *,
        int *, complex16 *, int *, complex16 *, complex16 *, int *) = zgemm;
#else
    double (*my_ddot)( int *, double *, int *, double *, int *) = ddot_;
    complex16 (*my_zdot)( int *, complex16 *, int *, complex16 *, int *) = zdotu_;
    void (*my_dgemm)(char *, char *, int *, int *, int *, double *, double *,
            int *, double *, int *, double *, double *, int *) = dgemm_;
    void (*my_zgemm)(char *, char *, int *, int *, int *, complex16 *, complex16 *,
            int *, complex16 *, int *, complex16 *, complex16 *, int *) = zgemm_;
#endif

    mwSize M, N, K, K2;
    mwSize nOmega1, nOmega2, nOmega;
    mwIndex i,j,k,m;
    double *U, *V, *output, *outputBLAS;
    double *U_imag, *V_imag, *output_imag; /* for complex data */
    complex16 *U_cplx, *V_cplx, temp_cplx; /* for complex data */
    double *omega, *omegaX, *omegaY;
    mwIndex *omegaI, *omegaJ;  /* for sparse version */
    int SPARSE = false;
    int COMPLEX = false;
    int USE_BLAS = false;
    int LARGE_BIT = false;
    int MODE_THREE = false;
    complex16 alpha_cplx, beta_cplx;
    complex16 *output_cplx;


    char transA = 'N', transB = 'T';
    mwSize LDA, LDB;
    double alpha, beta;
    double BLAS_CUTOFF;
    
    /* Check for proper number of input and output arguments */    
    if ( (nrhs < 3) || (nrhs > 5) ) {
        printUsage();
        mexErrMsgTxt("Three to five input argument required.");
    } 
    if(nlhs > 1){
        printUsage();
        mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Check data type of input argument  */
    if (!(mxIsDouble(prhs[0])) || !((mxIsDouble(prhs[1]))) ){
        printUsage();
        mexErrMsgTxt("Input arguments wrong data-type (must be doubles).");
    }   

    /* Get the size and pointers to input data */
    M  = mxGetM(prhs[0]);
    K  = mxGetN(prhs[0]);
    N  = mxGetM(prhs[1]);
    K2  = mxGetN(prhs[1]);
    if ( K != K2 ) {
        printUsage();
        mexErrMsgTxt("Inner dimension of U and V' must agree.");
    }
    COMPLEX = (( (mxIsComplex(prhs[0])) ) || (mxIsComplex(prhs[1])) );
    nOmega1 = mxGetM( prhs[2] );
    nOmega2 = mxGetN( prhs[2] );


    /* on 64-bit systems, these may be longs, but I really want
     * them to be ints so that they work with the BLAS calls */
    mxAssert(M<INT_MAX,"Matrix is too large for 32-bit FORTRAN");
    mxAssert(N<INT_MAX,"Matrix is too large for 32-bit FORTRAN");
    mxAssert(K<INT_MAX,"Matrix is too large for 32-bit FORTRAN");
    
    
    if ( (nOmega1 != 1) && (nOmega2 != 1) ) {
/*         mexErrMsgTxt("Omega must be a vector"); */
        /* Update:
         * if this happens, we assume Omega is really a sparse matrix
         * and everything is OK */
        if ( ( nOmega1 != M ) || ( nOmega2 != N ) || (!mxIsSparse( prhs[2] )) ){
            printUsage();
            mexErrMsgTxt("Omega must be a vector or a sparse matrix");
        }
        nOmega = mxGetNzmax( prhs[2] );
        SPARSE = true;
    } else {
        nOmega = nOmega1 < nOmega2 ? nOmega2 : nOmega1;
        if ( nOmega > N*M ) {
            printUsage();
            mexErrMsgTxt("Omega must have M*N or fewer entries");
        }
    }

    U = mxGetPr(prhs[0]);
    V = mxGetPr(prhs[1]);
    if (COMPLEX) {
        U_imag = mxGetPi(prhs[0]);
        V_imag = mxGetPi(prhs[1]);
        plhs[0] = mxCreateDoubleMatrix(nOmega, 1, mxCOMPLEX);
        output = mxGetPr(plhs[0]);
        output_imag = mxGetPi(plhs[0]);

        U_cplx = (complex16 *)mxMalloc( M*K*sizeof(complex16) );
        V_cplx = (complex16 *)mxMalloc( N*K*sizeof(complex16) );
        if ( (U_cplx == NULL) || ( V_cplx == NULL ) ) {
            /* this should be very rare */
            mexErrMsgTxt("Unable to allocate memory.  Matrix too large?");
        }

        for ( i=0 ; i < M*K ; i++ ){
            U_cplx[i].re = (double)U[i];
            U_cplx[i].im = (double)U_imag[i];
        }
        for ( i=0 ; i < N*K ; i++ ){
            V_cplx[i].re = (double)V[i];
            V_cplx[i].im = -(double)V_imag[i]; 
            /* minus sign, since adjoint, not transpose */
        }
        
    } else {
        plhs[0] = mxCreateDoubleMatrix(nOmega, 1, mxREAL);
        /* mxCreateDoubleMatrix automatically zeros out the entries */
        output = mxGetPr(plhs[0]);
    }

    /* Check for the optional input argument that specifies the cutoff
     * for level-3 BLAS */
    BLAS_CUTOFF = 0.2;  
    MODE_THREE = false;
    if ( nrhs > 4 ) {
        BLAS_CUTOFF = (double)*mxGetPr( prhs[4] );
        MODE_THREE = true;
    } else if (nrhs > 3 ) {
        /* decide if this third argument is BLAS_CUTOFF,
         * or omegaX */
        nOmega1 = mxGetM( prhs[3] );
        nOmega2 = mxGetN( prhs[3] );
        if ( (nOmega1 == 1) && (nOmega2 == 1) ) {
            BLAS_CUTOFF = (double)*mxGetPr( prhs[3] );
            MODE_THREE = false;
        } else {
            MODE_THREE = true;
        }
    }
    /* Decide if we'll make a level-3 BLAS call */
    USE_BLAS = ( nOmega >= BLAS_CUTOFF * (M*N) );



    /* 
     * We have 3 "modes":
     * Mode 1 
     *      if omega is a vector of linear indices
     * Mode 2
     *      if omega is given only by the nonzero elements of an input 
     *      sparse matrix Y
     * Mode 3
     *      if omega is given as a set of subscripts, i.e. omegaX, omegaY
     *      then the "for" loop is slightly different
     *
     *  October 29, 2009: changing this a bit.
     *      For any of the modes, we first check if length(omega) > .8*M*N
     *      If so, we make a level-3 BLAS call (dgemm), and then use this
     *      matrix as needed.
     *      
     */


    /* determine if we're on a 64-bit processor */
    LARGE_BIT = ( sizeof( size_t ) > 4 );



    if ( USE_BLAS ) {
        /* we need to compute A itself, so use level-3 BLAS */
        transA = 'N';
        transB = 'T';
        LDA = M;
        LDB = N;

        if (COMPLEX) {
            alpha_cplx.re = 1.0;  alpha_cplx.im = 0.0;
            beta_cplx.re = 0.0;   beta_cplx.im  = 0.0;
            /* need to make a new complex data structure */
            output_cplx = (complex16 *) mxMalloc( M*N*sizeof(complex16) );
            if ( output_cplx == NULL ) {
                /* we don't have enough RAM available */
                USE_BLAS = false;
            } else {
                my_zgemm(&transA,&transB,(int*)&M,(int*)&N,(int*)&K,
                        &alpha_cplx,U_cplx,(int*)&LDA,V_cplx,(int*)&LDB,&beta_cplx,output_cplx,(int *)&M );
            }
        } else {
            outputBLAS = (double *)mxMalloc( M*N*sizeof(double) );
            if ( outputBLAS == NULL ){
                /* we don't have enough RAM available */
                USE_BLAS = false;
            } else {
                alpha = 1.0;
                beta = 0.0;
                my_dgemm(&transA,&transB,(int*)&M,(int*)&N,(int*)&K,
                        &alpha,U,(int*)&LDA,V,(int*)&LDB,&beta,outputBLAS,(int*)&M );
            }
        }
    }



    /* --- MODE 1 ---*/
    if ( (!MODE_THREE) && (!SPARSE) ) {
         /* by default, make output the same shape (i.e. row- or column-
          * vector) as the input "omega" */
        mxSetM( plhs[0], mxGetM( prhs[2] ) );
        mxSetN( plhs[0], mxGetN( prhs[2] ) );
        omega = mxGetPr(prhs[2]);

        if ( !USE_BLAS ) {
            
            if ( (COMPLEX) && (LARGE_BIT) ) {
              for ( k=0 ; k < nOmega ; k++ ){
                i = (mwIndex) ( (mwIndex)(omega[k]-1) % M);
                j = (mwIndex) ( (mwIndex)(omega[k]-1)/ M);
                
                /* The following zdot call doesn't work for me on my 64-bit linux system */
                /*  This was a typo: N, instead of i
                temp_cplx = my_zdot( (int*) &K, U_cplx+N, (int*)&M, V_cplx+j, (int*)&N ); 
                temp_cplx = my_zdot( (int*) &K, U_cplx+i, (int*)&M, V_cplx+j, (int*)&N ); 
                */

                /* So, implement the BLAS call myself
                 * ZDOTU(N,ZX,INCX,ZY,INCY)  
                 * */
                temp_cplx.re = 0.0;
                temp_cplx.im = 0.0;
                for ( m=0 ; m < K ; m++ ){
                    temp_cplx.re += U_cplx[i+m*M].re * V_cplx[j+m*N].re
                        - U_cplx[i+m*M].im * V_cplx[j+m*N].im;
                    temp_cplx.im += U_cplx[i+m*M].im * V_cplx[j+m*N].re
                        + U_cplx[i+m*M].re * V_cplx[j+m*N].im;
                }
                output[k] = temp_cplx.re;
                output_imag[k] = temp_cplx.im;
                
              }
            } else if (COMPLEX) {
              for ( k=0 ; k < nOmega ; k++ ){
                i = (mwIndex) ( (mwIndex)(omega[k]-1) % M);
                j = (mwIndex) ( (mwIndex)(omega[k]-1)/ M);
                temp_cplx = my_zdot( (int*) &K, U_cplx+i, (int*)&M, V_cplx+j, (int*)&N ); 
                output[k] = temp_cplx.re;
                output_imag[k] = temp_cplx.im;
              }

            } else {
              for ( k=0 ; k < nOmega ; k++ ){
                /* don't forget that Matlab is 1-indexed, C is 0-indexed */
                i = (mwIndex) ( (mwIndex)(omega[k]-1) % M);
                j = (mwIndex) ( (mwIndex)(omega[k]-1)/ M);
                output[k] = my_ddot( (int*)&K, U+i, (int*)&M, V+j, (int*)&N );
              }
            }
        } else {  /* (USE_BLAS) is true */
            /* We already have the full matrix, so just find the entries */
            if (COMPLEX) {
                for ( k=0 ; k < nOmega ; k++ ){
                    output[k] = output_cplx[ (mwIndex)omega[k] -1 ].re;
                    output_imag[k] = output_cplx[ (mwIndex)omega[k] - 1 ].im;
                }
                mxFree( output_cplx );
            } else {
                for ( k=0 ; k < nOmega ; k++ ){
                    /*
                    i = (mwIndex) ( (mwIndex)(omega[k]-1) % M);
                    j = (mwIndex) ( (mwIndex)(omega[k]-1)/ M);
                    output[k] = outputBLAS( j*M + i );
                    */
                    /* This now simplifies a lot! */
                    output[k] = outputBLAS[ (mwIndex)omega[k] - 1 ];
                }
                mxFree( outputBLAS );
            }
        }

    } else {

        /* ----------- MODE 2 ------------------------------- */

        if (SPARSE) {
            /* sparse array indices in Matlab are rather confusing;
             * see the mxSetJc help file to get started.  The Ir index
             * is straightforward: it contains rows indices of nonzeros,
             * in column-major order.  But the Jc index is tricky... 
             * Basically, Jc (which has N+1 entries, not nnz entries like Ir)
             * tells you which Ir entries correspond to the jth row, thus fully 
             * specifying the indices.  Ir[ Jc[j]:Jc[J+1] ] are the rows
             * that correspond to column j. This works because Ir is
             * in column-major order.   For this to work (and match A(omega)),
             * we need omega to be sorted!  
             *Note: everything is already 0-based, not 1-based
             * */
            omegaI = mxGetIr( prhs[2] );
            omegaJ = mxGetJc( prhs[2] );
            
            if ( USE_BLAS ) {
                /* We already have the full matrix, so just find the entries */
                if (COMPLEX) {
                    for ( j=0 ; j < N ; j++ ){
                        for ( k = omegaJ[j] ; k < omegaJ[j+1] ; k++ ) {
                            i = omegaI[k];
                            output[k] = output_cplx[ j*M + i ].re;
                            output_imag[k] = output_cplx[ j*M + i ].im;
                        }
                    }
                    mxFree( output_cplx );
                } else {
                    for ( j=0 ; j < N ; j++ ){
                        for ( k = omegaJ[j] ; k < omegaJ[j+1] ; k++ ) {
                            i = omegaI[k];
                            output[k] = outputBLAS[ j*M + i ];
                        }
                    }
                    mxFree( outputBLAS );
                }

            } else if ((COMPLEX)&&(LARGE_BIT)) {
                for ( j=0 ; j < N ; j++ ){
                    for ( k = omegaJ[j] ; k < omegaJ[j+1] ; k++ ) {
                        i = omegaI[k];

                        temp_cplx.re = 0.0;
                        temp_cplx.im = 0.0;
                        for ( m=0 ; m < K ; m++ ){
                            temp_cplx.re += U_cplx[i+m*M].re * V_cplx[j+m*N].re
                                - U_cplx[i+m*M].im * V_cplx[j+m*N].im;
                            temp_cplx.im += U_cplx[i+m*M].im * V_cplx[j+m*N].re
                                + U_cplx[i+m*M].re * V_cplx[j+m*N].im;
                        }
                        output[k] = temp_cplx.re;
                        output_imag[k] = temp_cplx.im;
                    }
                }
            } else if (COMPLEX) {
                for ( j=0 ; j < N ; j++ ){
                    for ( k = omegaJ[j] ; k < omegaJ[j+1] ; k++ ) {
                        i = omegaI[k];
                        temp_cplx = my_zdot((int*) &K, U_cplx+i, (int*)&M, V_cplx+j, (int*)&N );
                        output[k] = temp_cplx.re;
                        output_imag[k] = temp_cplx.im;
                    }
                }
            } else {
                /* simple case */
                for ( j=0 ; j < N ; j++ ){
                    for ( k = omegaJ[j] ; k < omegaJ[j+1] ; k++ ) {
                        i = omegaI[k];
                        output[k] = my_ddot( (int*)&K, U+i, (int*)&M, V+j, (int*)&N );
                    }
                }
            }
        /* ----------- MODE 3 ------------------------------- */
        } else {
            /* we have omegaX and omegaY, the row and column indices */
            nOmega1 = mxGetM( prhs[3] );
            nOmega2 = mxGetN( prhs[3] );
            if ( (nOmega1 != 1) && (nOmega2 != 1) ) {
                printUsage();
                mexErrMsgTxt("OmegaY must be a vector");
            }
            nOmega1 = nOmega1 < nOmega2 ? nOmega2 : nOmega1;
            if ( nOmega1 != nOmega ) {
                printUsage();
    mexErrMsgTxt("In subscript index format, subscript vectors must be same length.");
            }
            omegaX = mxGetPr(prhs[2]);
            omegaY = mxGetPr(prhs[3]);

            if ( USE_BLAS ) {
                /* We already have the full matrix, so just find the entries */
                if (COMPLEX) {
                    for ( k=0 ; k < nOmega ; k++ ){
                        i = omegaX[k] - 1;
                        j = omegaY[k] - 1;
                        output[k] = output_cplx[ j*M + i ].re;
                        output_imag[k] = output_cplx[ j*M + i ].im;
                    }
                    mxFree( output_cplx );
                } else {
                    for ( k=0 ; k < nOmega ; k++ ){
                        i = omegaX[k] - 1;
                        j = omegaY[k] - 1;
                        output[k] = outputBLAS[ j*M + i ];
                    }
                    mxFree( outputBLAS );
                }

            } else if ((COMPLEX)&&(LARGE_BIT)) {
                for ( k=0 ; k < nOmega ; k++ ){
                    i = omegaX[k] - 1;
                    j = omegaY[k] - 1;
                    temp_cplx.re = 0.0;
                    temp_cplx.im = 0.0;
                    for ( m=0 ; m < K ; m++ ){
                        temp_cplx.re += U_cplx[i+m*M].re * V_cplx[j+m*N].re
                            - U_cplx[i+m*M].im * V_cplx[j+m*N].im;
                        temp_cplx.im += U_cplx[i+m*M].im * V_cplx[j+m*N].re
                            + U_cplx[i+m*M].re * V_cplx[j+m*N].im;
                    }
                    output[k] = temp_cplx.re;
                    output_imag[k] = temp_cplx.im;
                }
            } else if (COMPLEX) {
                for ( k=0 ; k < nOmega ; k++ ){
                    i = omegaX[k] - 1;
                    j = omegaY[k] - 1;
                    temp_cplx = my_zdot( (int*)&K, U_cplx+i, (int*)&M, V_cplx+j, (int*)&N );
                    output[k] = temp_cplx.re;
                    output_imag[k] = temp_cplx.im;
                }
            } else {
                for ( k=0 ; k < nOmega ; k++ ){
                    i = omegaX[k] - 1;
                    j = omegaY[k] - 1;
                    output[k] = my_ddot( (int*)&K, U+i, (int*)&M, V+j, (int*)&N );
                }
            }

        }
    }

}
