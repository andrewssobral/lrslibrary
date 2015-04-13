/* 
 * Stephen Becker, 10/30/09
 * Computes A(omega), where A = U*V' is never explicitly computed 
 * srbecker@caltech.edu
 *
 * This is a compantion to XonOmega.c
 * The difference is that this version expects U' and V' to be passed in,
 * not U and V.  This version should be faster
 *
 * See XonOmega.c for details on usage and compilation
 *
 * */

#include "mex.h"
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
    mexPrintf("XonOmegaTranspose.c: usage is\n\t b = XonOmegaTranspose(U.',V.',omega)\n");
    mexPrintf(" or b = XonOmegaTranspose(U.',V.',omega,cutoff)\n");
    mexPrintf("where A = U*V' and b = A(omega)\n");
    mexPrintf("\nAlternative usage is\n\t b = XonOmegaTranspose(U.',V.',omegaI,omegaJ),\n where [omegaI,omegaJ] = ind2sub(...,omega)\n");
    mexPrintf("\nor b = XonOmegaTranspose(U',V',omegaI,omegaJ,cutoff)\n");
    mexPrintf("\nAnother alternative usage is\n\t b = XonOmegaTranspose(U.',V.',OMEGA),\n where OMEGA is a sparse matrix with nonzeros on omega.\nThis will agree with the other forms of the command if omega is sorted\n");
    mexPrintf("Can also call this as b = XonOmegaTranspose(U.',V.',OMEGA, cutoff)\n\n");
    mexPrintf("'cutoff' affects how the computation is done\n");
    mexPrintf("if length(omega) > cutoff*numel(A), then the entire A matrix is first formed\n");
    mexPrintf("otherwise, A is never formed explicitly\n");
    mexPrintf("If U and V are complex, then make sure A=U*V' and not A=U*V.'\n");
    mexPrintf("If U and V are complex, make sure to pass in the variables U.' and V.',not U' and V''\n\n");
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

    int one = 1;

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
    /* TRANSPOSE VERSION: switch mxGetM and mxGetN */
    M  = mxGetN(prhs[0]);
    K  = mxGetM(prhs[0]);
    N  = mxGetN(prhs[1]);
    K2  = mxGetM(prhs[1]);
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
    BLAS_CUTOFF = 0.4;  
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
        /* TRANSPOSE VERSION: switch N and T below, and change the step size and LDA */
        /*
        transA = 'N';
        transB = 'T';
        LDA = M;
        LDB = N; */
        transA = 'T';
        transB = 'N';
        LDA = K;
        LDB = K;

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
            if ( outputBLAS == NULL ) {
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
              /* TRANSPOSE VERSION: this needs to be tested!!! */
              for ( k=0 ; k < nOmega ; k++ ){
                i = (mwIndex) ( (mwIndex)(omega[k]-1) % M);
                j = (mwIndex) ( (mwIndex)(omega[k]-1)/ M);
                
                /* implement the BLAS call myself, since BLAS isn't working for me
                 * ZDOTU(N,ZX,INCX,ZY,INCY)  
                 * */
                temp_cplx.re = 0.0;
                temp_cplx.im = 0.0;
                for ( m=0 ; m < K ; m++ ){
                    /*temp_cplx.re += U_cplx[i+m*M].re * V_cplx[j+m*N].re
                        - U_cplx[i+m*M].im * V_cplx[j+m*N].im;
                    temp_cplx.im += U_cplx[i+m*M].im * V_cplx[j+m*N].re
                        + U_cplx[i+m*M].re * V_cplx[j+m*N].im; */
                    temp_cplx.re += U_cplx[K*i+m].re * V_cplx[K*j+m].re
                        - U_cplx[K*i+m].im * V_cplx[K*j+m].im;
                    temp_cplx.im += U_cplx[K*i+m].im * V_cplx[K*j+m].re
                        + U_cplx[K*i+m].re * V_cplx[K*j+m].im;
                }
                output[k] = temp_cplx.re;
                output_imag[k] = temp_cplx.im;
                
              }
            } else if (COMPLEX) {
              /* TRANSPOSE VERSION: UPDATED */
              for ( k=0 ; k < nOmega ; k++ ){
                i = (mwIndex) ( (mwIndex)(omega[k]-1) % M);
                j = (mwIndex) ( (mwIndex)(omega[k]-1)/ M);
                /*temp_cplx = my_zdot( (int*) &K, U_cplx+i, (int*)&M, V_cplx+j, (int*)&N ); */
                temp_cplx = my_zdot( (int*) &K, U_cplx+i*K, (int*)&one, V_cplx+j*K, (int*)&one ); 
                output[k] = temp_cplx.re;
                output_imag[k] = temp_cplx.im;
              }

            } else {
              /* TRANSPOSE VERSION: UPDATED */
              for ( k=0 ; k < nOmega ; k++ ){
                /* don't forget that Matlab is 1-indexed, C is 0-indexed */
                i = (mwIndex) ( (mwIndex)(omega[k]-1) % M);
                j = (mwIndex) ( (mwIndex)(omega[k]-1)/ M);
                /*output[k] = my_ddot( (int*)&K, U+i, (int*)&M, V+j, (int*)&N );*/
                output[k] = my_ddot( (int*)&K, U+i*K, (int*)&one, V+j*K, (int*)&one );
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
                /* TRANSPOSE VERSION: this needs to be tested!!! */
                for ( j=0 ; j < N ; j++ ){
                    for ( k = omegaJ[j] ; k < omegaJ[j+1] ; k++ ) {
                        i = omegaI[k];

                        temp_cplx.re = 0.0;
                        temp_cplx.im = 0.0;
                        for ( m=0 ; m < K ; m++ ){
                            /*temp_cplx.re += U_cplx[i+m*M].re * V_cplx[j+m*N].re
                                - U_cplx[i+m*M].im * V_cplx[j+m*N].im;
                            temp_cplx.im += U_cplx[i+m*M].im * V_cplx[j+m*N].re
                                + U_cplx[i+m*M].re * V_cplx[j+m*N].im;*/
                            temp_cplx.re += U_cplx[K*i+m].re * V_cplx[K*j+m].re
                                - U_cplx[K*i+m].im * V_cplx[K*j+m].im;
                            temp_cplx.im += U_cplx[K*i+m].im * V_cplx[K*j+m].re
                                + U_cplx[K*i+m].re * V_cplx[K*j+m].im;
                        }
                        output[k] = temp_cplx.re;
                        output_imag[k] = temp_cplx.im;
                    }
                }
            } else if (COMPLEX) {
                /* TRANSPOSE VERSION: UPDATED */
                for ( j=0 ; j < N ; j++ ){
                    for ( k = omegaJ[j] ; k < omegaJ[j+1] ; k++ ) {
                        i = omegaI[k];
                        /* temp_cplx = my_zdot((int*) &K, U_cplx+i, (int*)&M, V_cplx+j, (int*)&N ); */
                        temp_cplx = my_zdot((int*) &K, U_cplx+i*K, (int*)&one, V_cplx+j*K, (int*)&one );
                        output[k] = temp_cplx.re;
                        output_imag[k] = temp_cplx.im;
                    }
                }
            } else {
                /* simple case */
                /* TRANSPOSE VERSION: UPDATED */
                for ( j=0 ; j < N ; j++ ){
                    for ( k = omegaJ[j] ; k < omegaJ[j+1] ; k++ ) {
                        i = omegaI[k];
                        /*output[k] = my_ddot( (int*)&K, U+i, (int*)&M, V+j, (int*)&N );*/
                        output[k] = my_ddot( (int*)&K, U+i*K, (int*)&one, V+j*K, (int*)&one );
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
                /* TRANSPOSE VERSION: this needs to be tested!!! */
                for ( k=0 ; k < nOmega ; k++ ){
                    i = omegaX[k] - 1;
                    j = omegaY[k] - 1;
                    temp_cplx.re = 0.0;
                    temp_cplx.im = 0.0;
                    for ( m=0 ; m < K ; m++ ){
                        /*temp_cplx.re += U_cplx[i+m*M].re * V_cplx[j+m*N].re
                            - U_cplx[i+m*M].im * V_cplx[j+m*N].im;
                        temp_cplx.im += U_cplx[i+m*M].im * V_cplx[j+m*N].re
                            + U_cplx[i+m*M].re * V_cplx[j+m*N].im; */
                        temp_cplx.re += U_cplx[K*i+m].re * V_cplx[K*j+m].re
                            - U_cplx[K*i+m].im * V_cplx[K*j+m].im;
                        temp_cplx.im += U_cplx[K*i+m].im * V_cplx[K*j+m].re
                            + U_cplx[K*i+m].re * V_cplx[K*j+m].im;
                    }
                    output[k] = temp_cplx.re;
                    output_imag[k] = temp_cplx.im;
                }
            } else if (COMPLEX) {
                /* TRANSPOSE VERSION: UPDATED */
                for ( k=0 ; k < nOmega ; k++ ){
                    i = omegaX[k] - 1;
                    j = omegaY[k] - 1;
                    /* temp_cplx = my_zdot( (int*)&K, U_cplx+i, (int*)&M, V_cplx+j, (int*)&N ); */
                    temp_cplx = my_zdot( (int*)&K, U_cplx+i*K, (int*)&one, V_cplx+j*K, (int*)&one );
                    output[k] = temp_cplx.re;
                    output_imag[k] = temp_cplx.im;
                }
            } else {
                /* TRANSPOSE VERSION: UPDATED */
                for ( k=0 ; k < nOmega ; k++ ){
                    i = omegaX[k] - 1;
                    j = omegaY[k] - 1;
                    /* output[k] = my_ddot( (int*)&K, U+i, (int*)&M, V+j, (int*)&N ); */
                    output[k] = my_ddot( (int*)&K, U+i*K, (int*)&one, V+j*K, (int*)&one );
                }
            }

        }
    }

}
