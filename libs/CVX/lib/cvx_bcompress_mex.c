#include "mex.h"
#include <stddef.h>

#if !defined(HAVE_OCTAVE) && ( !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 ) )
typedef int mwIndex;
typedef int mwSize;
#endif

/*
% Copyright 2005-2016 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
*/

#define FULL_COMPRESS    0
#define MAGNITUDE_ONLY   1
#define NO_NORMALIZATION 2

#define MY_EPS 1.776356839400251e-015
mwIndex *ss_ir, *ss_jc;
mwSize *ss_temp;
int ss_mode;
double *ss_pr;

static int mycomp( mwSize c1, mwSize c2 )
{
    mwIndex d1, d2;
    mwIndex *j1, *j2, *rb1, *rb2, *re1, *re2;
    double nm1, nm2;
    double *pb1, *pb2, *pe1;

    /*
     * Comparison #1: number of elements
     */

    j1 = ss_jc + c1;
    j2 = ss_jc + c2;
    d1 = j1[1] - j1[0];
    d2 = j2[1] - j2[0];
    if ( d1 < d2 ) return -1;
    if ( d1 > d2 ) return +1;
    if ( d1 == 0 ) return  0;
    
    /*
     * Comparison #2: Lexical pattern
     */

    rb1 = ss_ir + j1[0];
    re1 = rb1   + d1;
    rb2 = ss_ir + j2[0];
    re2 = rb2   + d2;
    while ( 1 ) {
        if ( *rb1 < *rb2  ) return -1;
        if ( *rb1 > *rb2  ) return +1;
        if ( ++rb1 == re1 )
            if ( ++rb2 == re2 ) break;
            else return -1;
        else if ( ++rb2 == re2 )
            return +1;
    }
    
    /*
     * Comparison #3: sort by normalized column values
     */

    pb1 = ss_pr + j1[0];
    pb2 = ss_pr + j2[0];
    switch ( ss_mode ) {
        case MAGNITUDE_ONLY:
            if ( *pb1 < 0 ) {
                if ( *pb2 > 0 )
                    return -1;
            } else if ( *pb2 < 0 )
                return +1;
        default:
        case FULL_COMPRESS:
            if ( d1 == 1 )
                return 0;
            break;
        case NO_NORMALIZATION:
            pe1 = pb1 + d1;
            while ( 1 ) {
                if ( *pb1 < *pb2  ) return -1;
                if ( *pb1 > *pb2  ) return +1;
                if ( ++pb1 == pe1 ) return 0;
                ++pb2;
            }
            break;
    }
    pe1 = pb1 + d1;
    nm1 = *pb2++;
    nm2 = *pb1++;
    while ( 1 ) {
        double b1 = nm1 * *pb1,
               b2 = nm2 * *pb2,
               bd = b1 - b2,
               bs = MY_EPS * ( ( b1 < 0 ? -b1 : +b1 ) + ( b2 < 0 ? -b2 : +b2 ) );
        if ( bd > bs || bd < -bs ) 
            return b1 < b2 ? -1 : +1;
        if ( ++pb1 == pe1 ) 
            return 0;
        ++pb2;
    }
    return 0;
}

void merge( mwSize *nb, mwSize *nm, mwSize *ne )
{
    mwSize *i = nb, *j = nm, *t = ss_temp;
    for ( ;; )
        if ( mycomp( *i, *j ) <= 0 ) {
            *t++ = *i++;
            if ( i == nm ) {
                while ( j != ne )
                    *t++ = *j++;
                break;
            }
        } else {
            *t++ = *j++;
            if ( j == ne ) {
                while ( i != nm )
                    *t++ = *i++;
                break;
            }
        }
    for ( i = nb, t = ss_temp ; i != ne ; *i++ = *t++ );
}

static void merge_sort( mwSize *nb, mwSize* ne )
{
    mwSize *nm;
    ptrdiff_t numel = ne - nb;
    if ( numel < 2 ) return;
    nm = nb + ( numel / 2 );
    merge_sort( nb, nm );
    merge_sort( nm, ne );
    merge( nb, nm, ne );
}

void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )

{
    mwSize n, k, col, nskip, lastcol, *ss_ndxs;
    double *map, *scl, lastnorm;
    int first = 1;
    
    n       = mxGetN(  prhs[0] );
    ss_ir   = mxGetIr( prhs[0] );
    ss_jc   = mxGetJc( prhs[0] );
    ss_pr   = mxGetPr( prhs[0] );
    ss_mode = nrhs < 2 ? FULL_COMPRESS : (int)mxGetScalar( prhs[1] );
    nskip   = nrhs < 3 ? 0 : (int)mxGetScalar(prhs[2]);
    ss_ndxs = mxCalloc( n, sizeof(mwSize) );
    ss_temp = mxCalloc( n, sizeof(mwSize) );
    
    /*
     * Sort the column indices
     */

    for ( col = 0 ; col != n ; ++col )
        ss_ndxs[col] = col;
    if ( nskip < n )
        merge_sort( ss_ndxs + nskip, ss_ndxs + n );

    /*
     * Determine which rows are unique, and which are scales of another row
     */

    plhs[0] = mxCreateDoubleMatrix( (mwSize)1, n, mxREAL );
    plhs[1] = mxCreateDoubleMatrix( (mwSize)1, n, mxREAL );
    if ( plhs[0] == 0 || plhs[1] == 0 )
        mexErrMsgTxt( "Unable to allocate output arguments" );
    map = mxGetPr( plhs[0] );
    scl = mxGetPr( plhs[1] );
    lastnorm = 0.0;
    for ( k = 0 ; k < n ; ++k ) {
        col = ss_ndxs[k];
        if ( col < nskip ) {
            map[col] = col + 1;
            scl[col] = 1.0;
        } else if ( ss_jc[col] == ss_jc[col+1] ) {
            map[col] = col + 1;
            scl[col] = 0;
        } else if ( first || mycomp( lastcol, col ) != 0 ) {
            lastcol  = col;
            lastnorm = ss_pr[ss_jc[col]];
            map[col] = lastcol + 1;
            scl[col] = 1.0;
            first = 0;
        } else {
            map[col] = lastcol + 1;
            scl[col] = ss_pr[ss_jc[col]] / lastnorm;
        }
    }
}
