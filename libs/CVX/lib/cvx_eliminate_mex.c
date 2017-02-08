#include "mex.h"
#include <stddef.h>
#include <math.h>

#if !defined(HAVE_OCTAVE) && ( !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 ) )
typedef int mwIndex;
typedef int mwSize;
#endif

static const double MAX_GROWTH = 16;

/*
% Copyright 2005-2016 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
*/

typedef struct {
    mwSize *values;
    mwSize *temp;
} col_sort_struct;

static int mycomp( mwIndex i, mwIndex j, col_sort_struct* ss )
{
    mwSize vi = ss->values[i], vj = ss->values[j];
    return vi > vj ? -1 : ( vi < vj ? +1 : ( i < j ? -1 : +1 ) );
}

static void merge( mwSize *nb, mwSize* nm, mwSize* ne, col_sort_struct* ss )
{
    mwSize *i = nb, *j = nm, *t = ss->temp;
    for ( ;; )
        if ( mycomp( *i, *j, ss ) <= 0 ) {
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
    for ( i = nb, t = ss->temp ; i != ne ; *i++ = *t++ );
}

void merge_sort( mwSize *nb, mwSize* ne, col_sort_struct* ss )
{
    mwSize *nm;
    ptrdiff_t numel = ne - nb;
    if ( numel < 2 ) return;
    nm = nb + ( numel / 2 );
    merge_sort( nb, nm, ss );
    merge_sort( nm, ne, ss );
    merge( nb, nm, ne, ss );
}

void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )

{
    mwIndex *row_index_A = mxGetIr( prhs[0] ),
            *col_count_A = mxGetJc( prhs[0] );
    mwSize  m            =  mxGetM( prhs[0] ),
            n            =  mxGetN( prhs[0] ),
            nobj         = (mwSize)*mxGetPr( prhs[1] ),
            nnzA         = (mwSize)col_count_A[n],
            nnzA1        = nnzA + 1,
            j, nQ, nnzQ,
            *row_counts  = (mwSize*)mxCalloc( m + n, sizeof(mwSize) ),
            *candidates  = row_counts + m;
    double  *value_A     = mxGetPr( prhs[0] ),
            *row_flags   = mxGetPr( plhs[0] = mxCreateNumericMatrix( m, (mwSize)1, mxDOUBLE_CLASS, mxREAL ) ),
            *col_flags   = mxGetPr( plhs[1] = mxCreateNumericMatrix( (mwSize)1, n, mxDOUBLE_CLASS, mxREAL ) ),
            *reserved    = mxGetPr( prhs[2] ),
            *c_reserved  = mxGetPr( prhs[3] );
    
    /*
     * This loop has two functions:
     * --- Locate one-variable equality constraints: xi == bj, and target 
     *     them for removal if the variable is free
     * --- Count the number of elements in each row
     */
    
    nQ = nnzQ = 0;
    reserved[0] = 1;
    row_flags[0] = -1;
    for ( j = nobj ; j != n ; ++j ) {
        mwIndex kBeg = col_count_A[j],
                kEnd = col_count_A[j+1], k;
        for ( k = kBeg ; k != kEnd ; ++k )
            ++row_counts[row_index_A[k]];
        if ( c_reserved[j] == 0 && kBeg != kEnd ) {
            mwIndex rb = row_index_A[kBeg];
            if ( kBeg - kEnd == ( rb == 0 ? 2 : 1 ) ) {
                rb = row_index_A[kEnd-1];
                if ( reserved[rb] == 0 ) {
                    candidates[nQ++] = j;
                    reserved[rb] = 1;
                    c_reserved[j] = 1;
                    col_flags[j] = rb + 1;
                    row_flags[rb] = 1;
                    ++nnzQ;
                }
            }
        }
    }
    for ( j = nobj ; j != n ; ++j ) {
        if ( c_reserved[j] == 0 ) {
            mwIndex kBeg  = col_count_A[j],
                    kEnd  = col_count_A[j+1], k, 
                    rb    = 0;
            double  c_max = 0;
            for ( k = kBeg ; k != kEnd ; ++k ) {
                mwIndex tr = row_index_A[k];
                if ( tr != 0 ) {
                    if ( row_flags[tr] && reserved[tr] ) {
                        c_max = 0;
                        break;
                    } else {
                        double temp = fabs(value_A[k]);
                        if ( c_max < temp ) 
                            c_max = temp;
                    }
                }
            }
            if ( c_max != 0 ) {
                mwIndex nc = kEnd - kBeg;
                c_max /= MAX_GROWTH;
                for ( k = kBeg ; k != kEnd ; ++k ) {
                    mwIndex tr = row_index_A[k];
                    if ( tr == 0 ) {
                        --nc;
                    } else if ( reserved[tr] == 0 ) {
                        mwSize nr = row_counts[tr],
                               n1 = ( nr - 1 ) * ( nc - 1 ),
                               n2 = nr + nc - 1;
                        if ( n1 <= n2 && row_flags[tr] <= 0 ) {
                            double temp = fabs(value_A[k]);
                            if ( temp > c_max || ( temp == c_max && row_flags[rb] && !row_flags[tr] ) ) {
                                rb = tr;
                                c_max = temp;
                            }
                        }
                    } else if ( row_flags[tr] ) {
                        rb = 0;
                        break;
                    }
                }
            }
            if ( rb != 0 ) {
                for ( k = kBeg ; k != kEnd ; ++k ) {
                    mwIndex tr = row_index_A[k];
                    if ( tr != 0 ) {
                        double fg = row_flags[tr];
                        if ( fg > 0 ) {
                            ++fg;
                            ++nnzQ;
                        } else
                            --fg;
                        row_flags[tr] = fg;
                    }
                }
                col_flags[j] = rb + 1;
                nnzQ += ( row_flags[rb] = - row_flags[rb] );
                candidates[nQ++] = j;
            }
        }
    }
    for ( j = 0 ; j != m ; ++j )
        if ( row_flags[j] < 0 )
            row_flags[j] = 0;
    if ( nnzQ > nQ ) {
        col_sort_struct ss;
        mwIndex i, iEnd;
        ss.values = (mwSize*)mxCalloc( 2 * n, sizeof(mwSize) );
        ss.temp   = ss.values + n;
        nnzQ = 0;
        for ( i = iEnd = 0 ; i != nQ ; ++i ) {
            mwIndex j    = candidates[i],
                    r    = (mwIndex)col_flags[j] - 1,
                    kBeg = col_count_A[j],
                    kEnd = col_count_A[j+1], k;
            ss.values[j] = row_flags[r] - 2;
            for ( k = kBeg ; k != kEnd ; ++k )
                if ( row_flags[row_index_A[k]] != 0 )
                    ++ss.values[j];
            if ( ss.values[j] != 0 ) {
                candidates[iEnd++] = j;
                nnzQ += ss.values[j];
            }
        }
        nnzQ /= 2;
        merge_sort( candidates, candidates + iEnd, &ss );
        for ( i = 0 ; i != iEnd && nnzQ > 0 ; ++i ) {
            mwSize   j = candidates[i];
            mwIndex  r    = col_flags[j] - 1,
                     kBeg = col_count_A[j],
                     kEnd = col_count_A[j+1], k;
            mwSize  found = 0;
            for ( k = kBeg ; k != kEnd ; ++k ) {
                mwIndex tr = row_index_A[k];
                if ( r != tr && row_flags[tr] != 0 ) {
                    --row_flags[tr];
                    ++found;
                }
            }
            if ( found ) {
                nnzQ -= found + row_flags[r] - 1;
                row_flags[r] = col_flags[j] = 0;
            }
        }
    }
}
