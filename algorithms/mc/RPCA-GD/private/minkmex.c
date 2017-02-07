/**************************************************************************
 * function [res loc] = minkmex(list, k)
 * Matlab C-Mex
 * Purpose: Same as MINK, i.e.,
 * Return in RES the K smallest elements of 2D matrix
 * LOC is Location of the smallest values
 *  - For full matrix, LOC contains linear indexing of the matrix.
 *              RES == LIST(LOC)
 *  - For sparse, location is returned in subindexes form by calling 
 *              [RES I J] = minkmex(list, k)
 *              RES == getsparse(list,I,J)
 * This MEX works on double array only, and output RES is unsorted column
 * Complexity O(n) where n is size of list
 * Note: Implementation of type "non-destructive", i.e., the original data 
 * will not be effectively swapped, but we keep track of a table of
 * permutation indexes.
 * Algorithm according to http://en.wikipedia.org/wiki/Selection_algorithm
 * Compilation: mex -O -v minkmex.c
 * Author Bruno Luong <brunoluong@yahoo.com>
 * Last update: 10-Jan-2010
 *              16-August-2010: change type to mwSignedIndex and
 *                              check nansplit>=0
 *              27-Aug-2011: correct bug for sparse array
 *				26-Apr-2013: fix memset warning and remove C++ comment style
 *************************************************************************/

#include "mex.h"
#include "matrix.h"

/* Define correct type depending on platform */
#if defined(_MSC_VER) || defined(__BORLANDC__)
typedef unsigned __int64 ulong64;
#elif defined(_LCC)
typedef long long long64;
typedef unsigned long long ulong64;
#else
typedef unsigned long long ulong64;
#endif

/* Global variables, used to avoid stacking them during recusive call since
   they do not change */
mwIndex k;
mwIndex *pos;
double *list;

#define MIDPOINT 0
#define MEDIAN3 1
#define MEDIANMEDIANS 2

/* Pivot Strategy, use one of the above */
#define PIVOT MIDPOINT

/*************************************************************************/
/*Find the index of the Median of the elements
of array that occur at every "shift" positions.*/
mwSize findMedianIndex(mwSize left, mwSize right, mwSize shift)
{
    mwSize tmp, groups, k;
    double minValue;
    mwSize *pi, *pj, *pk, *pright, *pminIndex;
    
    groups = (right-left)/shift + 1;
    pk = pos + (k = left + (groups/2)*shift);
    pright = pos + right;
    for (pi=pos+left; pi<=pk; pi+= shift)
    {
        pminIndex = pi;
        minValue = list[*pminIndex];
        
        for (pj=pi; pj<=pright; pj+=shift)
            if (list[*pj]<minValue) /* Comparison */
                minValue = list[*(pminIndex=pj)];
        /* Swap pos[i] with pos[minIndex] */
        tmp = *pi;
        *pi = *pminIndex;
        *pminIndex = tmp;
    }
    
    return k;
    
} /* findMedianIndex */

/*Computes the median of each group of 5 elements and stores
  it as the first element of the group (left). Recursively does this
  till there is only one group and hence only one Median */
mwSize findMedianOfMedians(mwSize left, mwSize right)
{
    mwSize i, shift, step, tmp;
    mwSize endIndex, medianIndex;
           
    if (left==right) return left;
   
    shift = 1;
    while (shift <= (right-left))
    {
        step=shift*5;
        for (i=left; i<=right; i+=step)
        {
            if ((endIndex=i+step-1)>=right)
                endIndex=right;
            medianIndex = findMedianIndex(i, endIndex, shift);
            /* Swap pos[i] with pos[medianIndex] */
            tmp = pos[i];
            pos[i] = pos[medianIndex];
            pos[medianIndex] = tmp;
        }
        shift = step;
    }
    return left;
} /* findMedianOfMedians */

/*************************************************************************/
/*Computes the median of three points (left,right,and mid) */
mwIndex findMedianThree(mwIndex left, mwIndex right)
{
    double vleft, vright, vmid;
    mwIndex mid;
    
    if (left==right) return left;
    
    vleft = list[pos[left]];
    vright = list[pos[right]];
    vmid = list[pos[mid = (left+right+1)/2]];
    
    if (vleft<vright)
    {
        if (vmid>vright)
            return right;
        else if (vmid<vleft)
            return left;
        else
            return mid;
        
    } else { /* (vleft>=vright) */
        
        if (vmid>vleft)
            return left;
        else if (vmid<vright)
            return right;
        else
            return mid;
        
    }       
} /* findMedianThree */

/* A quiet NaN is represented by any bit pattern 
   between X'7FF80000 00000000' and X'7FFFFFFF FFFFFFFF' or 
   between X'FFF80000 00000000' and X'FFFFFFFF FFFFFFFF'. */
#define NANmask 0x7ff8000000000000
#define ISNAN(x) ((*(ulong64*)(&x) & NANmask)  == NANmask)
#define MINF 0xfff0000000000000

/* Partitioning the list around the pivot NaN.
   After runing, at exit we obtain pindex satisfied: 
   l[left]...l[index] are regular numbers (might include Inf)
   l[index+1] ... l[right] are NaN
   where l[i] := list[pos[i]] for all i */                
mwIndex partNaN(mwIndex left, mwIndex right) {
    
    mwIndex *pleft, *pright, tmp;
    mwIndex *pfirst;
    
    pfirst = pleft = pos+left;
    pright = pos+right;
    
    for (;;) {
        while ((pleft<pright) && !ISNAN(list[*pleft]))
            pleft++;
        while ((pleft<pright) && ISNAN(list[*pright]))
            pright--;
        if (pleft<pright) {
            /* Swap left and right */
            tmp = *pleft;
            *pleft = *pright;
            *pright = tmp;
            pleft++, pright--;
        }
        else {
            if (pright>=pfirst && ISNAN(list[*pright]))
                pright--;
            return (pright-pos);
        }   
    } /* for-loop */
} /* partNaN */

/*************************************************************************/
/* Partitioning the list around pivot pivotValue := l[pivotIndex];
   After runing, at exit we obtain: 
   l[left]...l[index-1] < pivotValue <= l[index] ... l[right]
   where l[i] := list[pos[i]] for all i */
mwSize partition(mwSize left, mwSize right, mwSize pivotIndex) {
    
    double pivotValue;
    mwSize *pindex, *pi, *pright;
    mwSize tmp;
    
    pright=pos+right;
    pindex=pos+pivotIndex;
    pivotValue = list[tmp = *pindex];
    /* Swap pos[pivotIndex] with pos[right] */
    *pindex = *pright;
    *pright = tmp;
      
    pindex=pos+left;
    for (pi=pindex; pi<pright; pi++)
        /* Compare with pivotValue */
        if (list[*pi] < pivotValue) {
             /* if smaller; Swap pos[index] with pos[i] */
            tmp = *pindex;
            *pindex = *pi;
            *pi = tmp;           
            pindex++;
        }

     /* Swap pos[index] with pos[right] */
    tmp = *pindex;
    *pindex = *pright;
    *pright = tmp;  
    
    return (pindex-pos); /* Pointer arithmetic */
} /* Partition */

/* Partitioning the list around pivot 0;
 * After runing, at exit we obtain: 
   l[left]...l[index-1] < 0 <= l[index] ... l[right]
   where l[i] := list[pos[i]] for all i 
   Note: at return, index might be larger than right (if all elements are
         strictly greater than zero) */
mwIndex part0(mwIndex left, mwIndex right) {
    
    mwIndex *pindex, *pi, *pright;
    mwIndex tmp;
    
    pright=pos+right;   
    pindex=pos+left;
    for (pi=pindex; pi<=pright; pi++)
        /* Compare with pivotValue of zero */
        if (list[*pi] < 0.0) { /* Compare */
             /* if larger; Swap pos[index] with pos[i] */
            tmp = *pi;
            *pi = *pindex;
            *(pindex++) = tmp;
        }

    return (pindex-pos); /* Pointer arithmetic */
} /* part0 */

/* Recursive engine (partial quicksort) */
void findFirstK(mwIndex left, mwIndex right) {
    
    mwIndex pivotIndex;

    if (right > left) {
        
#if (PIVOT==MEDIANMEDIANS)
        pivotIndex = findMedianOfMedians(left, right);
#elif (PIVOT==MEDIAN3)
        pivotIndex = findMedianThree(left, right);
#else /* MIDPOINT */
        pivotIndex = (left+right+1)/2;
#endif
        
        pivotIndex = partition(left, right, pivotIndex);
        if (pivotIndex > k)
            findFirstK(left, pivotIndex-1);
        else if (pivotIndex < k)
            findFirstK(pivotIndex+1, right);
    }
    
    return;
} /* findFirstK */

/* Create the result contains k smallest values */
mxArray* MinMaxResult(mwIndex k, mwIndex p0, mwIndex nz,
                      mwIndex kout)
{
    mwIndex i;
    mwSize dims[2];
    mxArray* Result;
    double *data;
    
    /* Create the Matrix result (first output) */
    dims[0] = kout; dims[1] = 1;
    Result = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    if (Result == NULL)
        mexErrMsgTxt("Out of memory.");
    data = mxGetPr(Result);
    /* copy negative part (p0) */
    for (i=0; i<p0; i++) data[i]=list[pos[i]];
    
    if (nz>kout-p0)
        nz = kout-p0;
    /* Fill nz zeros */ 
    memset((void*)(data+p0), 0, sizeof(double)*nz);
    
    /* copy positive part (kout - (p0+nz)) */
    for (i=p0+nz; i<kout; i++) data[i]=list[pos[i-nz]];
    
    return Result;
} /* MinMaxResult */

/* Create the result contains the locatio of k smallest values */
mxArray* LocResult(mwIndex k, mwIndex p0, mwIndex nz,
                   mwIndex kout)
{
    mwIndex i;
    mwSize dims[2];
    mxArray* Result;
    double *data;
    
    dims[0] = kout; dims[1] = 1;
    Result = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    if (Result == NULL)
        mexErrMsgTxt("Out of memory.");
    data = mxGetPr(Result);
    
    /* index of negative part */
    for (i=0; i<p0; i++) data[i]=(double)(pos[i]+1); /* one-based indexing */
    
    if (nz>kout-p0)
        nz = kout-p0;
    /* Fill nz zeros */ 
    memset((void*)(data+p0), 0, sizeof(double)*nz);
    
    /* index of positive part */
    for (i=p0+nz; i<kout; i++) data[i]=(double)(pos[i-nz]+1);
    
    return Result;
} /* LocResult */

/* FindSPzero, find the location of zeros in sparse matrix */
void FindSPzero(const mxArray* S, mwIndex nz, double* I, double *J) {
    
    mwIndex nnz, i, ib;
    mwSize m, n, ai, aj, bi, bj;    
    mwIndex *irs, *jcs;
#if defined(_LCC)
    long64 nzS; /* For some reason LCC does not accept uint64 */
#else
    ulong64 nzS;
#endif
    
    /* Get size */
    m = mxGetM(S);
    n = mxGetN(S);  
    /* Number of non-zero of S */
    nnz = *(mxGetJc(S) + n);
    
    /* Number of zeros of S */
#if defined(_LCC) /* For some reason LCC does not accept uint64 */
    nzS = (long64)m*(long64)n;
    nzS = nzS - (long64)nnz;
#else
    nzS = (ulong64)m * (ulong64)n - (ulong64)nnz;
#endif

    /* Clip nz */
    if (nz>nzS) nz=(mwIndex)nzS;
    
    /* Get the sparse index pointers */
    irs = mxGetIr(S);
    jcs = mxGetJc(S);   
        
    /* i is index of I J */
    i = 0;
    /* (ai,aj) current subindex of zero element */
    ai = aj = 0;

    /* (bi,bj) current subindex of nonzero element */
    if ((ib=0)<nnz)
    {
        bi = irs[ib];
        bj = 0;
        while (jcs[bj+1]<=ib) bj++;
    }
    else bj = n;
    
    /* Loop until all output are filled */
    while (i<nz) {
        if ((aj<bj) || (aj==bj && ai<bi)) { /* (a < b) */
            I[i] = (double)(ai+1); /* Matlab index is one-based */
            J[i] = (double)(aj+1);
            i++;
            if (++ai==m) (ai = 0, ++aj); /* increment a and wraparound */
        }
        else {
            if ((aj==bj) && (ai==bi)) /* (a == b) */
                if (++ai==m) (ai = 0, ++aj); /* increment a and wraparound */
            if (++ib<nnz) /* increment b */
            {
                bi = irs[ib];
                while (jcs[bj+1]<=ib) bj++;
            }
            else bj = n;
        }
    } /* while loop */
    
    return;
} /* FINDSPZEROS */

/* Create the result contains the location of k smallest values
 for sparse matrix */
void SpLocResult(mwIndex k, mwIndex p0, mwIndex nz, mwIndex kout, 
                 const mxArray* S, mxArray** I, mxArray** J)
{
    mwIndex i, j, pi;
    mwSize dims[2], rows, columns;
    double *dataI, *dataJ;
    mwIndex *irs, *jcs;
    
    /* Get pointers and size */
    irs = mxGetIr(S);
    jcs = mxGetJc(S);
    rows = mxGetM(S);
    columns = mxGetN(S);
    
    /* Create array of rows and column index */
    dims[0] = kout; dims[1] = 1;
    
    *I = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    if (I == NULL)
        mexErrMsgTxt("Out of memory.");
    dataI = mxGetPr(*I);
    
    *J = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    if (J == NULL)
        mexErrMsgTxt("Out of memory.");
    dataJ = mxGetPr(*J);    
    
    /* index of negative part */
    for (i=0; i<p0; i++)
    {   
        pi = pos[i];
        /* Look for the column */
        j=0;
        while (jcs[j+1]<=pi) j++;

        dataI[i]=(double)(irs[pi]+1); /* one-based indexing */
        dataJ[i]=(double)(j+1); /* one-based indexing */
    }
    
    /* Clip the number of zeros */
    if (nz>kout-p0)
        nz = kout-p0;
    /* Find the place where zeros are hidden */
    FindSPzero(S, nz, dataI+p0, dataJ+p0);
    
    /* index of positive part */
    for (i=p0+nz; i<kout; i++)
    {
        pi = pos[i-nz];
        /* Look for the column */
        j=0;
        while (jcs[j+1]<=pi) j++;

        dataI[i]=(double)(irs[pi]+1); /* one-based indexing */
        dataJ[i]=(double)(j+1); /* one-based indexing */        
    }
    
    return;
} /* SpLocResult */

/* Gateway of minkmex */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    mwSignedIndex l, i, kout, nelem, p0, nz, nansplit, npos;
    mwSize columns, rows;
    int sparseflag; /* sparse */
    
    /* Check arguments */
    if (nrhs<2)
        mexErrMsgTxt("MINKMEX: Two input arguments required.");
    
    if (!mxIsNumeric(prhs[0]))
        mexErrMsgTxt("MINKMEX: First input LIST argument must be numeric.");
    
    sparseflag = mxIsSparse(prhs[0]); /* sparse */
    
    if (!mxIsNumeric(prhs[1]))                              
        mexErrMsgTxt("MINKMEX: Second input K must be numeric.");

    rows =  mxGetM(prhs[0]);
    columns = mxGetN(prhs[0]);
    nelem = rows*columns;
    /* Get the number of elements of the list of subindexes */
    if (sparseflag)
    {
        /* Search for the number of non-zero in sparse */      
        l = *(mxGetJc(prhs[0]) + columns);
    }
    else if (mxGetNumberOfDimensions(prhs[0])==2) /* Check for array */
        l = nelem;
    else
        mexErrMsgTxt("MINKMEX: First input LIST must be a 2D array.");
    
    if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS)
        mexErrMsgTxt("MINKMEX: First input LIST must be a double.");
    
    /* Get the number of elements of the list of subindexes */
    if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1)
        mexErrMsgTxt("MINKMEX: Second input K must be a scalar.");
    
    if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
        mexErrMsgTxt("MINKMEX: Second input K must be a double.");
    
    kout = k = (int)(*mxGetPr(prhs[1]));
    if (k<0)
        mexErrMsgTxt("MINKMEX: K must be non-negative integer.");
    
    /* Get a data pointer */
    list = mxGetPr(prhs[0]);

    /* Clip k */
    if (k>l) k=l;
    
    /* Clip kout */
    if (kout>nelem) kout=nelem;

    /* Clean programming */
    pos=NULL;
    
    /* Work for non-empty array */
    if (l>0) {
         /* Vector of index */
        pos = mxMalloc(sizeof(mwSize)*l);
        if (pos==NULL)
            mexErrMsgTxt("Out of memory.");
        /* Initialize the array of position (zero-based index) */
        for (i=0; i<l; i++) pos[i]=i;
        
        /* Call the recursive engine */
        k--; /* because we work on zero-based */
        nansplit = partNaN(0, l-1); /* Push NaN at the end */
        if (k<nansplit && nansplit>=0)
            findFirstK(0, nansplit);
        
        /* Look for the split of positive/negative numbers */
        if (sparseflag) {
            p0 = part0(0, k); /* number of strict negative elements */
            if (p0 < k) /* There are at least two positive elements */
            {
                /* Number of implicite zeros */
                nz = nelem-l;
                if (nz) /* in case the positive set is unordered */
                {
                    k -= nz;
                    findFirstK(p0, nansplit);
                    k += nz;
                }
            }
            /* ++ to restore one-based Matlab index */
            k++; 
        }
        else
            /* ++ to Restore one-based Matlab index */
            p0 = ++k;
    } /* if (l>0) */
    else p0 = 0;
    
    /* Number of implicite zero in (sparse) */
    nz = nelem-l;
    /* Create the Matrix result (first output) */
    plhs[0] = MinMaxResult(k, p0, nz, kout);
    
     /* Create the Matrix position (second output) */
    if (nlhs>=2)
    {
        if (sparseflag)
            SpLocResult(k, p0, nz, kout, prhs[0], &(plhs[1]), &(plhs[2]));
        else
            plhs[1] = LocResult(k, p0, nz, kout);
    }
    
    /* Free the array of position */
    if (pos) mxFree(pos);
    pos = NULL; /* clean programming */
    
    return;

} /* Gateway of minkmex.c */


