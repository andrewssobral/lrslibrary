/* C implementation of repmat.
 * by Tom Minka
 */
/*
mex -c mexutil.c
mex repmat.c mexutil.obj
to check for warnings:
gcc -Wall -I/cygdrive/c/Program\ Files/MATLAB/R2008a/extern/include -c repmat.c
*/
#include "mexutil.h"
#include <string.h>

/* ALWAYS_2D == 1 gives the usual matlab behavior where repmat(0,3) is
 * a 3x3 matrix.
 * ALWAYS_2D == 0 gives a 3x1 matrix instead.
 * repmat(x,sizes,1) is another way to get this behavior.
 */
#define ALWAYS_2D 1

/* Repeat a block of memory rep times.
 * The characters from dest[0] .. dest[chunk-1] are copied to
 * dest[chunk] .. dest[2*chunk-1] as well as
 * dest[2*chunk] .. dest[3*chunk-1] and so on up to
 * dest[(rep-1)*chunk] .. dest[(rep-1)*chunk + chunk-1]
 */
void memrep(char *dest, mwSize chunk, mwSize rep)
{
  /* printf("chunk = %d, rep = %d\n", chunk, rep); */
  if(chunk >= 1024) {
    /* big chunks */
    mwSize i;
    char *p = dest;
    for(i=1;i<rep;i++) {
      p += chunk;
      memcpy(p, dest, chunk);
      /* could also use: memcpy(p, p-chunk, chunk); */
    }
  } else {
    /* small chunks */
    if(rep == 1) return;
    memcpy(dest + chunk, dest, chunk); 
    if(rep & 1) {  /* odd number of repetitions */
      dest += chunk;
      memcpy(dest + chunk, dest, chunk);
    }
    /* now repeat using a block twice as big */
    memrep(dest, chunk<<1, rep>>1);
  }
}

void repmat(char *dest, const char *src, int ndim, mwSize *destdimsize, 
	    mwSize *dimsize, const mwSize *dims, mwSize *rep) 
{
  int d = ndim-1;
  mwSize i;
  mwSize chunk;
  /* copy the first repetition into dest */
  if(d == 0) {
    chunk = dimsize[0];
    memcpy(dest,src,chunk);
  }
  else {
    /* recursively repeat each slice of src */
    for(i=0;i<dims[d];i++) {
      repmat(dest + i*destdimsize[d-1], src + i*dimsize[d-1], 
	     ndim-1, destdimsize, dimsize, dims, rep);
    }
    chunk = destdimsize[d-1]*dims[d];
  }
  /* copy the result rep-1 times */
  memrep(dest,chunk,rep[d]);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  const mxArray *srcmat;
  int ndim, eltsize;
  mwSize *dimsize;
  const mwSize *dims;
  int ndimdest;
  mwSize *destdims, *destdimsize;
  char *src, *dest;
  mwSize *rep;
  int i,nrep;
  int extra_rep = 1;
  int empty;
	double *outp, *inp;
	mwSize m,n,numel;

  if(nrhs < 2) mexErrMsgTxt("Usage: repmat(A, [M N ...])");
  srcmat = prhs[0];
	if(0) {
		/* testing code, please ignore */
		/*
		plhs[0] = mxCreateNumericArrayE(ndim, dims, mxGetClassID(srcmat), 
																		mxIsComplex(srcmat)?mxCOMPLEX:mxREAL);
																		*/
		m = mxGetM(srcmat);
		n = mxGetN(srcmat);		
		plhs[0] = mxCreateDoubleMatrixE(m, n, mxREAL);
		outp = mxGetPr(plhs[0]);
		inp = mxGetPr(srcmat);
		numel = mxGetNumberOfElements(srcmat);
		memcpy(outp, inp, numel*sizeof(double));
		/*
		plhs[0] = mxCreateNumericMatrixE(dims[0], dims[1], mxGetClassID(srcmat), 
																		mxIsComplex(srcmat)?mxCOMPLEX:mxREAL);
		plhs[0] = mxCreateNumericMatrix(0, 0, mxGetClassID(srcmat), 
																		mxIsComplex(srcmat)?mxCOMPLEX:mxREAL);
		*/
		return;
	}

  if(!mxIsNumeric(srcmat) || mxIsSparse(srcmat) || mxIsCell(srcmat) || mxIsStruct(srcmat)) {
    /* call Matlab's repmat */
    mexCallMATLAB(nlhs,plhs,nrhs,(mxArray**)prhs,"xrepmat");return;
    /* mexErrMsgTxt("Sorry, can't handle sparse matrices yet."); */
  }
  ndim = mxGetNumberOfDimensions(srcmat);
  dims = mxGetDimensions(srcmat);
  eltsize = mxGetElementSize(srcmat);
	
  /* compute dimension sizes */
  dimsize = (mwSize*)mxCalloc(ndim, sizeof(mwSize));
  dimsize[0] = eltsize*dims[0];
  for(i=1;i<ndim;i++) dimsize[i] = dimsize[i-1]*dims[i];

  /* determine repetition vector */
  ndimdest = ndim;
  if(nrhs == 2) {
    /* prhs[1] is a vector of reps */
    nrep = mxGetN(prhs[1]);
    if(nrep > ndimdest) ndimdest = nrep;
    rep = (mwSize*)mxCalloc(ndimdest, sizeof(mwSize));
    for(i=0;i<nrep;i++) {
      double repv = mxGetPr(prhs[1])[i];
      rep[i] = (mwSize)repv;
    }
#if ALWAYS_2D
    if(nrep == 1) {
      /* special behavior */
      nrep = 2;
      rep[1] = rep[0];
    }
#endif
  }
  else {
    /* concatenate all prhs's */
    int ri=0;
    nrep = 0;
    for(i=0;i<nrhs-1;i++) {
      nrep += mxGetNumberOfElements(prhs[i+1]);
    }
    if(nrep > ndimdest) ndimdest = nrep;
    rep = (mwSize*)mxCalloc(ndimdest, sizeof(mwSize));
    for(i=0;i<nrhs-1;i++) {
      double *p = mxGetPr(prhs[i+1]);
      int j, sz = mxGetNumberOfElements(prhs[i+1]);
      for(j=0;j<sz;j++) rep[ri++] = (mwSize)p[j];
    }
  }
  for(i=nrep;i<ndimdest;i++) rep[i] = 1;

  /* compute output size */
  destdims = (mwSize*)mxCalloc(ndimdest, sizeof(mwSize));
  for(i=0;i<ndim;i++) destdims[i] = dims[i]*rep[i];
  for(;i<ndimdest;i++) { 
    destdims[i] = rep[i];
    extra_rep *= rep[i];
  }
  destdimsize = (mwSize*)mxCalloc(ndim, sizeof(mwSize));
  destdimsize[0] = eltsize*destdims[0];
  for(i=1;i<ndim;i++) destdimsize[i] = destdimsize[i-1]*destdims[i];

  /* for speed, array should be uninitialized */
  plhs[0] = mxCreateNumericArrayE(ndimdest, destdims, mxGetClassID(srcmat), 
				  mxIsComplex(srcmat)?mxCOMPLEX:mxREAL);

  /* if any rep[i] == 0, output should be empty array.
     Added by KPM 11/13/02.
  */
  empty = 0;
  for (i=0; i < nrep; i++) {
    if (rep[i]==0) 
      empty = 1;
  }
  if (empty) 
    return;

  src = (char*)mxGetData(srcmat);
  dest = (char*)mxGetData(plhs[0]);
  repmat(dest,src,ndim,destdimsize,dimsize,dims,rep);
  if(ndimdest > ndim) memrep(dest,destdimsize[ndim-1],extra_rep);
  if(mxIsComplex(srcmat)) {
    src = (char*)mxGetPi(srcmat);
    dest = (char*)mxGetPi(plhs[0]);
    repmat(dest,src,ndim,destdimsize,dimsize,dims,rep);
    if(ndimdest > ndim) memrep(dest,destdimsize[ndim-1],extra_rep);
  }
}
