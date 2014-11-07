extern "C" {
#include "mexutil.h"
}
#include "string.h"

void DoWork(double *outp, double *inp, int sz, int ntrials);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	mxArray *result, *addrs;
	double *paddr;
	int i,niter;
	mwSize m, n;
	double *outp, *inp;
	mwSize sz;
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	niter = (int)*mxGetPr(prhs[1]);
	if(nlhs > 1) {
		addrs = mxCreateDoubleMatrix(niter,1,mxREAL);
		paddr = mxGetPr(addrs);
	}
	for(i=0;i<niter;i++) {
		if(i > 0) mxDestroyArray(result);
		result = mxCreateDoubleMatrix(m,n,mxREAL);
		outp = mxGetPr(result);
		//inp = mxGetPr(prhs[0]);
		//sz = mxGetNumberOfElements(prhs[0]);
		//memcpy(outp, inp, sz*sizeof(double));
		//DoWork(outp, inp, sz, 1);
		if(nlhs > 1) paddr[i] = (double)(long)outp;
	}
	if(nlhs > 0) plhs[0] = result;
	if(nlhs > 1) plhs[1] = addrs;
}
