#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <random>
#include <map>
#include <omp.h>
#include <algorithm>
#include <random>
#include <queue>
#include "mex.h"

using namespace std;

// Sparse Matrix Multiplication Code
// Arguments:
// 1) U: nRows x r RHS of Matrix Multiplication
// 2) V: r x nCols LHS of Matrix Multiplication
// 3) Inds: n x 2 matrix of indices of entries to be computed
// 4) nThreads: Number of threads
// Return Value:
// 1) Values: n x 1 matrix of values

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int r = mxGetN(prhs[0]), m = mxGetM(prhs[0]), 
		n = mxGetN(prhs[1]), nTh = mxGetScalar(prhs[3]), 
		nObs = mxGetM(prhs[2]);

	double *uPr = mxGetPr(prhs[0]), *vPr = mxGetPr(prhs[1]), *indPr = mxGetPr(prhs[2]);

	plhs[0] = mxCreateDoubleMatrix(nObs, 1, mxREAL);

	double *retPr = mxGetPr(plhs[0]);

	omp_set_num_threads(nTh);

#pragma omp parallel for
	for (long long ind = 0; ind < nObs; ind++) {
		int rInd = indPr[ind] - 1, cInd = indPr[nObs + ind] - 1;
		for (int feat = 0; feat < r; feat++) {
			retPr[ind] += uPr[feat*m + rInd] * vPr[cInd*r + feat];
		}
	}
}