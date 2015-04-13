/*  Grassmann Averages
    Copyright (C) 2014  Søren Hauberg

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "trimmed_mean.h"
#include "mex.h"
#include <omp.h>
//#include <cstdio>

template<class FLOAT>
void estimate(const FLOAT *data, const mwIndex D, const mwIndex N, const int percent, FLOAT *output)
{
  // Perform the estimation
  #pragma omp parallel for num_threads(omp_get_num_procs())
  for (mwIndex d = 0; d < D; d++)
    {
      output[d] = trimmed_mean<FLOAT, mwIndex> (data, d*N, (d+1)*N, percent);
    }
}

extern "C"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check input
  if (nrhs < 2)
    {
      // XXX: What's the best way of raising errors?
      mexErrMsgTxt("trimmed_mean: not enough input arguments\n");
      return; // XXX: Do we need this return?
    }
  
  // Get data
  const int N = mxGetM(prhs[0]);
  const int D = mxGetN(prhs[0]);  
  const double dpercent = mxGetScalar(prhs[1]);
  const int percent = (int)((double)N*dpercent/100.0);

  const mxClassID matrix_class = mxGetClassID(prhs[0]);
  if (matrix_class == mxDOUBLE_CLASS)
    {
      const double *data = mxGetPr(prhs[0]);
  
      // Allocate output
      mxArray *retval = mxCreateDoubleMatrix(1, D, mxREAL);
      double *output = mxGetPr(retval);
 
      // Do the magic and return
      estimate(data, D, N, percent, output);
      plhs[0] = retval;
    }
  else if (matrix_class == mxSINGLE_CLASS)
    {
      const float *data = (float*)mxGetData(prhs[0]);
  
      // Allocate output
      mxArray *retval = mxCreateNumericMatrix(1, D, matrix_class, mxREAL);
      float *output = (float*)mxGetData(retval);
 
      // Do the magic and return
      estimate(data, D, N, percent, output);
      plhs[0] = retval;
    }
  else
    {
      mexErrMsgTxt("trimmed_mean: first input arugment must be either a double or a single precision matrx\n");
    }
}

