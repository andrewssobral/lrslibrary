/**************************************************************************
 *
 * File name: mexutils.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 15.8.2009
 *
 *************************************************************************/

#include "mexutils.h"
#include <math.h>



/* verify that the mxArray contains a double matrix */

void checkmatrix(const mxArray *param, char *fname, char *pname)
{
  char errmsg[100];
  sprintf(errmsg, "%.15s requires that %.25s be a double matrix.", fname, pname);
  if (!mxIsDouble(param) || mxIsComplex(param) || mxGetNumberOfDimensions(param)>2) {
    mexErrMsgTxt(errmsg);
  }
}


/* verify that the mxArray contains a 1-D double vector */

void checkvector(const mxArray *param, char *fname, char *pname)
{
  char errmsg[100];
  sprintf(errmsg, "%.15s requires that %.25s be a double vector.", fname, pname);
  if (!mxIsDouble(param) || mxIsComplex(param) || mxGetNumberOfDimensions(param)>2 || (mxGetM(param)!=1 && mxGetN(param)!=1)) {
    mexErrMsgTxt(errmsg);
  }
}


/* verify that the mxArray contains a double scalar */

void checkscalar(const mxArray *param, char *fname, char *pname)
{
  char errmsg[100];
  sprintf(errmsg, "%.15s requires that %.25s be a double scalar.", fname, pname);
  if (!mxIsDouble(param) || mxIsComplex(param) || mxGetNumberOfDimensions(param)>2 || 
      mxGetM(param)!=1 || mxGetN(param)!=1) 
  {
    mexErrMsgTxt(errmsg);
  }
}


/* verify that the mxArray contains a sparse matrix */

void checksparse(const mxArray *param, char *fname, char *pname)
{
  char errmsg[100];
  sprintf(errmsg, "%.15s requires that %.25s be sparse.", fname, pname);
  if (!mxIsSparse(param)) {
    mexErrMsgTxt(errmsg);
  }
}


/* verify that the mxArray contains a 1-dimensional cell array */

void checkcell_1d(const mxArray *param, char *fname, char *pname)
{
  char errmsg[100];
  sprintf(errmsg, "%.15s requires that %.25s be a 1-D cell array.", fname, pname);
  if (!mxIsCell(param) || mxGetNumberOfDimensions(param)>2 || (mxGetM(param)!=1 && mxGetN(param)!=1)) {
    mexErrMsgTxt(errmsg);
  }
}

