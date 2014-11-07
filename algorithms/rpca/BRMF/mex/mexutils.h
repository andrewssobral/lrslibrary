/**************************************************************************
 *
 * File name: mexutils.h
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 18.8.2009
 *
 * Utility functions for MEX files.
 *
 *************************************************************************/


#ifndef __MEX_UTILS_H__
#define __MEX_UTILS_H__

#include "mex.h"



/**************************************************************************
 * Function checkmatrix:
 *
 * Verify that the specified mxArray is real, of type double, and has 
 * no more than two dimensions. If not, an error message is printed
 * and the mex file terminates.
 * 
 * Parameters:
 *   param - the mxArray to be checked
 *   fname - the name of the function where the error occured (15 characters or less)
 *   pname - the name of the parameter (25 characters or less)
 *
 **************************************************************************/
void checkmatrix(const mxArray *param, char *fname, char *pname);


/**************************************************************************
 * Function checkvector:
 *
 * Verify that the specified mxArray is 1-D, real, and of type double. The
 * vector may be a column or row vector. Otherwise, an error message is
 * printed and the mex file terminates.
 * 
 * Parameters:
 *   param - the mxArray to be checked
 *   fname - the name of the function where the error occured (15 characters or less)
 *   pname - the name of the parameter (25 characters or less)
 *
 **************************************************************************/
void checkvector(const mxArray *param, char *fname, char *pname);


/**************************************************************************
 * Function checkscalar:
 *
 * Verify that the specified mxArray represents a real double scalar value. 
 * If not, an error message is printed and the mex file terminates.
 * 
 * Parameters:
 *   param - the mxArray to be checked
 *   fname - the name of the function where the error occured (15 characters or less)
 *   pname - the name of the parameter (25 characters or less)
 *
 **************************************************************************/
void checkscalar(const mxArray *param, char *fname, char *pname);


/**************************************************************************
 * Function checksparse:
 *
 * Verify that the specified mxArray contains a sparse matrix. If not,
 * an error message is printed and the mex file terminates.
 * 
 * Parameters:
 *   param - the mxArray to be checked
 *   fname - the name of the function where the error occured (15 characters or less)
 *   pname - the name of the parameter (25 characters or less)
 *
 **************************************************************************/
void checksparse(const mxArray *param, char *fname, char *pname);


/**************************************************************************
 * Function checkcell_1d:
 *
 * Verify that the specified mxArray is a 1-D cell array. The cell array 
 * may be arranged as either a column or a row. If not, an error message 
 * is printed and the mex file terminates.
 * 
 * Parameters:
 *   param - the mxArray to be checked
 *   fname - the name of the function where the error occured (15 characters or less)
 *   pname - the name of the parameter (25 characters or less)
 *
 **************************************************************************/
void checkcell_1d(const mxArray *param, char *fname, char *pname);


#endif

