// Calculate Sum1(X,Y,sLen) implmented in C++
// Output 1st is sum(sum(X./(Y+sLen))
// Output 2nd is sum(sum(X./(Y+sLen).^2)
// Copyright @ Guan Naiyang
#include "math.h"
#include "mex.h"   //--This one is required

void mexFunction(int nlhs, const mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //---Inside mexFunction---

	//Declarations
	const mxArray *xData;
    const mxArray *yData;
    const mxArray *zData;
	double *xValues, *yValues, *zValues, *outArray1, *outArray2;
	int index;
	int rowLen, colLen, len;
    double sum1st, sum2nd, temp;

	//Copy input pointer x
	xData = prhs[0];
    yData = prhs[1];
    zData = prhs[2];

	//Get matrix x
	xValues = mxGetPr(xData);
	rowLen = mxGetN(xData);
	colLen = mxGetM(xData);
    len = rowLen*colLen;
    yValues = mxGetPr(yData);
    zValues = mxGetPr(zData);

	//Allocate memory and assign output pointer
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //mxReal is our data-type
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL); //mxReal is our data-type

	//Get a pointer to the data space in our newly allocated memory
	outArray1 = mxGetPr(plhs[0]);
    outArray2 = mxGetPr(plhs[1]);

	//Calculate one by one
    sum1st = 0;
    sum2nd = 0;
	for(index=0; index<len; index++)
	{
        temp = xValues[index]/(yValues[index]+zValues[0]);
        sum1st = sum1st + temp;
        sum2nd = sum2nd + temp/(yValues[index]+zValues[0]);
	}
    
    outArray1[0] = sum1st;
    outArray2[0] = sum2nd;
    return;
}
