// Kullback-Leibler divergence implmented in C++
// Copyright @ Guan Naiyang
#include "math.h"
#include "mex.h"   //--This one is required

void mexFunction(int nlhs, const mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //---Inside mexFunction---

	//Declarations
	const mxArray *xData;
    const mxArray *yData;
	double *xValues, *yValues, *outArray;
	int index;
	int rowLen, colLen, len;
    double sumValues;

	//Copy input pointer x
	xData = prhs[0];
    yData = prhs[1];

	//Get matrix x
	xValues = mxGetPr(xData);
    yValues = mxGetPr(yData);
	rowLen = mxGetN(xData);
	colLen = mxGetM(xData);
    len = rowLen*colLen;

	//Allocate memory and assign output pointer
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //mxReal is our data-type

	//Get a pointer to the data space in our newly allocated memory
	outArray = mxGetPr(plhs[0]);

	//Calculate one by one
    sumValues = 0;
	for(index=0; index<len; index++)
	{
        sumValues = sumValues + xValues[index]*log(xValues[index]/yValues[index]) - xValues[index] + yValues[index];
	}
    
    outArray[0] = sumValues;
    return;
}
