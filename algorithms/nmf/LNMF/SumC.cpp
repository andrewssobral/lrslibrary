// Calculate Sum(X,Y,sLen) implmented in C++
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
    int colIndex = 0;
    int rowIndex = 0;
    
    //Copy input pointer x
    xData = prhs[0];
    yData = prhs[1];
    zData = prhs[2];
    
    //Get matrix x
    xValues = mxGetPr(xData);
    colLen = mxGetN(xData);
    rowLen = mxGetM(xData);
    yValues = mxGetPr(yData);
    zValues = mxGetPr(zData);
    
    //Allocate memory and assign output pointer
    plhs[0] = mxCreateDoubleMatrix(colLen, 1, mxREAL); //mxReal is our data-type
    plhs[1] = mxCreateDoubleMatrix(colLen, 1, mxREAL); //mxReal is our data-type
    
    //Get a pointer to the data space in our newly allocated memory
    outArray1 = mxGetPr(plhs[0]);
    outArray2 = mxGetPr(plhs[1]);
    
    //Calculate one by one
    for(colIndex=0; colIndex<colLen; colIndex++)
    {
        sum1st = 0;
        sum2nd = 0;
        for(rowIndex=0; rowIndex<rowLen; rowIndex++)
        {
            index = (colIndex*rowLen)+rowIndex;
            temp = yValues[index]+zValues[colIndex];
            sum1st = sum1st + xValues[index]/temp;
            sum2nd = sum2nd + xValues[index]/temp/temp;
        }
        
        //Output the values
        outArray1[colIndex] = sum1st;
        outArray2[colIndex] = sum2nd;
    }

    return;
}
