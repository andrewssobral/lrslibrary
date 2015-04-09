#include "mex.h"

#include <cmath>
#include <cstring>
#define IN_OLDTAU prhs[0]
#define IN_NEWTAU prhs[1]
#define IN_ALPHA prhs[2]
#define IN_RANDOMNUM prhs[3]
#define IN_ROW prhs[4]
#define IN_COL prhs[5]

#define OUT_TAU plhs[0]

double dummy = -1;
inline double differ(double x, double y1, double y2){
    return  abs(log(x) - log(y1)) + abs(log(x) - log(y2));
}

inline double* getPos(double* data, int x, int y, int row, int col){
    if (x < 0 || y < 0 ) return &dummy;
    return data + x + y * row;
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]){
    double *oldTau, *newTau, *randomNum;
    double alpha = mxGetScalar(IN_ALPHA);
    double row = mxGetScalar(IN_ROW), col = mxGetScalar(IN_COL);
    
    oldTau = mxGetPr(IN_OLDTAU);
    newTau = mxGetPr(IN_NEWTAU);
    randomNum = mxGetPr(IN_RANDOMNUM);
    
    OUT_TAU = mxCreateDoubleMatrix(row, col, mxREAL);
    double *tau = mxGetPr(OUT_TAU);
    double ac = 0;
    
    memcpy(tau, oldTau, sizeof(double) * row * col);
    for (int j = 0; j < row; ++j)
        for (int k = 0; k < col; ++k)
        {
            if (k == 0 || j == 0){
                *getPos(tau, j, k, row, col) = *getPos(newTau, j, k, row, col);
                ++ac;
                continue;
            }
            if ( alpha * (differ(*getPos(tau, j, k, row, col), *getPos(tau, j-1, k, row, col), *getPos(tau, j, k-1, row, col)) 
            - differ(*getPos(newTau, j, k, row, col), *getPos(tau, j-1, k, row, col), *getPos(tau, j, k-1, row, col)) )
            > *getPos(randomNum, j, k, row, col)  ){
                *getPos(tau, j, k, row, col) = *getPos(newTau, j, k, row, col);
                ++ac;
            }
        }

        //mexPrintf("AC Ratio = %.2lf\n", ac / row / col * 100);
}