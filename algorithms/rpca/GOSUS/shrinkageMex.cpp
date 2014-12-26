#include <mex.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
//#include <sys/time.h>
#include <matrix.h>
#include <string.h>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{ 
       
    double  *param, *prR, *prZ, *prZOutput;
    mwSize m, n, nzR;
    mwIndex *irR, *jcR, *irZ, *jcZ, *irZOutput, *jcZOutput;
    mxArray *Z;
    
    m=mxGetM( prhs[0] );
    n=mxGetN( prhs[0] );
    irR = mxGetIr( prhs[0] );
    jcR = mxGetJc( prhs[0] );
    prR = mxGetPr( prhs[0] );
    nzR = mxGetNzmax( prhs[0] );
     
    param = mxGetPr(prhs[1]);
    double eta = param[0];
    double rho = param[1];
    double threshold = 1 ;
    
    if ( rho>0 ) 
        threshold = eta/rho;
    
    Z = mxCreateSparse(m, n, nzR, mxREAL);
    irZ = mxGetIr( Z );
    jcZ = mxGetJc( Z );
    prZ = mxGetPr( Z );
    
    
    int indexR = 0;
    int indexZ = 0;
    jcZ[0] = 0;
    for( int i=0;i<n; i++){
         
        // Compute norm
        
        double rNorm = 0;
        
        int indexStart = indexR; // starting index for thresholding
        
        for(int j=0; j<jcR[i+1]-jcR[i]; j++){
            
            rNorm += prR[indexR]*prR[indexR];
            
            // cout << "( " << irR[indexR]+1  <<", "<< i+1 <<") : "<< prR[indexR]<<endl;
            indexR++;
        }
        
        rNorm = sqrt(rNorm);
        
        if (rNorm > threshold){
            
            // update jcZ: copy column number
            jcZ[i+1] = jcZ[i] + jcR[i+1]-jcR[i];
            
           double scale = (rNorm - threshold)/rNorm;
                
            for(int j=0; j<jcR[i+1]-jcR[i]; j++){
                              
                // update irZ: copy row number
                irZ[indexZ] = irR[indexStart]; 
                
                // update  prZ
                 prZ[indexZ] = prR[indexStart]*scale;
                 
               // cout << "( " << ir[index]+1  <<", "<< i+1 <<") : "<< pr[index]<<endl;
                indexZ++;
                indexStart++; //?????
            }
        
        }else{
          
            // update jcZ when z_i is set to zero after thresholding
            
            jcZ[i+1] = jcZ[i] ;
        }
            
        
         
    }
 
    // Prepare output
    
    plhs[0] =mxCreateSparse(m, n, indexZ, mxREAL);
    irZOutput = mxGetIr(  plhs[0] );
    jcZOutput = mxGetJc(  plhs[0] );
    prZOutput = mxGetPr(  plhs[0] );
    
    memcpy(irZOutput, irZ, indexZ * sizeof(mwIndex) );
    memcpy(jcZOutput, jcZ, (n+1)* sizeof(mwIndex) );
    memcpy(prZOutput, prZ, indexZ * sizeof(double) );
    
    //clean up;
    mxDestroyArray(Z);
    
     
}
