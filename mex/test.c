/*=================================================================
 *      Function that does interesting stuff.
 *
 *      The calling syntax is:
 *
 *        [matrices] = test(thetaArray)
 *
 *      This is a MEX-file for MATLAB.
 *=================================================================*/
#include <math.h>
#include "matrix.h"
#include "mex.h"

#if !defined(MAX)
#define    MAX(A, B)    ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define    MIN(A, B)    ((A) < (B) ? (A) : (B))
#endif

#define    ANGLES_IN     prhs[0]
#define    MATRICES_OUT  plhs[0]

static void doProduct(double anglesIn[], double matricesOut[], unsigned int noOfAnglesIn)
{
    unsigned int j;
    double M[4],N[4];
    
    for(j = 0; j < noOfAnglesIn; j++){
        double theta = anglesIn[j];
        // Consecutive matrices are spaced out by 4   
        M[1] = cos(theta);
        M[2] = -sin(theta);
        M[3] = sin(theta);
        M[4] = cos(theta);
        N[1] = theta * theta;
        N[2] = 5 * theta;
        N[3] = -1;
        N[4] = 0;
                
        matricesOut[j*4 + 0] = M[1] * N[1] + M[3] * N[2];
        matricesOut[j*4 + 1] = M[2] * N[1] + M[4] * N[2];
        matricesOut[j*4 + 2] = M[1] * N[3] + M[3] * N[4];
        matricesOut[j*4 + 3] = M[2] * N[3] + M[4] * N[4];
    }    
    return;
}

unsigned int errorCheckReturnLength(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    mwSize noOfColumns, noOfRows;
    
    if (nrhs != 1) {
        mexErrMsgTxt("I receive one array with the angles.");
    } else if (nlhs != 1) {
        mexErrMsgTxt("I return one variable.");
    }
    noOfColumns = mxGetM(ANGLES_IN);
    noOfRows = mxGetN(ANGLES_IN);
   
    if (MIN(noOfColumns, noOfRows) != 1) {
        mexErrMsgTxt("I receive one ARRAY with the angles.");
    }
    return MAX(noOfColumns, noOfRows);
}

//main entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *angles, *outputMatrices;
    unsigned int lengthOfInput;
    mwSignedIndex dims[3] = {2,2,0};
    
	lengthOfInput = errorCheckReturnLength(nlhs, plhs, nrhs, prhs);
    
    /* Create a 3D matrix for the return argument */
    dims[2] = lengthOfInput;
    MATRICES_OUT = mxCreateNumericArray(3, dims,  mxDOUBLE_CLASS, mxREAL);
   
    /* Assign pointers to the various parameters */
    outputMatrices = mxGetPr(MATRICES_OUT);
    angles = mxGetPr(ANGLES_IN);
   
    /* Do the actual computations in a subroutine */
    doProduct(angles, outputMatrices, lengthOfInput);

    return;
}


