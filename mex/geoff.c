/*=================================================================
 *    test - cosine
 *
 * The calling syntax is:
 *
 *        [cx] = mycos(x)
 *
 * 
 *
 * This is a MEX-file for MATLAB.
 * Copyleft!!!!!!!!
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"


#if !defined(MAX)
#define    MAX(A, B)    ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define    MIN(A, B)    ((A) < (B) ? (A) : (B))
#endif

/* Input Arguments */

#define    X_IN    prhs[0]

#if !defined(MAX)
#define    MAX(A, B)    ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define    MIN(A, B)    ((A) < (B) ? (A) : (B))
#endif
/* Output Arguments */

#define    Y_OUT    plhs[0]



//CHANGE HERE!!!!
//this function calculates the derivatives
static void my_cosine(double X[], double Y[], unsigned int size)
{
    for(unsigned int i=0; i < size; i++){
	Y[i] = cos(X[i]);
    }    
    return;
}

//you have to check i/o sizes below, see err msgs
//main entry point
void mexFunction(int nlhs, mxArray *plhs[],
int nrhs, const mxArray*prhs[] )

{
    double *dy;
    double *t,*y;
    mwSize m,n;
   
    /* Check for proper number of arguments */
   
    if (nrhs != 1) {
        mexErrMsgTxt("I receive one array with the angles");
    } else if (nlhs > 1) {
        mexErrMsgTxt("too many output arrays");
    }
   
    m = mxGetM(X_IN);
    n = mxGetN(X_IN);
   
    /* Create a matrix for the return argument */
    Y_OUT = mxCreateDoubleMatrix(1, MAX(m,n), mxREAL);
   
    /* Assign pointers to the various parameters */
    y = mxGetPr(Y_OUT);
    x = mxGetPr(X_IN);
   
    /* Do the actual computations in a subroutine */
    my_cosine(x,y,MAX(m,n));

    return;
   
}


