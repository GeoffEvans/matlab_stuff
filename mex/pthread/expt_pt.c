#include <math.h>
#include "mex.h"
#include "pthread.h"
#include "stdio.h"
 
#define XR(i,j) xr[i+4*j]

char *filename;
double *data;
size_t numElementsExpected;
size_t numElementsWritten;
mxArray *x, *pa, *xx;
double handle;

static void fill_array( double	*xr, int n )
{
    int i,j;
    /* Fill real and imaginary parts of array. */
    for (j = 0; j < 4; j++) {
        for (i = 0; i <= j; i++) {
            XR(i,j) = 4 + i - j;
            XR(j,i) = XR(i,j) + n;
        }
    }
}

void *thread_run(void *p)
{  
    int k;
    double *v;
    v = mxMalloc(1* sizeof(double));
    v[0] = 0.5;
    xx = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetPr(xx, v);
    
    for (k = 0; k < 10; k++)
    {        
        fill_array(mxGetPr(x), k);
        mexCallMATLAB(0, NULL, 1, &xx, "pause");
		mexSet(handle, "CDATA", x);
    }
    mxDestroyArray(x);
    mxDestroyArray(xx);
    pthread_exit(NULL);
}

void mexFunction(int nlhs,       mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    pthread_t thread;
    mwSize m;
    pa = prhs[0];
    handle = mxGetScalar(pa);

    m = 4;    
    x =  mxCreateDoubleMatrix(m, m, mxREAL);
    
    if (pthread_create(&thread, NULL, thread_run, NULL))
        mexErrMsgIdAndTxt("YMA:MexIO:threadFailed", "Thread creation failed");
}
