#include "mex.h"
#include "pthread.h"
#include "stdio.h"
 
char *filename;
double *data;
size_t numElementsExpected;
size_t numElementsWritten;
 
/* thread compute function */ 
void *thread_run(void *p)
{
    /* Open the file for binary output */
    FILE *fp = fopen(filename, "wb");
    if (fp == NULL)
        mexErrMsgIdAndTxt("YMA:MexIO:errorOpeningFile", "Could not open file %s", filename);
 
    /* Write the data to file */
    numElementsWritten = (size_t) fwrite(data, sizeof(double), numElementsExpected, fp);
    fclose(fp);
 
    /* Ensure that the data was correctly written */
    if (numElementsWritten != numElementsExpected)
        mexErrMsgIdAndTxt("YMA:MexIO:errorWritingFile",
                "Error writing data to %s: wrote %d, expected %d\n", 
                filename, numElementsWritten, numElementsExpected);
 
    /* Cleanup */
    pthread_exit(NULL);
}
 
/* The MEX gateway function */
void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    pthread_t thread;
 
    /* Check for proper number of input and output arguments */
    if (nrhs != 2)
        mexErrMsgIdAndTxt("YMA:MexIO:invalidNumInputs", "2 input args required: filename, data");
    if (nlhs > 0)
        mexErrMsgIdAndTxt("YMA:MexIO:maxlhs", "Too many output arguments");
    if (!mxIsChar(prhs[0]))
        mexErrMsgIdAndTxt("YMA:MexIO:invalidInput", "Input filename must be of type string");
    if (!mxIsDouble(prhs[1]))
        mexErrMsgIdAndTxt("YMA:MexIO:invalidInput", "Input data must be of type double");
 
    /* Get the inputs: filename & data */
    filename = mxArrayToString(prhs[0]);
    data = mxGetPr(prhs[1]);  
    numElementsExpected = mxGetNumberOfElements(prhs[1]);
 
    /* Launch a new I/O thread using default attributes */
    if (pthread_create(&thread, NULL, thread_run, NULL))
        mexErrMsgIdAndTxt("YMA:MexIO:threadFailed", "Thread creation failed");
}