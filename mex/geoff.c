/*=================================================================
 *    Function that does interesting stuff.
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
#include <pcap.h>


#if !defined(MAX)
#define    MAX(A, B)    ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define    MIN(A, B)    ((A) < (B) ? (A) : (B))
#endif

/* Input Arguments */
#define    X_IN    prhs[0]


//CHANGE HERE!!!!
//this function calculates the derivatives
static void my_cosine(double X[], double Y[], unsigned int m, unsigned int n)
{
    unsigned int i, j;
    
    for(j = 0; j < n; j++){
        for(i = 0; i < m; i++){
            Y[m*i + j] = cos(X[m*i + j]);
        }
    }    
    return;
}

//you have to check i/o sizes below, see err msgs
//main entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{
    double *x,*y;
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
    plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
   
    /* Assign pointers to the various parameters */
    y = mxGetPr(plhs[0]);
    x = mxGetPr(X_IN);
   
    /* Do the actual computations in a subroutine */
    my_cosine(x, y, m, n);

    pcap_if_t *alldevs;
    pcap_if_t *d;
    int i=0;
    char errbuf[PCAP_ERRBUF_SIZE];
    
    /* Retrieve the device list from the local machine */
    if (pcap_findalldevs_ex(PCAP_SRC_IF_STRING, NULL /* auth is not needed */, &alldevs, errbuf) == -1)
    {
        fprintf(stderr,"Error in pcap_findalldevs_ex: %s\n", errbuf);
        exit(1);
    }
    
    /* Print the list */
    for(d= alldevs; d != NULL; d= d->next)
    {
        printf("%d. %s", ++i, d->name);
        if (d->description)
            printf(" (%s)\n", d->description);
        else
            printf(" (No description available)\n");
    }
    
    if (i == 0)
    {
        printf("\nNo interfaces found! Make sure WinPcap is installed.\n");
        return;
    }

    /* We don't need any more the device list. Free it */
    pcap_freealldevs(alldevs);
    
    return;
   
}


