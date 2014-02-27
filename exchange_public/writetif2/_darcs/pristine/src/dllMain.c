#include <windows.h>
#include <math.h>
#include <Rts.h>
#include "mex.h"
#include "Wtifc_stub.h"
#define __ADDER_DLL_EXPORT
#define ADDER_API _declspec(dllexport)
#define MAX_FILENAME    300
extern void __stginit_Wtifc(void);

static char* args[] = { "ghcDll", NULL };
/* N.B. argv arrays must end with NULL */

BOOL STDCALL DllMain( HANDLE hModule, 
		      DWORD  ul_reason_for_call, 
		      LPVOID lpReserved
                      ){
  return TRUE;
}

ADDER_API BOOL adder_Begin(){
  startupHaskell(1, args, __stginit_Wtifc);
  return HS_BOOL_TRUE;
}

ADDER_API void adder_End(){
  shutdownHaskell();
}


/* Input Arguments */

#define	T_IN	prhs[1]
#define	XY_IN	prhs[2]
#define	Y_IN	prhs[0]


/* Output Arguments */

#define	YP_OUT	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

static	double	mu = 1/82.45;
static	double	mus = 1 - 1/82.45;

// See %MATLAB%\extern\mx\mxisfinit.c
static int dtoi32(double d)
{
/* Function that converts double to int32 */
    int i=0;
    
    if(mxIsFinite(d)) {
	if(d < (double)INT_MAX && d > (double)INT_MIN) {
	    i = (int) d;
	} else {
	    i =  ((d > 0) ? INT_MAX : INT_MIN);
	}	    
    } else if(mxIsInf(d)) {
	i = ( (d > 0) ? INT_MAX : INT_MIN);
	/* NOTE: Test for NaN is here for illustration only.  If a double is
	   not finite and is not infinity, then it is a NaN */
    } else if(mxIsNaN(d)) {
	mexWarnMsgTxt("dtoi32: NaN detected.  Translating to 0.\n");
	i = 0;
    }
    return i;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
/*     double *yp;  */
    char t[MAX_FILENAME];
/*     double *t,*y;  */
    mwSize m,n; 

    adder_Begin();
    
/*     /\* Check for proper number of arguments *\/ */
    
    if (nrhs != 3) {
	mexErrMsgTxt("Three input arguments required.");
    }

/*     /\* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. *\/  */
    
    m = mxGetM(Y_IN);
    n = mxGetN(Y_IN);
    if (!mxIsChar(T_IN))
	mexErrMsgTxt("Missing filename");
    if (!mxIsUint8(Y_IN))
	mexErrMsgTxt("Invalid input arguments: first argument should be a M-by-N uint8 (grayscale image) array");
    
    double* pr = mxGetPr(XY_IN); 

    /* float xyf[2]; */
    /* int i; */

    /* for(i=0; i < 2; i++) { */
    /* 	xyf[i] = pr[i]; // Ok */
    /* } */

/*     /\* Assign pointers to the various parameters *\/  */
/*     yp = mxGetPr(YP_OUT); */
    
    mxGetString(T_IN, t, MAX_FILENAME - 1);
    unsigned char* y = (unsigned char*) mxGetData(Y_IN);
        
    /* Do the actual computations in a subroutine */
/*     yprime(yp,t,y);  */

    /* wtifc(t, y, xyf, n, m); */
    wtifc(t, y, pr, n, m);

// if this line presents MATLAB does not crash on ``clear functions'' but allows to 
// evaluate MEX-file only once
// Second eval causes crash..
// if there's no adder_End it's possibe to call MEX several times though 
// MATLAB crashes on ``clear functions''
// PS ``clear functions'' effectively unloads every MEX-file ..

/*     adder_End(); */


    return;
    
}


/* Local Variables:  */
/* compile-command:"\"c:/Program Files/Microsoft Visual Studio 9.0/VC/BIN/nmake\"" */
/* End: */
