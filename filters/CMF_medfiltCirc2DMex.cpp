#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdexcept>
#include "CMF_library.h"

#define yIn    prhs[0]
#define rIn    prhs[1]
#define tIn    prhs[2]
#define uOut   plhs[0]

/*
 * Mex-Wrapper for the arc distance median filter for circle-valued images
 * with non-quantized data
 *
 * written by Martin Storath, 2016
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int M = mxGetM(yIn);
    int N = mxGetN(yIn);
    double* y = mxGetPr(yIn);
    int T = (int)mxGetScalar(tIn);
    int R = (int)mxGetScalar(rIn);
    uOut = mxCreateDoubleMatrix(M, N, mxREAL); // create output       
    double* u = mxGetPr(uOut); // get output array
    medfiltCirc2D(u, y, M, N, R, T);
}
