#include <math.h>
#include <matrix.h>
#include <mex.h>
#include "CMF_library.h"

#define yIn    prhs[0]
#define rIn    prhs[1]
#define tIn    prhs[2]
#define vIn    prhs[3]
#define uOut   plhs[0]

/*
 * Mex-Wrapper for the arc distance median filter for circle-valued images
 * with quantized data
 *
 * written by Martin Storath, 2016
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int M = mxGetM(yIn);
    int N = mxGetN(yIn);
    double* y = mxGetPr(yIn);
    double* v = mxGetPr(vIn);
    int S = mxGetM(vIn);
    int T = (int)mxGetScalar(tIn);
    int R = (int)mxGetScalar(rIn);
    uOut = mxCreateDoubleMatrix(M, N, mxREAL);
    double* u = mxGetPr(uOut);
    medfiltCirc2DQuant(u, y, M, N, R, T, v, S);
}