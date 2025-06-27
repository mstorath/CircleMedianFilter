#include <math.h>
#include "matrix.h"
#include <mex.h>
#include "CMF_library.h"

#define yIn         prhs[0]
#define rIn         prhs[1]
#define tIn         prhs[2]
#define nIterIn     prhs[3]
#define stopTolIn   prhs[4]
#define uOut        plhs[0]

/*
 * Mex-Wrapper for the geometric median filter in R^N
 *
 * written by Martin Storath, 2016
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int RD = (int)mxGetScalar(rIn);
    int TD = (int)mxGetScalar(tIn);
    int nIter = (int)mxGetScalar(nIterIn);
    double stopTol = mxGetScalar(stopTolIn);
    int nDimNum = mxGetNumberOfDimensions(yIn);
    const int* pDims = mxGetDimensions(yIn);
    double* y = mxGetPr(yIn);
    uOut = mxCreateNumericArray(nDimNum, pDims, mxDOUBLE_CLASS, mxREAL);
    double* u = mxGetPr(uOut);
    int M = pDims[0];
    int N = pDims[1];
    int K = pDims[2];
    medfiltGeoRN2D(u, y, M, N, K, RD, TD, nIter, stopTol);
}
