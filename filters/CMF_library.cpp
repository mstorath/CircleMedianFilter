#include <math.h>
#include <iostream>
#include "time.h"
#include <stdio.h>

/*
 * Reference implementation of the arc distance median filters
 * for quantized and non-quantized data
 * described in the paper
 * Martin Storath, Andreas Weinmann 
 * "Fast median filtering for phase or orientation data"
 * IEEE Transactions on Pattern Analysis and Machine Intelligence, 2017
 *
 * written by Martin Storath, 2016
 */

/*
 *precomputation of 2\pi
 */
const double M_PI2 = 2.0 * M_PI;

/*
 * compute antipodal point of point a
 * a must be in the interval (-pi, pi]
 */
inline double getAntipodalPoint(double a) {
    if (a <= 0) {
        return a + M_PI;
    } else {
        return a - M_PI;
    }
}

/*
 * distance of the circle between a and b
 * a and b must be in the interval (-pi, pi]
 */
inline double dist(double a, double b) {
    double aux = fabs(a - b);
    return (aux >  M_PI) ? M_PI2 - aux : aux;
}

/*
 * get mirror boundary conditions as value
 */
inline double getMB(double* y, int i, int j, int M, int N){
    if (i >= M) {
        i = 2 * M - i - 1;
    } else if (i < 0) {
        i =  -i-1;
    }
    if (j >= N) {
        j = 2 * N -j -1;
    } else if (j < 0) {
        j = -j-1;
    }
    return y[j * M + i];
}

/*
 * get mirror boundary conditions as index
 */
inline int getMBInd(int j, int N){
    if (j >= N) {
        j = 2 * N -j -1;
    } else if (j < 0) {
        j = -j-1;
    }
    return j;
}

/*
 * check for valid dimensions
 */
inline void checkDims(int M, int N, int R, int T) {
    bool filterOdd = ((R % 2) == 1) && ((T % 2) == 1);
    if (!filterOdd) {
        throw std::invalid_argument( "Invalid input: filter size must be odd.");
    } 
    bool filterValid = (M >= R) && (N >= T);
    if (!filterValid) {
        throw std::invalid_argument( "Invalid input: filter cannot be larger than the image.");
    }
}


/*
 * Circle median filter for (non-quantized) data y
 * 
 * y: data in M x N, all elements must be in the interval (-pi, pi]
 *      (the implementation does not check that!)
 * u: output array M x N,
 * M, N: size of data
 * RD, TD: size of filter mask (must be odd integers)
 *
 */
void medfiltCirc2D(double* u, double* y, int M, int N, int RD, int TD) {  
    // check dimensions
    checkDims(M, N, RD, TD); 
 
    // init
    int R = (RD-1)/2;
    int T = (TD-1)/2;
    
    // create auxiliary values and arrays
    double  dAux, yAux;
    double  yAux2, minAux, minY, fAux,
            yInAux, yOutAux,
            yPlus1, yPlus2, yMinus1, yMinus2;
    double* Z = new double[RD * TD];
    
    // index steps for G
    double* G = new double[N * RD * TD];
    int g1 = RD * TD;
    int g2 = TD;
    int kG, lG;
    
    // compute G^00
    for (int k = -R; k <= R; k++) {
        kG = k + R;
        for (int l = -T; l <= T; l++) {
            dAux = 0;
            yAux = getMB(y, k, l, M, N);
            for (int r = -R; r <= R; r++) {
                for (int t = -T; t <= T; t++) {
                    yAux2 = getMB(y, r, t, M, N);
                    dAux += dist(yAux, yAux2);
                }
            }
            lG = l + T;
            G[kG * g2 + lG] = dAux;
        }
    }
    
    // compute G^0n for n = 1,...,N
    for (int n = 1; n < N; n++) {
        for (int k = -R; k <= R; k++) {
            kG = k + R;
            for (int l = -T; l <= T-1; l++) {
                dAux = 0;
                yAux = getMB(y, k, n+l, M, N);
                for (int r = -R; r <= R; r++) {
                    yInAux = getMB(y, r, n+T, M, N);
                    yOutAux = getMB(y, r, n-T-1, M, N);
                    dAux += dist(yAux, yInAux) - dist(yAux, yOutAux);
                }
                lG = l + T;
                G[n * g1 + kG * g2 + lG] = G[ (n-1) * g1 + kG * g2 + (lG+1)] + dAux;
            }
        }
        
        // for k= -R,.., R and l=T
        lG = 2 * T;
        for (int k = -R; k <= R; k++) {
            dAux = 0;
            yAux = getMB(y, k, n+T, M, N);
            for (int r = -R; r <= R; r++) {
                for (int t = -T; t <= T; t++) {
                    yAux2 = getMB(y, r, n+t, M, N);
                    dAux += dist(yAux, yAux2);
                }
            }
            kG = k + R;
            G[n * g1 + kG * g2 + lG] = dAux;
        }
    }
    
    // compute medians for first row
    for (int n = 0; n < N; n++) {
        minAux = INFINITY;
        minY = 0;
        for (int k = -R; k <= R; k++) {
            kG = k + R;
            for (int l = -T; l <= T; l++) {
                // min of G
                lG = l + T;
                fAux = G[n * g1 + kG * g2 + lG];
                if (fAux < minAux) {
                    minY = getMB(y, k, n+l, M, N);
                    minAux = fAux;
                }
            }
        }
        u[n*M] = minY;
    }
    
    for (int m = 1; m < M; m++) {
        
        // first column
        int n = 0;
        
        // init Z
        for (int k = -R; k <= R; k++) {
            kG = k + R;
            for (int l = -T; l <= T; l++) {
                dAux = 0;
                yAux = getMB(y, m+k, n+l, M, N);
                for (int t = -T; t <= T; t++) {
                    yInAux = getMB(y, m+R, n+t, M, N);
                    yOutAux = getMB(y, m-R-1, n+t, M, N);
                    dAux += dist(yAux, yInAux) - dist(yAux, yOutAux);
                }
                lG = l + T;
                Z[kG * g2 + lG] = dAux;
            }
        }
        
        // init G for n = 0 (first column)
        for (int l = -T; l <= T; l++) {
            lG = l + T;
            for (int k = -R; k <= R-1; k++) {
                kG = k + R;
                G[n*g1 + kG*g2 + lG] = G[n*g1 + (kG+1)*g2 + lG] + Z[kG*g2 + lG];
            }
            dAux = 0;
            yAux = getMB(y, m+R, n+l, M, N);
            for (int r = -R; r <= R; r++) {
                for (int t = -T; t <= T; t++) {
                    dAux += dist( getMB(y, m+r, n+t, M, N), yAux);
                }
            }
            kG = 2*R;
            G[n*g1 + kG*g2 + lG] = dAux;
        }
        
        // process row m and columns n=1,..., N-1
        for (int n = 1; n < N; n++) {
            
            int nG = n*g1;
           
            // update Z 
            yPlus1 = getMB(y, m+R, n+T, M, N);
            yPlus2 = getMB(y, m-R-1, n-T-1, M, N);
            yMinus1 = getMB(y, m-R-1, n+T, M, N);
            yMinus2 = getMB(y, m+R, n-T-1, M, N);
            
            for (int k = -R; k <= R; k++) {
                kG = k + R;
                for (int l = -T; l <= T-1; l++) {
                    lG = l + T;
                    yAux = getMB(y, m+k, n+l, M, N);
                    Z[kG  * g2 + lG] = Z[kG  * g2 + lG+1]
                            + dist(yAux, yPlus1)  + dist(yAux, yPlus2)
                            - dist(yAux, yMinus1) - dist(yAux, yMinus2);
                }
                dAux = 0;
                yAux = getMB(y, m+k, n+T, M, N);
                for (int t = -T; t <= T; t++) {
                    dAux +=   dist(yAux, getMB(y, m+R, n+t, M, N)) 
                            - dist(yAux, getMB(y, m-R-1, n+t, M, N));
                }
                lG = 2 * T;
                Z[kG * g2 + lG] = dAux;
            }

            // compute G for k = 1,...,R-1, and l = 1,..., T
            for (int l = -T; l <= T; l++) {
                lG = l + T;
                for (int k = -R; k <= R-1; k++) {
                    kG = k + R;
                    G[nG + kG*g2 + lG] = G[n*g1 + (kG+1)*g2 + lG] + Z[kG*g2 + lG];
                }
            }
            
            // last line of G
            kG = 2*R;
            for (int l = -T; l <= T-1; l++) {
                yAux = getMB(y, m+R, n+l, M, N);
                dAux = 0;
                for (int r = -R; r <= R; r++) {
                    dAux +=   dist(yAux, getMB(y, m+r, n+T, M, N) )
                             -dist(yAux, getMB(y, m+r, n-T-1, M, N) );
                }
                lG = l + T;
                G[nG + kG*g2 + lG] = G[(n-1)*g1 + kG*g2 + lG+1] + dAux;
            }
            
            // last element of G by simple summation
            dAux = 0;
            yAux = getMB(y, m+R, n+T, M, N);
            for (int r = -R; r <= R; r++) {
                for (int t = -T; t <= T; t++) {
                    yAux2 = getMB(y, m+r, n+t, M, N);
                    dAux += dist(yAux, yAux2);
                }
            }
            G[nG + 2*R*g2 + 2*T] = dAux;
        }

        // compute medians for row m and columns n = 0,..., N-1
        for (int n = 0; n < N; n++) {
            minAux = INFINITY;
            minY = 0;
            for (int k = -R; k <= R; k++) {
                kG = k + R;
                for (int l = -T; l <= T; l++) {
                    // min of G
                    lG = l + T;
                    fAux = G[n * g1 + kG * g2 + lG];
                    if (fAux < minAux) {
                        minY = getMB(y, m+k, n+l, M, N);
                        minAux = fAux;
                    }
                }
            }
            u[n*M + m] = minY;
        }
    }
    
    // delete arrays
    delete [] G;
    delete [] Z;
}

/*
 * Circle median filter for quantized data y
 * 
 * y: data in M x N, (all elements must be in the interval (-pi, pi])
 * u: output array M x N,
 * M, N: size of data
 * RD, TD: size of filter mask (must be odd integers)
 * v: possible values of y
 * S: number of elements of V
 *
 */
void medfiltCirc2DQuant(double* u, double* y, int M, int N, int RD, int TD, double* v, int S) { 
    // check dimensions
    checkDims(M, N, RD, TD); 

    
    // init
    int R = (RD-1)/2;
    int T = (TD-1)/2;
    
    // create auxiliary values and arrays
    double dAux, vAux, minAux, minY;
    double* F = new double[N * S];
    double* Z = new double[S];
    
    // compute F^00
    for (int s = 0; s < S; s++) {
        dAux = 0;
        vAux = v[s];
        for (int r = -R; r <= R; r++) {
            for (int t = -T; t <= T; t++) {
                dAux += dist(vAux, getMB(y, r, t, M, N));
            }
        }
        F[s] = dAux;
    }
    
    // compute F^0n
    for (int n = 1; n < N; n++) {
        for (int s = 0; s < S; s++) {
            dAux = 0;
            vAux = v[s];
            for (int r = -R; r <= R; r++) {
                dAux +=   dist(vAux, getMB(y, r, n+T, M, N))
                - dist(vAux, getMB(y, r, n-T-1, M, N));
            }
            F[n * S + s] = dAux + F[(n-1) * S + s];
        }
    }
    
    // compute circle medians for first row
    for (int n = 0; n < N; n++) {
        minAux = F[n * S];
        minY = v[0];
        for (int s = 1; s < S; s++) {
            if (F[n * S + s] < minAux) {
                minY = v[s];
                minAux = F[n * S + s];
            }
        }
        u[n*M] = minY;
    }
    
    // process other rows
    for (int m = 1; m < M; m++) {
        // n = 0
        for (int s = 0; s < S; s++) {
            dAux = 0;
            vAux = v[s];
            for (int t = -T; t <= T; t++) {
                dAux +=   dist(vAux, getMB(y, m+R, t, M, N))
                        - dist(vAux, getMB(y, m-R-1, t, M, N));
            }
            Z[s] = dAux;
            F[s] = F[s] + Z[s];
        }
        
        // n = 1,...N-1
        for (int n = 1; n < N; n++) {
            for (int s = 0; s < S; s++) {
                vAux = v[s];
                Z[s] +=   dist(vAux, getMB(y, m+R, n+T, M, N))
                        + dist(vAux, getMB(y, m-R-1, n-T-1, M, N))
                        - dist(vAux, getMB(y, m+R, n-T-1, M, N))
                        - dist(vAux, getMB(y, m-R-1, n+T, M, N));
                
            }
            for (int s = 0; s < S; s++) {
                F[n*S + s] = F[n*S + s] + Z[s];
            }
        }
        
        // compute circle medians
        for (int n = 0; n < N; n++) {
            minAux = F[n*S];
            minY = v[0];
            for (int s = 1; s < S; s++) {
                if (F[n * S + s] < minAux) {
                    minY = v[s];
                    minAux = F[n * S + s];
                }
            }
            u[n*M + m] = minY;
        }
    }
    
    // delete arrays
    delete [] F;
    delete [] Z;
}


/*
 * Geometric median filter in R^K using Weiszfeld algorithm
 * y: data in M x N x K,
 * u: output array M x N x K,
 * M, N, K: size of data
 * RD, TD: size of filter mask (must be odd integers)
 *
 * Warning: This implementation assumes that data is given on the unit circle.
 * Application on general vector-valued data may give NaNs 
 * when the data contains zero vectors
 *
 */
void medfiltGeoRN2D(double* u, double* y, int M, int N, int K, int RD, int TD, int maxIter, double stopTol) {
    // check dimensions
    checkDims(M, N, RD, TD); 

    // init
    int R = (RD-1)/2;
    int T = (TD-1)/2;
    double* med = new double[K];
    double* medOld = new double[K];
    double d, den;
    double relChangeDen, relChangeNum;
    double* num = new double[K];
    int rMB, tMB;
    
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            // init with 0
            for (int k = 0; k < K; k++) {
                med[k] = 0.0;
            }
            // Weiszfeld iteration
            for (int i = 0; i < maxIter; i++) {
                // init numerator and denominator
                for (int k = 0; k < K; k++) {
                    num[k] = 0.0;
                }
                den = 0.0;
                // sum over all values in mask
                for (int r = -R; r <= R; r++) {
                    for (int t = -T; t <= T; t++) {
                        // mirror boundary conditions
                        rMB = getMBInd(m+r, M);
                        tMB = getMBInd(n+t, N);
                        // 2-norm of med - y
                        d = 0.0;
                        for (int k = 0; k < K; k++) {
                            d += pow(med[k] - y[M*N*k + M*tMB  + rMB], 2.0);
                        }
                        d = sqrt(d);
                        // accumulate numerator and denominator
                        den += 1.0/d;
                        for (int k = 0; k < K; k++) {
                            num[k] += y[M*N*k + M*tMB  + rMB] / d;
                        }
                        
                    }
                }
                // break iteration if denominator is not finite
                if (!isfinite(den)) { 
                        break;
                }
                // new iterate
                for (int k = 0; k < K; k++) {
                     medOld[k] = med[k];
                     med[k] = num[k]/den;
                }
                
                // stop criterion (see Fritz et al., Comp. Stat., 2012), 
                // (default stopTol should be 10^-9)
                relChangeNum = 0;
                relChangeDen = 0;
                for (int k = 0; k < K; k++) {
                    relChangeNum += fabs(medOld[k] - med[k]);
                    relChangeDen += fabs(med[k]);
                }
                if (relChangeNum <= (stopTol * relChangeDen)) {
                    break;
                }
            }
            
            // write median to output array u
            for (int k = 0; k < K; k++) {
                    u[M*N*k + M*n + m] = med[k];
            }
        }
    }
    
    // delete arrays
    delete [] med;
    delete [] medOld;
    delete [] num;
    
}
